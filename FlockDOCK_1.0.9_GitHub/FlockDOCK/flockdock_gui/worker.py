import os, subprocess, time, sys
from PySide6.QtCore import QThread, Signal

class PipelineWorker(QThread):
    output = Signal(str)
    stage = Signal(str)
    done = Signal(int, float)
    failed = Signal(str)
    def __init__(self, project, scripts_dir, parent=None):
        super().__init__(parent)
        self.project, self.scripts_dir, self.proc = project, scripts_dir, None
        self._stop_requested = False
    def run(self):
        start = time.time()
        env = os.environ.copy()
        env.update({"FLOCKDOCK_PROJECT": self.project.project_dir, "UNIDOCK_PROJECT": self.project.project_dir,
                    "GENDOCK_PROJECT": self.project.project_dir,
                    "ONEDOCK_PROJECT": self.project.project_dir, "PYTHON": sys.executable})
        try:
            cmd = [sys.executable, "-u", str(self.scripts_dir / "flockdock_pipeline.py")]
            self.output.emit(f"Python environment: {sys.executable}")
            self.proc = subprocess.Popen(cmd, cwd=self.project.project_dir, env=env, stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                                         bufsize=1, creationflags=getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0))
            assert self.proc.stdout is not None
            for line in self.proc.stdout:
                line = line.rstrip("\r\n")
                self.output.emit(line)
                if line.startswith("STAGE:"): self.stage.emit(line[6:].strip())
                if self._stop_requested: self._terminate_tree(); break
            code = self.proc.wait()
            if self._stop_requested: self.output.emit("Run stopped by user.")
            self.done.emit(code if not self._stop_requested else 130, time.time() - start)
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}")
    def _terminate_tree(self):
        if not self.proc or self.proc.poll() is not None: return
        try:
            if os.name == "nt": subprocess.run(["taskkill", "/PID", str(self.proc.pid), "/T", "/F"],
                                                 capture_output=True, timeout=10)
            else: self.proc.terminate()
        except Exception: pass
    def stop(self):
        self._stop_requested = True
        self._terminate_tree()
