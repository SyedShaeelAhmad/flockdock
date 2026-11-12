from dataclasses import dataclass, asdict, field
from pathlib import Path
import json
import os
import tempfile

from .config import APP_VERSION, PROJECT_SCHEMA_VERSION

@dataclass
class DockingSettings:
    search_mode: str = 'fast'
    exhaustiveness: int = 8
    num_modes: int = 9
    pocket_mode: str = 'automatic'
    drugs_only: bool = False
    pocket_pad: float = 5.0
    min_box: float = 18.0
    max_box: float = 30.0

@dataclass
class Project:
    name: str
    project_dir: str
    ligand_file: str
    proteins: list[str] = field(default_factory=list)
    settings: DockingSettings = field(default_factory=DockingSettings)
    pocket_configs: dict = field(default_factory=dict)
    uniprot_ids: dict = field(default_factory=dict)

    def save(self, path):
        """Atomically write a versioned project file to avoid partial JSON on a crash."""
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = asdict(self)
        payload["version"] = APP_VERSION
        payload["project_schema_version"] = PROJECT_SCHEMA_VERSION
        fd, tmp_name = tempfile.mkstemp(prefix=f".{target.name}_", suffix=".tmp", dir=target.parent)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2, allow_nan=False)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(tmp_name, target)
        finally:
            if os.path.exists(tmp_name):
                try:
                    os.remove(tmp_name)
                except OSError:
                    pass
