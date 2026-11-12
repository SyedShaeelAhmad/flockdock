import sys
import subprocess
import re
from pathlib import Path

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QFileDialog, QFormLayout,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QListWidget,
    QMainWindow, QMessageBox, QPlainTextEdit, QProgressBar, QScrollArea,
    QPushButton, QSpinBox, QDoubleSpinBox, QSplitter, QVBoxLayout, QWidget
)

from .config import APP_NAME, APP_VERSION, PROJECT_FILE
from .project import Project, DockingSettings
from .worker import PipelineWorker


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.worker = None
        self._library_compounds = 0
        self.setWindowTitle(f"{APP_NAME} {APP_VERSION}")
        self.resize(1250, 800)
        self.build()
        self.gpu()

    def build(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)

        title = QLabel("FlockDOCK")
        title.setObjectName("title")
        root.addWidget(title)
        root.addWidget(QLabel("High-throughput structure-based virtual screening"))

        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter, 1)

        # ---------------- Left panel ----------------
        left = QWidget()
        layout = QVBoxLayout(left)

        project_box = QGroupBox("Project")
        form = QFormLayout(project_box)
        self.name = QLineEdit("My FlockDOCK Project")
        form.addRow("Name", self.name)

        project_row = QHBoxLayout()
        self.pdir = QLineEdit()
        browse = QPushButton("Browse")
        browse.clicked.connect(self.browse_project)
        project_row.addWidget(self.pdir)
        project_row.addWidget(browse)
        form.addRow("Directory", project_row)
        layout.addWidget(project_box)

        library_box = QGroupBox("Compound library")
        library_row = QVBoxLayout(library_box)

        file_row = QHBoxLayout()
        self.lig = QLineEdit()
        self.lig.setPlaceholderText("CSV containing compound ID + SMILES")
        browse_library = QPushButton("Browse")
        browse_library.clicked.connect(self.browse_lig)
        file_row.addWidget(self.lig)
        file_row.addWidget(browse_library)
        library_row.addLayout(file_row)

        self.library_info = QLabel("No library loaded")
        library_row.addWidget(self.library_info)
        layout.addWidget(library_box)

        target_box = QGroupBox("Target proteins")
        target_layout = QVBoxLayout(target_box)

        self.pin = QLineEdit()
        self.pin.setPlaceholderText("Gene symbol, e.g. AKT1")
        self.pin.returnPressed.connect(self.add)

        add_button = QPushButton("+ Add target")
        add_button.clicked.connect(self.add)

        import_targets_button = QPushButton("Import protein list (CSV)")
        import_targets_button.clicked.connect(self.browse_target_csv)
        self.target_csv_info = QLabel("Or import multiple targets from a CSV with gene symbols and optional UniProt IDs.")
        self.target_csv_info.setWordWrap(True)
        self.protein_uniprot = {}

        self.target_count = QLabel("0 targets selected")
        self.targets = QListWidget()
        self.targets.currentItemChanged.connect(self.load_pocket_editor)

        # Per-target pocket editor. The selected target is shown on the left
        # and its pocket definition is edited on the right. This matters when
        # several proteins are selected because each protein can have its own
        # coordinates/residue set.
        self.pocket_configs = {}
        pocket_box = QGroupBox("Binding pocket")
        pocket_form = QFormLayout(pocket_box)
        self.pocket_target_label = QLabel("Select a target")
        self.pocket_mode = QComboBox()
        self.pocket_mode.addItems([
            "Automatic (FlockDOCK logic)",
            "Custom coordinates",
            "Custom residues",
        ])
        self.pocket_mode.currentIndexChanged.connect(self.update_pocket_fields)

        self.coord_widget = QWidget()
        coord_form = QFormLayout(self.coord_widget)
        coord_form.setContentsMargins(0, 0, 0, 0)
        self.cx = self._pocket_spin(0.0)
        self.cy = self._pocket_spin(0.0)
        self.cz = self._pocket_spin(0.0)
        self.sx = self._pocket_spin(20.0, minimum=1.0, maximum=100.0)
        self.sy = self._pocket_spin(20.0, minimum=1.0, maximum=100.0)
        self.sz = self._pocket_spin(20.0, minimum=1.0, maximum=100.0)
        coord_form.addRow("Center X", self.cx)
        coord_form.addRow("Center Y", self.cy)
        coord_form.addRow("Center Z", self.cz)
        coord_form.addRow("Box X", self.sx)
        coord_form.addRow("Box Y", self.sy)
        coord_form.addRow("Box Z", self.sz)

        self.residue_widget = QWidget()
        residue_form = QFormLayout(self.residue_widget)
        residue_form.setContentsMargins(0, 0, 0, 0)
        self.pocket_chain = QLineEdit()
        self.pocket_chain.setPlaceholderText("A")
        self.pocket_residues = QLineEdit()
        self.pocket_residues.setPlaceholderText("151,154,157 or A:151,A:154,A:157")
        self.pocket_pad = self._pocket_spin(5.0, minimum=0.0, maximum=30.0)
        residue_form.addRow("Chain", self.pocket_chain)
        residue_form.addRow("Residues", self.pocket_residues)
        residue_form.addRow("Padding (Å)", self.pocket_pad)

        self.pocket_hint = QLabel("Uses FlockDOCK's automatic cocrystal/P2Rank pocket selection.")
        self.pocket_hint.setWordWrap(True)

        pocket_form.addRow("Target", self.pocket_target_label)
        pocket_form.addRow("Mode", self.pocket_mode)
        pocket_form.addRow("", self.pocket_hint)
        pocket_form.addRow("", self.coord_widget)
        pocket_form.addRow("", self.residue_widget)

        # Keep the pocket editor compact while allowing all controls to be
        # reached on smaller screens. The editor itself scrolls vertically.
        pocket_scroll = QScrollArea()
        pocket_scroll.setWidgetResizable(True)
        pocket_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        pocket_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        pocket_scroll.setMinimumWidth(300)
        pocket_scroll.setWidget(pocket_box)

        target_content = QHBoxLayout()
        target_content.addWidget(self.targets, 1)
        target_content.addWidget(pocket_scroll, 1)

        remove_button = QPushButton("Remove selected")
        remove_button.clicked.connect(self.remove)

        target_layout.addWidget(self.pin)
        target_layout.addWidget(add_button)
        target_layout.addWidget(import_targets_button)
        target_layout.addWidget(self.target_csv_info)
        target_layout.addWidget(self.target_count)
        target_layout.addLayout(target_content, 1)
        target_layout.addWidget(remove_button)
        layout.addWidget(target_box, 1)

        self.update_pocket_fields()

        docking_box = QGroupBox("Docking")
        docking_form = QFormLayout(docking_box)

        self.search = QComboBox()
        self.search.addItems(["fast", "balance", "detail"])

        self.ex = QSpinBox()
        self.ex.setRange(1, 128)
        self.ex.setValue(8)

        self.modes = QSpinBox()
        self.modes.setRange(1, 20)
        self.modes.setValue(9)

        self.drugs = QCheckBox("Drug-like ligands only")

        docking_form.addRow("Search mode", self.search)
        docking_form.addRow("Exhaustiveness", self.ex)
        docking_form.addRow("Binding modes", self.modes)
        docking_form.addRow("", self.drugs)
        layout.addWidget(docking_box)

        gpu_box = QGroupBox("GPU")
        gpu_layout = QVBoxLayout(gpu_box)
        self.gpu_label = QLabel("Checking...")
        gpu_layout.addWidget(self.gpu_label)
        layout.addWidget(gpu_box)

        splitter.addWidget(left)

        # ---------------- Right panel ----------------
        right = QWidget()
        right_layout = QVBoxLayout(right)

        pipeline_box = QGroupBox("Pipeline")
        pipeline_layout = QVBoxLayout(pipeline_box)

        self.stage = QLabel("Ready")
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.hide()

        button_row = QHBoxLayout()
        self.start = QPushButton("START DOCKING")
        self.start.setObjectName("start")
        self.start.clicked.connect(self.start_run)

        self.stop = QPushButton("STOP")
        self.stop.setEnabled(False)
        self.stop.clicked.connect(self.stop_run)

        button_row.addWidget(self.start)
        button_row.addWidget(self.stop)

        pipeline_layout.addWidget(self.stage)
        pipeline_layout.addWidget(self.progress)
        pipeline_layout.addLayout(button_row)
        right_layout.addWidget(pipeline_box)

        log_box = QGroupBox("FlockDOCK log")
        log_layout = QVBoxLayout(log_box)
        self.out = QPlainTextEdit()
        self.out.setReadOnly(True)
        log_layout.addWidget(self.out)
        right_layout.addWidget(log_box, 1)

        splitter.addWidget(right)
        splitter.setSizes([450, 800])

    def gpu(self):
        try:
            result = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=name,memory.total",
                    "--format=csv,noheader",
                ],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0 and result.stdout.strip():
                gpu = result.stdout.strip().splitlines()[0]
                self.gpu_label.setText("✓ " + gpu)
            else:
                self.gpu_label.setText("⚠ NVIDIA GPU not detected")
        except Exception:
            self.gpu_label.setText("⚠ NVIDIA GPU not detected")

    def browse_project(self):
        path = QFileDialog.getExistingDirectory(
            self, "Choose project directory"
        )
        if path:
            self.pdir.setText(path)

    def browse_lig(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Choose compound library",
            "",
            "CSV (*.csv);;All files (*)",
        )
        if path:
            self.lig.setText(path)
            self.load_library(path)

    def load_library(self, path):
        """Validate that the library contains compounds.

        The Gene/Target column is deliberately optional and is never used
        to decide which compounds are docked.
        """
        try:
            df = pd.read_csv(path)
            df.columns = [str(c).strip() for c in df.columns]

            aliases = {
                "id": [
                    "LOTUS_ID",
                    "Ligand_LOTUS_ID",
                    "ligand_id",
                    "compound_id",
                    "molecule_id",
                    "id",
                ],
                "smiles": [
                    "Canonical_SMILES",
                    "Ligand_SMILES",
                    "SMILES",
                    "smiles",
                    "canonical smiles",
                ],
            }

            lower = {c.lower(): c for c in df.columns}

            def resolve(names):
                for name in names:
                    if name.lower() in lower:
                        return lower[name.lower()]
                for name in names:
                    for column in df.columns:
                        if name.lower() in column.lower():
                            return column
                return None

            id_col = resolve(aliases["id"])
            smiles_col = resolve(aliases["smiles"])

            if id_col is None or smiles_col is None:
                self.library_info.setText(
                    "⚠ Missing compound ID or SMILES column"
                )
                QMessageBox.warning(
                    self,
                    "Invalid compound library",
                    "FlockDOCK needs a compound ID column and a SMILES column.\n\n"
                    f"Detected columns:\n{', '.join(df.columns)}",
                )
                return False

            valid = df[id_col].notna() & df[smiles_col].notna()
            n_valid = int(valid.sum())

            if n_valid == 0:
                self.library_info.setText("⚠ No valid compounds found")
                QMessageBox.warning(
                    self,
                    "Empty compound library",
                    "No rows containing both a compound ID and SMILES were found.",
                )
                return False

            self._library_compounds = n_valid
            self.library_info.setText(
                f"✓ {n_valid:,} compounds ready • "
                "target list is independent of library"
            )
            return True

        except Exception as exc:
            self.library_info.setText("⚠ Could not read library")
            QMessageBox.warning(
                self,
                "Library error",
                f"Could not read the CSV:\n\n{exc}",
            )
            return False

    def browse_target_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Import protein target list", "", "CSV files (*.csv);;All files (*)"
        )
        if path:
            self.load_target_csv(path)

    def load_target_csv(self, path):
        """Import and de-duplicate protein targets from a CSV."""
        try:
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
            frame.columns = [str(c).strip() for c in frame.columns]
            norm = {c.lower().replace("_", "").replace(" ", ""): c for c in frame.columns}

            def find_col(candidates):
                for candidate in candidates:
                    key = candidate.lower().replace("_", "").replace(" ", "")
                    if key in norm:
                        return norm[key]
                return None

            gene_col = find_col(["HGNCSymbol", "HGNC", "GeneSymbol", "Gene", "Symbol", "Target"])
            uniprot_col = find_col(["UniProt", "UniProtID", "UniProtAccession", "Accession", "ProteinID"])
            if gene_col is None:
                raise ValueError("Could not identify a gene-symbol column. Include HGNCSymbol, Gene, GeneSymbol or Symbol.")

            imported = {}
            for idx, row in frame.iterrows():
                gene = str(row[gene_col]).strip().upper()
                accession = str(row[uniprot_col]).strip().upper() if uniprot_col else ""
                if not gene or gene.lower() in {"nan", "none", "null"}:
                    continue
                if accession and not re.fullmatch(r"[A-Z0-9]{6,10}(?:-[0-9]+)?", accession):
                    raise ValueError(f"Invalid UniProt accession '{accession}' on CSV row {idx + 2}.")
                if gene in imported and imported[gene] and accession and imported[gene] != accession:
                    raise ValueError(f"Conflicting UniProt IDs for {gene}: {imported[gene]} and {accession}. Resolve duplicate rows first.")
                if gene not in imported or (not imported[gene] and accession):
                    imported[gene] = accession

            if not imported:
                raise ValueError("No valid gene symbols were found in the CSV.")

            added = 0
            existing = {self.targets.item(i).text().split("  [UniProt:")[0].upper() for i in range(self.targets.count())}
            for gene, accession in imported.items():
                if gene not in existing:
                    self.targets.addItem(gene)
                    added += 1
                if accession:
                    self.protein_uniprot[gene] = accession
                item = next((self.targets.item(i) for i in range(self.targets.count())
                             if self.targets.item(i).text().split("  [UniProt:")[0].upper() == gene), None)
                if item and accession:
                    item.setToolTip(f"UniProt: {accession}")
                    item.setText(f"{gene}  [UniProt: {accession}]")

            self.target_csv_info.setText(
                f"✓ Imported {len(imported)} unique proteins ({added} newly added); "
                f"{sum(bool(v) for v in imported.values())} UniProt IDs mapped."
            )
            self.update_target_count()
            if self.targets.count():
                self.targets.setCurrentRow(0)
            QMessageBox.information(
                self, "Protein list imported",
                f"Detected {len(imported)} unique proteins.\nAdded {added} new targets.\n"
                f"UniProt IDs mapped: {sum(bool(v) for v in imported.values())}.\n\n"
                "All imported targets will be processed in the docking pipeline."
            )
            return True
        except Exception as exc:
            QMessageBox.warning(self, "Protein CSV error", str(exc))
            self.target_csv_info.setText("⚠ Protein CSV could not be imported.")
            return False

    def add(self):
        target = self.pin.text().strip().upper()
        if not target:
            return

        existing = [
            self.targets.item(i).text().split("  [UniProt:")[0].upper()
            for i in range(self.targets.count())
        ]

        if target not in existing:
            self.targets.addItem(target)
            self.pocket_configs[target] = {"mode": "automatic"}
            self.targets.setCurrentRow(self.targets.count() - 1)

        self.pin.clear()
        self.update_target_count()

    def update_target_count(self):
        count = self.targets.count()
        suffix = "" if count == 1 else "s"
        self.target_count.setText(f"{count} target{suffix} selected")

    def remove(self):
        for item in self.targets.selectedItems():
            target = item.text().split("  [UniProt:")[0].strip().upper()
            self.pocket_configs.pop(target, None)
            self.targets.takeItem(self.targets.row(item))
        self.update_target_count()
        self.load_pocket_editor(self.targets.currentItem(), None)

    def _pocket_spin(self, value=0.0, minimum=-1000.0, maximum=1000.0):
        box = QDoubleSpinBox()
        box.setRange(minimum, maximum)
        box.setDecimals(3)
        box.setSingleStep(0.1)
        box.setValue(value)
        return box

    def _current_target(self):
        item = self.targets.currentItem()
        return item.text().split("  [UniProt:")[0].strip().upper() if item else None

    def _read_pocket_editor(self):
        target = self._current_target()
        if not target:
            return
        mode_index = self.pocket_mode.currentIndex()
        if mode_index == 0:
            cfg = {"mode": "automatic"}
        elif mode_index == 1:
            cfg = {
                "mode": "coordinates",
                "center": [self.cx.value(), self.cy.value(), self.cz.value()],
                "size": [self.sx.value(), self.sy.value(), self.sz.value()],
            }
        else:
            cfg = {
                "mode": "residues",
                "chain": self.pocket_chain.text().strip(),
                "residues": self.pocket_residues.text().strip(),
                "pad": self.pocket_pad.value(),
            }
        self.pocket_configs[target] = cfg

    def load_pocket_editor(self, current, previous=None):
        # Save the previous target before loading the newly selected one.
        if previous is not None:
            prev_target = previous.text().split("  [UniProt:")[0].strip().upper()
            self._read_pocket_editor_for(prev_target)

        target = current.text().split("  [UniProt:")[0].strip().upper() if current else None
        if not target:
            self.pocket_target_label.setText("Select a target")
            self.pocket_mode.blockSignals(True)
            self.pocket_mode.setCurrentIndex(0)
            self.pocket_mode.blockSignals(False)
            self.update_pocket_fields()
            return

        cfg = self.pocket_configs.get(target, {"mode": "automatic"})
        self.pocket_target_label.setText(target)
        self.pocket_mode.blockSignals(True)
        mode = cfg.get("mode", "automatic")
        self.pocket_mode.setCurrentIndex({
            "automatic": 0, "coordinates": 1, "residues": 2
        }.get(mode, 0))
        self.pocket_mode.blockSignals(False)

        center = cfg.get("center", [0.0, 0.0, 0.0])
        size = cfg.get("size", [20.0, 20.0, 20.0])
        for spin, value in zip((self.cx, self.cy, self.cz), center):
            spin.setValue(float(value))
        for spin, value in zip((self.sx, self.sy, self.sz), size):
            spin.setValue(float(value))
        self.pocket_chain.setText(str(cfg.get("chain", "")))
        self.pocket_residues.setText(str(cfg.get("residues", "")))
        self.pocket_pad.setValue(float(cfg.get("pad", 5.0)))
        self.update_pocket_fields()

    def _read_pocket_editor_for(self, target):
        if not target:
            return
        mode_index = self.pocket_mode.currentIndex()
        if mode_index == 0:
            cfg = {"mode": "automatic"}
        elif mode_index == 1:
            cfg = {"mode": "coordinates",
                   "center": [self.cx.value(), self.cy.value(), self.cz.value()],
                   "size": [self.sx.value(), self.sy.value(), self.sz.value()]}
        else:
            cfg = {"mode": "residues",
                   "chain": self.pocket_chain.text().strip(),
                   "residues": self.pocket_residues.text().strip(),
                   "pad": self.pocket_pad.value()}
        self.pocket_configs[target] = cfg

    def update_pocket_fields(self):
        mode = self.pocket_mode.currentIndex()
        self.coord_widget.setVisible(mode == 1)
        self.residue_widget.setVisible(mode == 2)
        if mode == 0:
            self.pocket_hint.setText("Uses FlockDOCK's automatic cocrystal-ligand logic, with P2Rank as fallback.")
        elif mode == 1:
            self.pocket_hint.setText("Enter the Uni-Dock box center and dimensions in Å.")
        else:
            self.pocket_hint.setText("Enter residues that define the pocket. The box is centered on those residues and expanded by the padding value.")
        if self._current_target():
            self._read_pocket_editor()

    def validate_pocket_configs(self):
        # Commit the currently displayed target before validation.
        self._read_pocket_editor()
        for i in range(self.targets.count()):
            target = self.targets.item(i).text().split("  [UniProt:")[0].strip().upper()
            cfg = self.pocket_configs.get(target, {"mode": "automatic"})
            if cfg.get("mode") == "residues" and not cfg.get("residues", "").strip():
                QMessageBox.warning(self, "Missing pocket residues",
                                     f"Enter at least one residue for {target}, or choose Automatic/Custom coordinates.")
                self.targets.setCurrentRow(i)
                return False
            if cfg.get("mode") == "residues" and not cfg.get("chain", "").strip():
                QMessageBox.warning(self, "Missing chain",
                                     f"Enter the protein chain for the selected residues of {target}.")
                self.targets.setCurrentRow(i)
                return False
        return True

    def start_run(self):
        # Automatically add a typed target if the user forgot to click Add.
        if self.pin.text().strip():
            self.add()

        project_text = self.pdir.text().strip()
        ligand_text = self.lig.text().strip()

        if not project_text:
            QMessageBox.warning(
                self,
                "Missing project directory",
                "Choose a project directory first.",
            )
            return

        if not ligand_text:
            QMessageBox.warning(
                self,
                "Missing compound library",
                "Choose your compound CSV first.",
            )
            return

        if not Path(ligand_text).is_file():
            QMessageBox.warning(
                self,
                "Compound library not found",
                "The selected CSV does not exist on this Windows computer:\n\n"
                f"{ligand_text}\n\n"
                "Click Browse and select the actual compound library.",
            )
            return

        if not self.load_library(ligand_text):
            return

        if self.targets.count() == 0:
            QMessageBox.warning(
                self,
                "Missing target",
                "Add at least one target protein.",
            )
            return

        project_dir = Path(project_text).expanduser().resolve()
        try:
            project_dir.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            QMessageBox.critical(self, "Project directory error", f"Could not create or access the project directory:\n\n{project_dir}\n\n{exc}")
            return

        targets = [
            self.targets.item(i).text().split("  [UniProt:")[0].strip().upper()
            for i in range(self.targets.count())
        ]

        if not self.validate_pocket_configs():
            return

        settings = DockingSettings(
            self.search.currentText(),
            self.ex.value(),
            self.modes.value(),
            "automatic",
            self.drugs.isChecked(),
        )

        project = Project(
            name=self.name.text().strip() or "FlockDOCK Project",
            project_dir=str(project_dir),
            ligand_file=self.lig.text().strip(),
            proteins=targets,
            settings=settings,
            pocket_configs={t: self.pocket_configs.get(t, {"mode": "automatic"}) for t in targets},
            uniprot_ids={t: self.protein_uniprot.get(t, "") for t in targets if self.protein_uniprot.get(t, "")},
        )

        try:
            project.save(project_dir / PROJECT_FILE)
        except Exception as exc:
            QMessageBox.critical(self, "Could not save project", f"FlockDOCK could not safely write the project file:\n\n{exc}")
            return

        self.out.clear()
        self.out.appendPlainText("Starting FlockDOCK...")

        base_dir = Path(
            getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[1])
        )
        scripts = base_dir / "flockdock_scripts"

        self.worker = PipelineWorker(project, scripts, self)
        self.worker.output.connect(self.out.appendPlainText)
        self.worker.stage.connect(self.stage.setText)
        self.worker.done.connect(self.finished)
        self.worker.failed.connect(self.failed)
        self.worker.start()

        self.start.setEnabled(False)
        self.stop.setEnabled(True)
        self.progress.show()
        self.stage.setText("Starting")

    def stop_run(self):
        if self.worker:
            self.worker.stop()
            self.stage.setText("Stopping...")

    def finished(self, code, elapsed):
        self.start.setEnabled(True)
        self.stop.setEnabled(False)
        self.progress.hide()

        if code == 0:
            self.stage.setText(f"Completed in {elapsed / 60:.1f} min")
        else:
            self.stage.setText(f"Failed (exit {code})")

    def failed(self, error):
        self.start.setEnabled(True)
        self.stop.setEnabled(False)
        self.progress.hide()
        self.stage.setText("Error")
        self.out.appendPlainText("ERROR: " + error)


def main():
    app = QApplication(sys.argv)
    app.setStyleSheet(
        "#title { font-size: 30px; font-weight: 700; } "
        "QPushButton#start { font-weight: 700; padding: 10px; }"
    )
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
