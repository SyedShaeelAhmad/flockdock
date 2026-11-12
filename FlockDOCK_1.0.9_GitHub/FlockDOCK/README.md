# FlockDOCK

FlockDOCK is a Windows application for preparing small-molecule libraries, selecting protein binding pockets, and running GPU-accelerated molecular docking. It uses RDKit, Meeko, P2Rank, and UniDock in WSL2.

## Requirements

- Windows 10 or 11, 64-bit
- Python 3.11, 3.12, or 3.13
- WSL2 with Ubuntu
- An NVIDIA GPU and a Windows driver that supports CUDA in WSL2
- Java 17 or later for automatic pocket prediction with P2Rank

The setup script installs the FlockDOCK Python dependencies. On first use, FlockDOCK can install P2Rank and UniDock in their respective environments.

## Install and run

1. Download or clone this repository.
2. Run `setup_FlockDOCK.bat` and wait for dependency installation to finish.
3. Run `FlockDOCK.bat`.
4. Choose a compound library CSV, a project folder, and one or more protein targets.
5. Select automatic pocket detection or provide a pocket definition, then start docking.

The compound CSV must include a compound identifier and a SMILES string for each molecule. Target lists can be entered in the app or imported from a CSV with gene symbols and optional UniProt accessions.

## Pocket selection

FlockDOCK can use a co-crystallised ligand, a P2Rank prediction, user-specified coordinates, or selected protein residues to define a docking box. Some structures do not provide a usable pocket. In that case, select a suitable structure or define the pocket manually; the software does not invent a binding site.

## Results

Each project folder contains the input library and generated ligand, protein, pocket, docking, and log files. The `07_results` folder contains docking poses and scores. The `08_logs` folder records stage output and target status. Keep these files with the project when documenting or reproducing an analysis.

## Scientific use

Docking scores are approximate ranking measures. They are not measured binding affinities and do not establish biological activity. Inspect structures and poses, use appropriate controls, and validate important findings experimentally. Predicted structures and pockets require particular review.

## License

FlockDOCK is distributed under the MIT License. See [LICENSE](LICENSE).
