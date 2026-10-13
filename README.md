# Extreme Protein Detection RKHS

A Python project for detecting extreme proteins in MSA using a reproducing kernel Hilbert space (RKHS) approach.

## Overview

This repository implements an end-to-end workflow for analyzing protein detection data and identifying extreme or anomalous protein patterns using kernel-based methods. It is designed for research and exploratory analysis prior to phylogenetic analysis.

## Features

- Kernel-based protein detection modeling
- RKHS feature representation and similarity analysis
- Experimental pipeline for classification and anomaly detection
- Reproducible project structure for research workflows

## Project Structure

```text
.
├── README.md
├── .gitignore
├── requirements.txt
├── detect_outliers.py
├── src/
│   └── ...
├── notebooks/
│   └── ...
├── example/
│   └── ...
└── docs/
    └── ...

```

## Installation

Clone the repository:

```bash
git clone https://github.com/ycharlie-odt/exteme-protein-detection-rkhs.git
cd exteme-protein-detection-rkhs
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

Run the analysis pipeline or model scripts from the project source directory:

```bash
python3 detect_outliers.py -i example/sample_alignment.fasta -o example/results.csv --kernel id --gamma 1.0 --delta 0.05 --clean-fasta  --verbose
```

For notebook-based exploration:

```bash
jupyter notebook
```

## Requirements

This project typically relies on:

- Python 3.9+
- NumPy
- SciPy
- scikit-learn
- pandas
- matplotlib
- Jupyter

If your project uses different dependencies, update `requirements.txt` accordingly.

## Contributing

Contributions are welcome. Please open an issue or submit a pull request with a brief description of the change and relevant tests.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.

## Acknowledgments

Built for computational biology and machine learning research in protein detection and RKHS-based modeling.
