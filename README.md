# Efficient Learning of Differential Networks in Multi-Source Non-Paranormal Graphical Models

This repository provides the source code and data for reproducing the results in our paper:

> **Efficient Learning of Differential Networks in Multi-Source Non-Paranormal Graphical Models**  
> Mojtaba Nikahd, Seyed Abolfazl Motahari  
> [Manuscript link](https://doi.org/10.48550/arXiv.2410.02496)

---

## 📂 Repository Structure

```
.
├── SPDtrace/
│   ├── Library.R             # Utility functions
│   └── SolutionPath.cpp      # Core implementation of the proposed SPD-trace method
├── Simulation/
│   ├── CrossFDTL.cpp         # C++ implementation of CrossFDTL
│   ├── editedDTrace.R        # Modified version of the D-trace method from DiffGraph
│   ├── First_Scenario.R      # First synthetic experiment: accuracy & runtime comparison
│   └── Second_Scenario.R     # Second synthetic experiment: heterogeneous datasets
└── Overian_cancer/
    ├── Ovarian.R             # Main script for real-data application (ovarian cancer)
    ├── Data/                 # Gene expression data (preprocessed)
    └── Results/              # Saved outputs from the real data analysis
```

---

## 🧪 Experiments Overview

### 🔬 Synthetic Experiments
- **First Scenario:** Compares SPD-trace with APGD D-trace and CrossFDTL on runtime and accuracy.
- **Second Scenario:** Evaluates robustness of SPD-trace using heterogeneous vs. homogeneous datasets.
- Metrics: Precision, recall, and computation time.

### 🧬 Real-World Application: Ovarian Cancer (TCGA)
- Inference of differential networks between platinum-sensitive and platinum-resistant tumors.
- Model selection via StARS.
- Biological validation via enrichment tests using known gene sets.

---

## ⚙️ Installation & Requirements

Before running the code, ensure you have the following:

### 📦 Required Software
- **R** (>= 4.0.0)
- **C++ compiler** (with Rcpp support)

### 📦 Required R Packages

Install the necessary packages using the following in R:

```r
install.packages(c("Rcpp", "igraph", "Matrix"))
# You may also need 'DiffGraph' for reference:
# install.packages("DiffGraph")
```

> Note: `editedDTrace.R` is a modified version of the D-trace method from the `DiffGraph` package.

---

## ▶️ How to Run

Clone the repository and launch R or use terminal commands as shown below.

### 1️⃣ Synthetic Experiments

**First Scenario: Performance Comparison**
```bash
Rscript Simulation/First_Scenario.R
```

**Second Scenario: Heterogeneous Datasets**
```bash
Rscript Simulation/Second_Scenario.R
```

### 2️⃣ Real Data Experiment (Ovarian Cancer)

```bash
Rscript Overian_cancer/Ovarian.R
```

The output including the inferred differential network will be saved to `Overian_cancer/Results/`.

---

## 📖 Citation

If you find this work useful, please cite it using the following BibTeX entry:

```bibtex
@article{Nikahd2025Differential,
  title     = {Efficient Learning of Differential Networks in Multi-Source Non-Paranormal Graphical Models},
  author    = {Nikahd, Mojtaba and Motahari, Seyed Abolfazl},
  year      = {2025},
  note      = {Manuscript under review},
  url       = {https://doi.org/10.48550/arXiv.2410.02496}
}
```

---

## 📫 Contact

For questions, feedback, or bug reports, feel free to:

- Open an issue on GitHub  
- Contact: [nikahd@ce.sharif.ir](mailto:nikahd@ce.sharif.ir)

---

## 📝 License

This project is released for academic and research use only.  
Please cite the paper if you use any part of this code or data.
