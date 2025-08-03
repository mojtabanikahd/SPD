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
  - **Windows:** Install Rtools (https://cran.r-project.org/bin/windows/Rtools/)
  - **macOS:** Install Xcode Command Line Tools: `xcode-select --install`
  - **Linux:** Install build-essential: `sudo apt-get install build-essential`

### 📦 Required R Packages

Install the necessary packages using the following in R:

```r
# Core packages for C++ compilation
install.packages("Rcpp")

# Statistical and data manipulation packages
install.packages(c("stats", "dplyr", "stringr"))

# Visualization packages
install.packages(c("ggplot2", "ggpubr"))

# Network analysis and graphical models
install.packages(c("igraph", "GGMselect"))

# Data I/O packages
install.packages(c("rmatio", "readxl"))

# Additional dependencies (if not already installed)
install.packages(c("Matrix", "methods"))
```

### 📋 Package Dependencies with Version Placeholders

The following R packages are required with their minimum versions:

```r
# Core dependencies
Rcpp >= 1.0.0
stats >= 4.0.0

# Data manipulation and analysis
dplyr >= 1.0.0
stringr >= 1.4.0

# Visualization
ggplot2 >= 3.3.0
ggpubr >= 0.4.0

# Network analysis and graphical models
igraph >= 1.2.0
GGMselect >= 0.1.0

# Data I/O
rmatio >= 0.15.0
readxl >= 1.3.0
```

---

## 🔧 C++ Build Instructions

### Compiling C++ Code

The project contains C++ files that need to be compiled:

1. **SPDtrace/SolutionPath.cpp** - Core SPD-trace implementation
2. **Simulation/CrossFDTL.cpp** - CrossFDTL algorithm implementation

#### Method 1: Using Rcpp::sourceCpp()

Add the following to your R scripts to compile C++ code:

```r
library(Rcpp)

# Compile SolutionPath.cpp
sourceCpp("SPDtrace/SolutionPath.cpp")

# Compile CrossFDTL.cpp
sourceCpp("Simulation/CrossFDTL.cpp")
```

#### Method 2: Using R Package Structure

For more complex projects, you can create a proper R package structure:

1. Create `src/` directory and move C++ files there
2. Create `R/` directory for R functions
3. Use `Rcpp::compileAttributes()` to generate R wrappers:

```r
library(Rcpp)
compileAttributes()
```

#### Method 3: Using Makevars (Advanced)

Create a `~/.R/Makevars` file for custom compilation flags:

```makefile
# For Windows (in ~/.R/Makevars.win)
CXX14 = g++ -std=c++14
CXX14STD = -std=c++14
CXX14FLAGS = -O2
CXX14PICFLAGS = -fpic
CXX14SHLIB = g++ -shared
```

---

## ▶️ How to Run

### Prerequisites Check

First, verify all dependencies are installed:

```r
# Check if all required packages are available
required_packages <- c("Rcpp", "stats", "dplyr", "stringr", "ggplot2", 
                      "ggpubr", "igraph", "GGMselect", "rmatio", "readxl")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if(length(missing_packages) > 0) {
  cat("Missing packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Please install them using: install.packages(c(", 
      paste0('"', missing_packages, '"', collapse = ", "), "))\n")
}
```

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

## 🔍 Troubleshooting

### Common Issues

1. **C++ compilation errors:**
   - Ensure Rtools (Windows) or Xcode (macOS) is properly installed
   - Check that `Rcpp` package is installed and loaded
   - Verify C++ compiler is in your system PATH

2. **Package loading errors:**
   - Update R to latest version
   - Install packages from CRAN: `install.packages("package_name")`
   - For Bioconductor packages: `BiocManager::install("package_name")`

3. **Memory issues:**
   - Increase R memory limit: `memory.limit(size = 8000)`
   - Use garbage collection: `gc()`

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
