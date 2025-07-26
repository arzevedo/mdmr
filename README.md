# 📦 MDMr: Bayesian Network Modeling for Dynamic Multivariate Time Series

**MDMr** is an R package for learning the structure and estimating the dynamic parameters of Bayesian networks from multivariate time series. It integrates structure learning algorithms — including `bnlearn::hc` and `GOBNILP` — with Kalman filtering and smoothing to estimate time-varying parameters for each node.

---

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/arzevedo/mdmr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/arzevedo/mdmr/actions)

## 🚀 Installation

```r
# Install dependencies (if not already installed)
install.packages(c("bnlearn", "ggplot2", "ggstream", "magick", "Rcpp", "reticulate"))

# Then install MDMr from source
devtools::install_github("arzevedo/mdmr")
```

---

## 📘 Example Usage

This walkthrough demonstrates how to use `MDMr` to learn a Bayesian network and visualize the results.

### 1. Load the package and sample data

```r
library(mdmr)

# Example multivariate time series data [T x N]
head(y)
```

### 2. Fit the MDM model

```r
# Estimate structure and dynamics with default settings
res <- mdm(y)
```

- `res` is an object of class `"mdm"` containing:
  - The inferred DAG structure
  - Filtering and smoothing estimates
  - Local scores and optimization metadata

### 3. Print summary of results

```r
srt(res)
```

---

## 📊 Visualizations

### DAG Structure

```r
# Plot the estimated DAG using bnlearn::graphviz.plot
plot_dag(res)
```

> Displays the structure as a directed acyclic graph with nodes and directed edges.

---

### Arcs Over Time

```r
# Highlight which arcs were selected and their local scores
plot_arcs(res)
```

> Helps identify dominant parent-child relationships and score-driven decisions.

---

### Dynamic Heatmap Animation

```r
# Animated heatmap of posterior estimates across time
plot_heatmap_animation(res)
```

> Each tile represents the magnitude of a dynamic parameter at each time step. This animation captures temporal variation in the network structure.

---

### Streamplot for a Node

```r
# Plot contribution of parents to the dynamic evolution of a target node
plot_stream(res, child_node = 2)
```

> Useful for assessing how different parents dynamically influence a given node.

---

### Marginal Posterior for a Node

```r
# Plot the marginal posterior means and confidence bands
plot_marginal(res, target_node = 2)
```

> Shows how the coefficients associated with a node evolve over time (filtered or smoothed).

---

## 📚 Citation

If you use this package in your research, please cite:

```bibtex
@manual{mdmr2025,
  title = {MDMr: Bayesian Dynamic Structure Learning for Multivariate Time Series},
  author = {Arthur R Azevedo},
  year = {2025},
  url = {https://github.com/arzevedo/mdmr}
}
```


---

## 🔧 Structure Learning Backends

By default, the `mdm()` function uses the hill-climbing algorithm from the `bnlearn` package to learn the structure of the Bayesian network.

If you prefer to use the Integer Programming Approach (IPA) via GOBNILP, make sure to:

1. Properly install **SCIP** on your machine.
2. Compile and install the **GOBNILP** binary.
3. Set the binary path in R using:

```r
options(gobnilp.path = "/full/path/to/gobnilp")
```

This enables `mdmr` to use `GOBNILP` for structure optimization via the `run_gobnilp()` interface.


---

## ⚙️ Advanced

To use `GOBNILP` as the structure optimizer, ensure that:
- GOBNILP is compiled and its binary is accessible
- The `gobnilp.path` option is set via `options(gobnilp.path = "path/to/gobnilp")`

The function `run_gobnilp()` handles the interface.
