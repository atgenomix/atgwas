# ATGWAS: Spark-Powered Shiny App for GWAS Visualization

ATGWAS is an interactive R Shiny application designed for **Genome-Wide Association Study (GWAS)** result exploration. It connects directly to a **Spark-based catalog**, enabling scalable analysis of large GWAS summary statistics.

## 🚀 Features

- 📊 **Interactive Manhattan Plot** using Plotly  
- 📈 **QQ Plot** for visualizing p-value distributions  
- 🧬 **GWAS Table Viewer** with download support  
- 🔍 **Dynamic database selector** based on Spark catalog  
- ⚡ **High-performance Spark backend** with `sparklyr`  
- 🎨 Customizable genome-wide significance threshold  

## 📦 Required R Packages

The app imports the following packages:

```r
shiny, qqman, DT, sparklyr, DBI, shinycssloaders, plotly, scales, dplyr
````

Install them via:

```r
install.packages(c("shiny", "qqman", "DT", "DBI", "shinycssloaders", "plotly", "scales", "dplyr"))
# For sparklyr
install.packages("sparklyr")
```

## 🔧 Usage

```r
# Install from GitHub
remotes::install_github("atgenomix/atgwas")
# Launch the viewer
atgwas::gwasViewer()
```

### 🖥️ UI Overview

* **Sidebar**

  * `Genome-wide threshold`: Customize p-value significance cutoff
  * `Database Selector`: Choose Spark SQL database to analyze

* **Tabs**

  * `Table`: Displays the GWAS result table with CSV export
  * `QQ Plot`: Observed vs expected p-values (supports PNG export)
  * `Manhattan`: Interactive Manhattan plot with zoom, tooltips

## 📂 Data Requirements

The selected Spark database table must contain at least:

* `CHR`: Chromosome identifier
* `BP`: Base-pair position
* `P`: p-value
* (Optional) `SNP`: Variant ID

## 📈 Visualization Details

### Manhattan Plot

Significant variants (e.g., `P < 5e-8`) are highlighted interactively. The app preprocesses and subsets data to optimize rendering for large datasets.

### QQ Plot

Displays p-value inflation/deflation to detect population stratification or biases.

## 🧪 Developer Modules

* `dbBrowserUI()` / `dbBrowserServer()` — Spark database selection
* `sparkConnectionUI()` / `sparkConnectionServer()` — Custom Spark connection
* `plot_manhattan()` — Static ggplot2-based Manhattan plot
* `prep_manhattan()` — Data preprocessing for Manhattan plot

## 🛠 System Requirements

* **R** ≥ 4.1
* **Spark** ≥ 3.0
* Working Spark cluster or local instance
* Internet connection for package installation

## 📤 Export Options

* **CSV Download** — Filtered GWAS result table
* **PNG Download** — QQ Plot

## 👨‍🔬 Example Workflow

1. Start the app and connect to your Spark master
2. Select a database from the sidebar
3. View and download the GWAS summary table
4. Switch to the QQ plot tab and export as PNG
5. Explore the Manhattan plot interactively

---

## 📄 License

This project is licensed under the **Apache License 2.0**.

Copyright 2025 atgenomix, Inc. (Charles Chuang)

See the [LICENSE](./LICENSE) file for details.
## 🤝 Acknowledgements

* [`sparklyr`](https://spark.rstudio.com/)
* [`plotly`](https://plotly.com/r/)
* [`DT`](https://rstudio.github.io/DT/)
* [`shiny`](https://shiny.rstudio.com/)
