# Single-cell RNA-Seq Differential Expression Analysis Pipeline

This Snakemake pipeline performs differential expression analysis on single-cell RNA-Seq data. It includes steps for quality control, data preprocessing, and generating key plots such as PCA and volcano plots. 

## Pipeline Overview

1. **Quality Control and Data Preprocessing**: Each sample undergoes quality control, and a Seurat object is created for further analysis.
2. **Differential Expression Analysis**: The pipeline calculates differentially expressed genes between specified conditions and generates output files for analysis and visualization, including a volcano plot and a PCA plot.

### Input Data

The pipeline is designed for single-cell RNA-Seq data available from GEO. Specifically, we used data from the following dataset: [GSE132771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771). 

**Dataset Description**:
- **Organism**: human samples
- **Conditions**: 
  - For human samples: Data includes sorted cells from patients with idiopathic pulmonary fibrosis (IPF) and scleroderma, as well as normal controls. 
- **Samples for Testing**: I used data from two scleroderma lungs and two normal lungs.

## Installation

To run this pipeline, you'll need:
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [R](https://www.r-project.org/) with the following packages:
  - `Seurat`
  - `dplyr`
  - `ggplot2`

Make sure R and the necessary packages are installed before running the pipeline.

## Installation and Environment Setup

To run this pipeline, you need to set up a Conda environment with the required dependencies. Follow the steps below to install and activate the environment.

### 1. Install Conda

If you don't already have Conda installed, download and install either [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html). Follow the installation instructions for your operating system.

### 2. Create a Conda Environment

Once Conda is installed, navigate to the directory containing the `environment.yml` file and create the environment by running the following command:

```bash
conda env create -f environment.yml


## Pipeline Structure

- `data/`: Folder containing the input data files.
- `results/`: Folder where results such as differential expression results and plots are saved.
- `scripts/`: Folder containing R scripts for quality control and differential expression analysis.

## Usage

To run the pipeline, navigate to the directory containing the `Snakefile` and execute:

```bash
snakemake --latency-wait 100
