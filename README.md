# methylscaper

`methylscaper` is an R package for visualizing data that jointly profile endogenous methylation and chromatin accessibility (MAPit, NOMe-seq, scNMT-seq, nanoNOMe, etc.). The package offers pre-processing for single-molecule data and accepts input from Bismark (or similar alignment programs) for single-cell data. A common interface for visualizing both data types is constructed by generating ordered representational methylation-state matrices. The package provides a Shiny app to allow for interactive and optimal ordering of the individual DNA molecules to discover methylation patterns, nucleosome positioning, and transcription factor binding.

# Webserver

`methylscaper` is available for use on a webserver located at [www.methylscaper.com](http://www.methylscaper.com) and clicking on "Start Shiny App". Also located on the website is easy access to the vignette and example datasets to use.

# Installation

For local use of methylscaper, it can be installed into R from [Bioconductor](http://bioconductor.org/packages/devel/bioc/html/methylscaper.html) (using R version > 4.1.0): 

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("methylscaper")
```

Alternatively, the version specified by the `ref` parameter below only requires R >= 4.0 (current stable release).

```{r}
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("rhondabacher/methylscaper", ref="R4.0")
```

## Development Setup

To contribute to this package, follow these steps after cloning the repository:

Run this command in bash
`git config core.hooksPath .hooks`

Verify if the hooks are installed
`ls -la .hooks`

If you encounter any issues with the pre-commit hooks
- Ensure you've run the `git config` command as described above.
- Check that the .hooks/pre-commit file is executable. If not, run:
`chmod +x .hooks/pre-commit`

# Dependencies

Note: on Ubuntu, users may need to install libgsl via:
`sudo apt-get install libgsl-dev`

The following packages are required for methylscaper. If installation fails, you may need to manually install the dependencies using the function 'install.packages' for CRAN packages or 'BiocManager::install' for Bioconductor packages.

* CRAN packages: shiny, shinyjs, shinyFiles, seriation, seqinr, Rfast, data.table
* Bioconductor packages: BiocParallel, Biostrings, SummarizedExperiment


For any other installation issues/questions please leave a message on our [Github Issues](https://github.com/rhondabacher/methylscaper/issues).

# Manuscript

The methylscaper manuscript is now available at [Bioinformatics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab438/6298588). 

# Running methylscaper

To load the package into R:

```{r}
library(methylscaper)
```

To use methylscaper as a Shiny App:
```{r}
methylscaper::methylscaper()
```

Many more details and examples on how to use the Shiny App or functions directly in R are located in the vignette:

* [Vignette](http://bioconductor.org/packages/devel/bioc/vignettes/methylscaper/inst/doc/methylScaper.html)


# Contact/Maintainer

* [Rhonda Bacher](https://www.rhondabacher.com) (rbacher@ufl.edu)  
Department of Biostatistics, University of Florida
