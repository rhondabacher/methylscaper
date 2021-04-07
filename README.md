# methylscaper

`methylscaper` is an R package for visualizing data that jointly profile endogenous methylation and chromatin accessibility (MAPit, NOMe-seq, scNMT-seq, nanoNOMe, etc.). The package offers pre-processing for single-molecule data and accepts input from Bismark (or similar alignment programs) for single-cell data. A common interface for visualizing both data types is done by generating ordered representational methylation-state matrices. The package provides a Shiny app to allow for interactive and optimal ordering of the individual DNA molecules to discover methylation patterns and nucleosome positioning.

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
devtools::install_github("rhondabacher/methylscaper", build_vignettes=TRUE, ref="R4.0")
```

# Dependencies

Note: on Ubuntu, users may need to install libgsl via:
`sudo apt-get install libgsl-dev`

The following packages are required for methylscaper. If installation fails, you may need to manually install the dependencies using the function 'install.packages' for CRAN packages or 'BiocManager::install' for Bioconductor packages.

* CRAN packages: shiny, shinyjs, shinyFiles, seriation, seqinr, Rfast, data.table
* Bioconductor packages: BiocParallel, Biostrings, SummarizedExperiment


For any other installation issues/questions please leave a message on our [Github Issues](https://github.com/rhondabacher/methylscaper/issues).

# Manuscript

A preprint of the methylscaper manuscript is now available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.11.13.382465v1).

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
