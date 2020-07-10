# methylscaper

`methylscaper` is a Shiny app for visualizing methylation data. Install the package with 

```{r}
devtools::install_github("pknight24/methylscaper", build_vignettes=TRUE)
```

Then run the Shiny app with

```{r}
library(methylscaper)
methylscaper::methylscaper()
```

To view the vignette run:
```{r}
browseVignettes('methylscaper')
```

The package also provides functions for ordering and plotting methylation data without needing to run the app.
