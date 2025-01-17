
## Spatial-statistical downscaling

This is the R package **spds** (currently developer's version) for the paper:

Zheng, X., Cressie, N., Clarke, D. A., McGeoch, M. A., and Zammit-Mangion, A. (2025). 
"Spatial-statistical downscaling with uncertainty quantification in biodiversity modelling". To appear in *Methods in Ecology and Evolution*.

You can install the package with **devtools**
```
devtools::install_github("xzheng42/corgi-examples/", subdir = "spds")
library(spds)
```

The main functions of the package are `gpbau` and `predict.gpbau`:

- `gpbau` fits a GPB (Gaussian process over Basic Areal Units) via Markov chain Monte Carlo (MCMC).
- `predict.gpbau` (or simply `predict`) generates posterior predictive samples of the downscaled covariate given a gpbau object.

Detailed guidelines for using the functions are referred to their help pages in R. R scripts to reproduce results in the paper are available in 
[*data-examples/*](https://github.com/xzheng42/corgi-examples/tree/main/data-examples) 
with instructions available in [*corgi-examples/*](https://github.com/xzheng42/corgi-examples).

Notes:

- The current version was tested on macOS 10.15.7 under R version 4.2.2 and on Fedora Linux 38 under R versions 4.3.3.

- Instead of installing the package, you can download the package and load it in R via devtools::load_all("~/path/to/spds/").

### Reporting bugs and issues

If you find a bug or issue, please use the Github Issues (i.e., click the `Issues` tab, click the green `New issue` button and create a new issue). It is much appreciated if you could include the following in the issue: 

- A description of the bug or issue;
- A minimal reproducible example (e.g., source code and data files) with your comments.
