# Romeb: Robust Median-Based Bayesian Approach

**Romeb** is an R package for robust median-based Bayesian modeling for longitudinal data with missingness (MAR, MNAR, MCAR). It uses JAGS and MCMC methods to estimate latent growth curves.

## 📦 Installation

You can install the development version of `Romeb` from GitHub:

```r
install.packages("devtools")
devtools::install_github("DandanTang0/Romeb", dependencies = TRUE)
```

## 🚀 Usage

```r
library(Romeb)

# Example data
data <- matrix(rnorm(300*5), 300, 5)

# Run model
result <- Romeb("MAR", data = data, Time = 5, seed = 123)

result
```

## 📖 Vignette

Installation with vignette

```r
devtools::install_github(
  "DandanTang0/Romeb",
  build_vignettes = TRUE,   
  dependencies    = TRUE,
  force           = TRUE  # replace previous installation
)
```

After installation, view the vignette:

```r
browseVignettes("Romeb")
```

Users can read or rerun the full, documented example.

## 🔧 Supported Missing Data Mechanisms

- MAR (Missing At Random)
- MNAR (Missing Not At Random)
- MCAR (Missing Completely At Random)
- No missing

## 📚 License

MIT License
