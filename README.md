# Romeb: Robust Median-Based Bayesian Approach

**Romeb** is an R package for robust median-based Bayesian modeling for longitudinal data with missingness (MAR, MNAR, MCAR). It uses JAGS and MCMC methods to estimate latent growth curves.

## 📦 Installation

You can install the development version of `Romeb` from GitHub:

```r
install.packages("devtools")
devtools::install_github("DandanTang0/Romeb")
```

## 🚀 Usage

```r
library(Romeb)

# Example data
data <- matrix(rnorm(300*5), 300, 5)

# Run model
result <- Romeb("MAR", data = data, seed = 123)

result
```

## 📖 Vignette

After installation, view the vignette:

```r
browseVignettes("Romeb")
```

## 🔧 Supported Missing Data Mechanisms

- MAR (Missing At Random)
- MNAR (Missing Not At Random)
- MCAR (Missing Completely At Random)
- No missing

## 📚 License

MIT License
