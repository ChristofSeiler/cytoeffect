# cytoeffect

<!-- badges: start -->

[![Build Status](https://app.travis-ci.com/ChristofSeiler/cytoeffect.svg?branch=master)](https://app.travis-ci.com/ChristofSeiler/cytoeffect)
<!-- 
[![R build status](https://github.com/ChristofSeiler/cytoeffect/workflows/R-CMD-check/badge.svg)](https://github.com/ChristofSeiler/cytoeffect/actions)
-->
[![Codecov test coverage](https://codecov.io/gh/ChristofSeiler/cytoeffect/branch/master/graph/badge.svg)](https://codecov.io/gh/ChristofSeiler/cytoeffect?branch=master)

<!-- badges: end -->

## Goal

Regression analysis with multivariate outcomes for mass cytometry experiments.

## Installation

Before you can install `cytoeffect`, you will need to install the package `rstan` including the C++ Toolchain. Here are instructions for all main platforms:

https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

Then, you can install `cytoeffect`:

```r
install.packages("devtools")
devtools::install_github("ChristofSeiler/cytoeffect")
```
