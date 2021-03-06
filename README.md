KRSA Manuscript Repo
================

# Introduction

This repo is for the KRSA application note.

# Analysis Script

The full analysis script is found in krsa\_manuscript\_analysis.Rmd

# Data

The datasets files are found in the data folder

## Parameters Used for this analysis:

-   chip type: “STK”

-   minimum signal: 5 (to filter out peptides with low signals)

-   *R*<sup>2</sup>: 0.9 (to filter out peptides with weak linear fit)

-   Log2 Fold change cutoffs: (0.2,0.3,0.4)

-   byChip: FALSE (across chip analysis)

-   KRSA function parameters (for upstream kinase analysis):

    -   mapping file: KRSA\_coverage\_STK\_PamChip\_87102\_v1 (built in
        the KRSA package)

    -   coverage file: KRSA\_coverage\_STK\_PamChip\_87102\_v1 (built in
        the KRSA package)

    -   iterations: 2000

    -   seed number: 123

-   Z score cutoff: 2
