# apoptosisQuantification

The package `apoptosisQuantification` quantifies apoptosis and necroptosis
using predefined markers in proteomics datasets.

## Installation 

To install the package, enter

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install("tnaake/apoptosisQuantification")
```

## Methods to quantify apoptosis/necroptosis

The package contains three ways to quantify apoptosis/necroptosis:

1. by comparing the distribution of marker loadings and non-marker 
   loadings (from PCA),
2. by looking at the contribution of mitochrondrial proteins
   (this might be misleading when looking at a proteomics data set),
3. by running `decoupleR` and calculating a score for 
   apoptosis/necroptosis.

