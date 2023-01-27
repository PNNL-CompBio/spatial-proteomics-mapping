# spatial-proteomics-mapping
Here we are evaluating the ability to measure spatially resolved proteomics in human spleen tissue. We are leveraging our previous work mapping [human pancreas]() to carry out similar evaluation and measurement.

## Data 

We will first collect the data and determine how many measurements we have.

## Analysis

### Spatial distribution of proteomics measurements

We will leverage two techniques to evaluate the distribution of proteomic and phosphoproteomic measurements across the tissue:

1. we will use PCA of both proteomics and phospho data
2. we will use k-means to identify putative clusters in the data
3. we will use the [BayesSpace]() package to determine how the samples cluster

### Tissue signatures
We will leverage the existing signatures to identify which proteins predict regions of the spleen so they can be leveraged on the spatial data. Here we will use elastic net to identify features that best predict the tissue types, then apply this model to the spatial data.

### Functional enrichment maps
There are multiple fronts to functional enrichment. We will likely leverage the leapR package for its increased flexibility.

1. functional enrichment stats (GSEA/Fishers/KSTAR) across predicted tissue types
2. functional enrichment stats across each region (this sounds unlikely to work)
3. networks combining significant proteins/kinases (PCSF)


## Figures
Ultimately we will use this repository to generate figures for the manuscript. These will be TBD later on.
