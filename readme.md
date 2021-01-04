# Larger fish ontogeny shapes ecological network structure

### Authors: Korryn Bodner, Chris Brimacombe, Marie-Josée Fortin and Péter K. Molnár
\
## About the Project
This repository contains code to recreate the analyses for "Why body size matters: how larger fish ontogeny shapes ecological network structure". The purpose of these analyses were to construct and compare ecological interaction networks, specifically, non-stage-structured (adults only) and stage-structured (adults and juveniles) freshwater stream fish interaction networks. Data for the analyses can be found in the repository's data folder and is sourced from NEON's [Fish electrofishing, gill netting, and fyke netting counts](https://data.neonscience.org/data-products/DP1.20107.001). There are three main sections to the analysis: Network Construction, Network Comparison (3D Plot and Matrices & GCD-11 and GCM-11), and Dissimilarity Analysis.  

## Getting Started

Code for Network Construction, 3D Plot and Matrices, and Dissimilarity Analysis were run in R 4.0.3 while GCD-11 and GCM-11 Calculations were performed in Python 2.7.13. Please make sure you have access to R 4.0 or above and Python 2.7 to run the analyses below.  

## Network Construction

### *Overview:* 
To construct inferred adult and stage-structured species interaction networks, we use the packages [EMtree](https://rmomal.github.io/EMtree/) (Momal et al. 2020) which draws on [PLNmodels](http://julien.cremeriefamily.info/PLNmodels/) (Chiquet et al. 2018, 2019) to create networks. We also construct random networks using a rewiring procedure found in [igraph](https://igraph.org/) (Csardi and Nepusz. 2006). Code to construct networks can be found in "Network_Construction_PLN_EMtree_Rewiring.R" and is organized into the following four steps:

1. organize data
2. build PLNmodels
3. construct resampled networks with EMtree
4. create randomly rewired networks

#### *Code:*<span style="font-weight:normal"> "Network_Construction_PLN_EMtree_Rewiring.R" </span> 

#### *Data:*<span style="font-weight:normal"> "NEON_aggregate_adultOnly.csv" &   "NEON_aggregate_stageStructured.csv" </span>

## Network Comparison

### GCD-11 and GCM-11 Calculations

### *Overview:*

Graphlet correlation distance-11 (GCD-11) and graphlet correlation matrix-11 (GCM-11) are methods adopted from "Revealing the Hidden Language of Complex Networks" by Yaveroglu et al. (2014). Python code to calculate GCD-11 and the components of GCM-11 can be found [here](http://www0.cs.ucl.ac.uk/staff/natasa/GCD/).

### 3D Plot and Matrices (for GCD-11 and GCM-11)

### *Overview:*
Here  GCD-11 and GCM-11 outputs are used as inputs to create summary measures and figures. The code creates a 3D Plot using GCD-11 scores reduced to three dimensions via Metric-Multidimensional Scaling with the [car](https://cran.r-project.org/web/packages/car/index.html) package (Fox and Weisberg 2019). To create GCM-11 figures, .ndump2 files are read in, spearman correlations are calculated, and correlation coefficients and significance levels are plotted.


#### *Code:*<span style="font-weight:normal"> "Network_Comparison_MDS_GCM-11.R" </span> 

#### *Data:*<span style="font-weight:normal"> "rewire_adult" & "rewire_fifty" with ".ndump2" extension </span> 

## Dissimilarity Analysis

### *Overview:*
The following code creates a Jaccard dissimilarity index of connections with other species between "larger" adults and juveniles, and creates a linear regression and corresponding plot of adult and juvenile size differences with their jaccard dissimilarity scores.

#### *Code:*<span style="font-weight:normal"> "Dissimilarity_analysis.R"</span> 
#### *Data:*<span style="font-weight:normal"> "NEON_fish_inds.csv"</span>  
<div style="line-height:10%;">
    <br>
</div>

## References

Chiquet, J., M. Mariadassou, and S. Robin. 2018. Variational inference for probabilistic poisson PCA. Annals of Applied Statistics 12:2674–2698.

Chiquet, J., M. Mariadassou, and S. Robin. 2019. Variational inference of sparse network from count data. Pages 1988–1997 36th International Conference on Machine Learning, ICML 2019.

Csardi, G., and T. Nepusz. 2006. The igraph software package for complex network research. InterJournal Complex Systems 1695.

Fox, J., and S. Weisberg. 2019. An R Companion to Applied Regression. SAGE Publications Inc., Thousand Oaks, CA, USA.

Momal, R., S. Robin, and C. Ambroise. 2020. Tree-based inference of species interaction networks from abundance data. Methods in Ecology and Evolution 11:621–632.

National Ecological Observatory Network 2020. Data product DP1.20107.001, fish electrofishing, gill netting and fyke netting counts. – Provisional data accessed 1 Nov. 2020.

Yaveroğlu, Ö. N., et al. 2014. Revealing the hidden Language of complex networks. Scientific Reports 4:1–9.



