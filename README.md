# *Hapi* - Chromosome-length haplotype inference using genotypic data of single gamete cells

#### Please cite the following publication if the data or code in the repostiory is used in your study:
Li, R., Qu, H., Chen, J. et al. (2020). Inference of Chromosome-Length Haplotypes Using Genomic Data of Three or a Few More Single Gametes. *[Molecular Biology and Evolution](https://doi.org/10.1093/molbev/msaa176)*


## Introduction

`Hapi` is a novel easy-to-use and high-efficient algorithm that only requires **3 to 5 gametes** to reconstruct **accurate and high-resolution haplotypes** of an individual. The gamete genotype data may be generated from **various platforms including genotyping arrays and sequencing even with low-coverage**. `Hapi` simply takes genotype data of known hetSNPs in single gamete cells as input and report the high-resolution haplotypes as well as the confidence level of each phased hetSNPs. The package also includes a module allowing **downstream analyses and visualization of identified crossovers** in the gametes. 


## Manual and R script
The comprehensive manual of `Hapi` is available here: [Hapi Manual](http://htmlpreview.github.io/?https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/Hapi_manual.html).  
R code of the workflow is available here: [Hapi Workflow](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/Hapi_workflow.R)


## Installation

```R
install.packages('Hapi')
```
