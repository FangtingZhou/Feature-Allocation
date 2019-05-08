# Bayesian Double Feature Allocation

# Fangting Zhou

## Data and Method

### Abstract

The data used for the analysis consists of 12 infant samples with 6 breast feeding and 6 formula feeding. 
The data set contains the information of 162 operational taxonomic units aligned using rapid annotation using subsystem technology against the SEED subsystem database. 

### Availability

The data is not publicly available at present.

### Description

To explain the variation of microbial compositions across infants, Bayesian double feature allocation is proposed to deal with the count data matrix. 
The model will infer latent features that are associated with both operational taxonomic units and infants. 
At the same time, the result can be regarded as overlapping clustering for operational taxonomic units and infants simultaneously. 
To account for the compositional nature of data, the multinomial distribution will be incorporated into the model.

## Code

### Abstract

Analysis for this report is done with R. 
The corresponding code is provided to take exploratory data analysis, conduct various preprocessing steps, fit a Bayesian model via Markov chain Monte Carlo methods and generate descriptive plots used in the report.

### Description

All of the R scripts used in the report are available in a public repository on GitHub [https://github.com/FangtingZhou/Feature-Allocation].

### Optional Information

R version 3.5.3 is used for the analysis in this report. The necessary R libraries for the code used for data processing and analysis are

* ggplot2, version 3.1.0 (http://ggplot2.tidyverse.org)

* reshape2, version 1.4.3 (https://github.com/hadley/reshape)

* truncnorm, version 1.0-8 (https://github.com/olafmersmann/truncnorm)

## Instructions for Use

### Reproducibility

All data preparation and analyses are reproduced, as well as all Figures in the report.
