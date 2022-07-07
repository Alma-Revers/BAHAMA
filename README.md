# BAHAMA
This repository contains a R-package, tutorial and the Bayesian model as presented in the paper: BAHAMA: A Bayesian hierarchical model for the detection of MedDRA coded adverse events in RCTs

## Structure of the repository:
- The BAHAMA R package
- Simulation study as presented in the paper
- Stan code data file of the BAHAMA model

## How to install the BAHAMA package
Install devtools (if needed)
```
install.packages("devtools")
```

Load the devtool package
```
library(devtools)
```

Now you can install the BAHAMA package by using:
```
install_github("Alma-Revers/BAHAMA", subdir="BAHAMA")
```

## How to use the BAHAMA package
Load the BAHAMA package after installation
```
library(BAHAMA)
```
The BAHAMA method needs 3 dataset, (1)the adverse event counts, (2)the medDRA tree, and (3)data from the patient/cases.
Example data: 
1) AE_data, a vector with AEs
2) tree, a tree with the medDRA structure of the example AE data
3) sampleData, a vector with a treatment group of the cases

1) Adverse event data needs to be a cases by AE matrix. If counts are ordered in a long format, use the function:
```
counts <- aeCountsFromVector(AE_data)
```
2) A matrix of MedDRA IDs. The rows correspond to the MedDRA structure (SOC, HLGT, HLT, PT), each row correspond to a single AE. If multiaxiallity is present, a single AE is represented in multiple rows.
```
meddra_dataset <- medDRADataFromMatrix(tree[,tree[4,] %in% colnames(counts)])
```

3) A data.frame with information about the cases/patients, which treatment group (x)

Create a BAHAMA dataset which is the input for the Bayesian model in stan. The function BAHAMADataSet can be used. The threshold for the incidence of the AE and the threshold on the strucutre can be set here as parameters.
```
BAHAMA_dataSet <- BAHAMADataSet(aeCounts = counts, medDRAData = meddra_dataset, sampleData = as.data.frame(sampleData),
                                thresholdIncidence = 5, thresholdStructure = 5)
```

Next is to drawn samples of the posterior. 
```
BAHAMA_results <- BAHAMA(BAHAMA_dataSet)
```
In case of a warning about convergency issues, we recommend the following;
- Change the number of iteration and other settings of STAN, see stan documentation
- Change the thresholds

Lastly, results can be shown within a shiny app. 
```
plotBAHAMA(BAHAMA_results$tidy_rr_pt, BAHAMA_results$tidy_rr_hlt, BAHAMA_results$tidy_rr_hlgt, BAHAMA_results$tidy_rr_soc)
```


