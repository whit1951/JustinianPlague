# Justinianic Plague

* `PlagueModelFunctions.R` -contains various ODE functions for hypothesized transmission routes of plague
* `PlottingTimeCourses.R`- plots time courses and comparison time courses for different ODE models

## Sensitivity Analysis

* `LHSNonuniform.R`- creates uniform and non-uniform LHS sampling space using `lhs` package
* `GlobalSensitivityAnalysis.R`- calculates PRCC values from LHS sampling; dependent on `LHSNonuniform.R`

## Rmd files
* `JustinianicPlagueModellingFigures.Rmd`- produces time courses of ODE models from `PlagueModelFunctions.R`; produces Fig. 1., Fig. S1 and Fig. S2
* `UniformLHSPRCC.Rmd`- Using uniform LHS distributions, creates boxplots, scatterplots, and PRCC graphs with confidence intervals using  `LHSNonuniform.R` and `GlobalSensitivityAnalysis.R`; produces Fig. S3
* `NonUniformLHSPRCC.Rmd`-Using nonuniform LHS distributions, creates boxplots, scatterplots, and PRCC graphs with confidence intervals using `LHSNonuniform.R` and `GlobalSensitivityAnalysis.R`; produces Fig. 2 and Figs S4-S7
