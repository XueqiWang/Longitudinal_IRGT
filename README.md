The current folder includes R code for reproducing all of the tables and figures in the article “Designing individually randomized group treatment trials with repeated outcome measurements” by Wang et al.

For questions or comments about the code please contact Xueqi Wang at xueqi.wang@duke.edu.

I. List of Supporting Files: These supporting files are sourced in the main files that reproduce the numbers in the submitted manuscript.

1. LIRGT_SampleSize.R = function to calculate sample size for longitudinal IRGT trials with a continuous or binary outcome;
2. LIRGT_Power.R = function to calculate power for longitudinal IRGT trials with a continuous or binary outcome;
3. conGEN.R = function to simulate continuous outcome data in longitudinal IRGT trials;
4. conMAEE.R = function to calculate the bias-corrected sandwich variances for GEE/MAEE analyses with continuous outcomes under Model 1-3;
5. conMAEE_F.R = function to calculate the bias-corrected sandwich variances for GEE/MAEE analyses with continuous outcomes under Model 4-5;
6. binGEN.R = function to simulate binary outcome data in longitudinal IRGT trials;
7. binMAEE.R = function to calculate the bias-corrected sandwich variances for logistic-binomial GEE/MAEE analyses under Model 1-3;
8. binMAEE_F.R = function to calculate the bias-corrected sandwich variances for logistic-binomial GEE/MAEE analyses under Model 4-5;
9. geesmv_new.R = function to calculate the bias-corrected sandwich variances for GEE analyses with the independence working correlation model with continuous or binary outcomes.

II. List of Main Files: These main files are used to reproduce the results in the submitted manuscript.

10. conScenarios_Analytical.R = reproduce predicted power results for continuous outcomes;
11. conScenarios_1.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 1 with continuous outcomes;
12. conScenarios_2.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 2 with continuous outcomes;
13. conScenarios_3.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 3 with continuous outcomes;
14. conScenarios_4.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 4 with continuous outcomes;
15. conScenarios_5.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 5 with continuous outcomes;
16. conDATA.R = collect power/size results from raw data after running the above simulation programs for continuous outcomes;
17. conRES.R = reproduce numbers and the format of the table and figures for continuous outcomes;
18. binScenarios_Analytical.R = reproduce predicted power results for binary outcomes;
19. binScenarios_1.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 1 with binary outcomes;
20. binScenarios_2.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 2 with binary outcomes;
21. binScenarios_3.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 3 with binary outcomes;
22. binScenarios_4.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 4 with binary outcomes;
23. binScenarios_5.R = reproduce simulation results of GEE analyses using an independence working correlation matrix and GEE/MAEE analyses using the true working correlation structure, for Model 5 with binary outcomes;
24. binDATA.R = collect power/size results from raw data after running the above simulation programs for binary outcomes;
25. binRES.R = reproduce numbers and the format of the table and figures for binary outcomes;
26. SYV.R = reproduce results of the application;

III. Folder

27. conResults = folder to save power/size data and results for continuous outcomes;
28. binResults = folder to save power/size data and results for binary outcomes.

IV. Software

Analyses were conducted with R, version 4.2.2 (https://www.r-project.org/). The calculations used R packages mvtnorm (version 1.1-3), MASS (version 7.3-58.1), gee (version 4.13-25), geesmv (version 1.3), openxlsx (version 4.2.5.1), reshape2 (version 1.4.4), ggplot2 (version 3.4.1), ggpubr (version 0.5.0), scales (version 1.2.1), pracma (version 2.4.2), directlabels (version 2021.1.13), and cowplot (version 1.1.1).

V. R commands for the installation of R packages

install.packages(c("mvtnorm", "MASS", "gee", "geesmv", "openxlsx", "reshape2", "ggplot2", "ggpubr", "scales", "pracma", "directlabels", "cowplot"))

NOTE: Make sure the current working directory is the current folder before running each program.
