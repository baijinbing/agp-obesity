1. Univariate analysis of BMI, diet, exercise, and weight change separately

1.1 Go to folder ‘univariate_BMI’
1.2 Run the shell script ‘qsub_step1_LDM.sh’ in a linux cluster, which will submit 42 jobs with different seed values. Make sure you create a subfolder “tmp_files” for storing intermediate results. Note that each job will take about 2-3 hours. 
1.3 Run the shell script ‘qsub_step2_summary.sh’ in the linux cluster, which reads in the intermediate results.
1.4 repeat 1.1-1.3 for folders ‘univariate_diet’, ‘univariate_exercise’, and ‘univariate_wt’


2. Multivariate analysis of all predictors

2.1 Go to folder ‘multivariate_all_predictors’
2.2 Same as 1.2
2.3 Same as 1.3
2.4 Crop results for detected OTUs associated with BMI from ‘LDM_summary.Rout’ and save them in ‘LDM_multi_BMI.txt’ (which exists in the folder). Similar procedure for diet and exercise. 
2.5 Run ‘Enrichment_detectedOTUs_BMI.r’ for enrichment analysis of OTUs found to be associated with BMI based on results saved in ’LDM_multi_BMI.txt’. Similar procedure for diet and exercise.  
 