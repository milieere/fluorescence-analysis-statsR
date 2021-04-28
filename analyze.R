#'
#' 28/04/2021
#' SIMONA BARANKOVA
#'
#' This script serves to analyze relative protein expression from fluorescent images, 
#' and performing basic statistics. 
#' 
#' Naming of images convention → Protein_otherinfo, i.e. Sec6_GFP_e1000.tif
#' 
#' Takes as input the folder with csv files - results from Python analysis. One csv = 1 sample
#' (one day of imaging). For information about data extraction pipeline see github:
#' https://github.com/simonabarankova/fluorescence-analysis.
#' 
#' This script normalizes the data, calculates means of the samples and then mean of
#' the population with standard error. 
#' 
#' Outputs plots of the mean intensity (boxplot) with error bars - standard error.

####################### SET HERE ########################################################

## Set the result_dir to the path to the directory with your results csv's
## Set the script_dir to the path of the functions_fluorescence.R script

results_dir <- 'C:/Users/Simona/Documents/discoscience/doctorate_upf/bioinfo/protein_expression/results'
script_file <- 'C:/Users/Simona/Documents/discoscience/doctorate_upf/bioinfo/protein_expression/analysis/functions_fluorescence.R'

####################### SET HERE ########################################################

#Load the functions
source(script_file, local=T)

#Load the data
results <- load_data(results_dir, normalize = TRUE)

#Check distribution of samples
distr <- checkdistr(results, normalized=TRUE)

#Compute the statistics 
summary <- summarize_data(results, normalized = TRUE)

#Compute final statistics (mean of the replicate means)
summary_fin <- summarize_conditions(summary)

#Plot the data, mean intensities and distributions of cell areas (indicative for growth rate)
plots_meanint <- plot_data_meanint(results, normalized = TRUE)

#Make pairwise comparisons, Tukey HSD to find out what samples are significantly different
stat <- make_stats(summary)

plots_areadistr <- plot_data_area(results)




