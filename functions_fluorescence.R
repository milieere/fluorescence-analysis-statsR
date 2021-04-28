#'
#'This is  set of functions used by analyze.R to summarize mean intensity 



###############################################################################################################


#' load_data takes the directory with csv files from replicates outputted by Python pipeline
#' with YeastSpotter.
#' @param results_dir takes the directory with csv files 
#' @param normalize if set to true, max value per replicate will be choosen and the rest
#' of values will be divided with the max value (to normalize between replicates with different exposure 
#' and compare relatively).
#' @return dataset that has concatenated data from all replicates 
load_data <- function(results_dir, normalize = FALSE){
  library(dplyr)
  
  list <- list.files(results_dir)
  dataset <- data.frame()
  
  for(i in 1:length(list)){
    if (length(list)==1){
      temp_file <- read.csv(file=list[i], header=TRUE, sep=',', na.strings='NA', stringsAsFactors = FALSE)
      dataset <- temp_file
    } else {
      temp_file <- read.csv(file=list[i], header=TRUE, sep=',', na.strings='NA', stringsAsFactors = FALSE)
      temp_file$rep_no <- c(rep(list[i], nrow(temp_file)))
      dataset <- bind_rows(dataset, temp_file)
      dataset <- dataset %>% group_by(rep_no, filename) %>% mutate(no_cells = n())
    }
  }
  
  #Normalize all the mean intensity values (0:1) if normalize is set to TRUE 
  if(normalize == TRUE){
    dataset <- dataset %>% group_by(rep_no) %>% dplyr::mutate(maximum=max(mean_intensity))
    dataset$normalized_mean_int <- dataset$mean_intensity / dataset$maximum
  }
    
  return(dataset)
}

###############################################################################################################

#' Look on the distribution of data in each sample, if they are homogeneous so we can compare them statistically
#' @param results takes the results dataframe outputted by load_data
#' @return distributions - plots of individual distributions of each sample
checkdistr <- function(results, normalized=FALSE){
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  
  if(normalized == TRUE){
    meanint <- 'normalized_mean_int'
  } else {
    meanint <- mean_intensity
  }

  makedistr <- function(results, meanint){
    distributions <- ggplot(results, aes_string(x=meanint))+
      geom_density(aes(group=interaction(filename, rep_no), color=rep_no))+
      theme_bw()+
      facet_grid(filename~.)
    
    return(distributions)
  }
  
  distributions <- makedistr(results, meanint)
  return(distributions)
  
}

###############################################################################################################

#' This function computes means of each protein per replicate and checks data (ab)normality hehe
#' @param results takes the dataframe outputted by load_data
#' @return summary with means per protein (condition) and replicate.
summarize_data <- function(results, normalized = FALSE){
  library(dplyr)
  library(ggpubr)
  
  if(normalized == TRUE){
    summary_replicates <- results %>% group_by(rep_no, filename) %>% 
      summarise(mean=mean(normalized_mean_int), 
                sd = sd(normalized_mean_int), 
                stat = shapiro.test(normalized_mean_int)$statistic, 
                p = shapiro.test(normalized_mean_int)$p.value, 
                n=mean(no_cells))
  } else {
    summary_replicates <- results %>% group_by(rep_no, filename) %>% 
      summarise(mean=mean(mean_intensity), 
                sd = sd(mean_intensity), 
                stat = shapiro.test(mean_intensity)$statistic, 
                p = shapiro.test(mean_intensity)$p.value, 
                n=mean(no_cells))
  }
  
  return(summary_replicates)
}

###############################################################################################################

#' This function computes means of means from replicates and their standard deviation.
#' @param summary takes summary file outputted by summarize_data
#' @return summary with means per protein (condition) and replicate.
summarize_conditions <- function(summary){
  library(dplyr)
  #Makes summary from means of replicates per condition
  summary_conditions <- summary %>% group_by(filename) %>% summarise(mean_replicates = mean(mean), sd=sd(mean))
  
  return(summary_conditions)
}

###############################################################################################################

#' Performs Tukey HSD test for pairwise comparisons between proteins.
#' @param summary takes summary file outputted by summarize_data
#' @return stat.test with pairwise comparisons
make_stats <- function(summary){
  #Perform statistical test (TukeyHSD)
  library(rstatix)  # https://github.com/kassambara/rstatix
  
  stat.test <- aov(mean ~ filename, data = summary) %>%
    tukey_hsd()
  
  return(stat.test)
}

###############################################################################################################

#' Plots the data.
#' @param results takes the results dataframe outputted by load_data
#' @return stat.test with pairwise comparisons
plot_data_meanint <- function(results, normalized = FALSE){
  library(ggplot2)
  library(dplyr)
  
  if(normalized==FALSE){
    plot <- ggplot(results, aes(x=filename, y=mean_intensity, color = rep_no)) + 
      geom_boxplot(colour = "grey50") +
      geom_point(position='jitter', alpha=0.3, size = 0.7) + 
      stat_summary(fun=mean, geom="point", shape=18,
                   size=6, aes(color=rep_no))
  } else {
    plot <- ggplot(results, aes(x=filename, y=normalized_mean_int, color = rep_no)) + 
      geom_boxplot(colour = "grey50") +
      geom_point(position='jitter', alpha=0.3, size = 0.7) + 
      stat_summary(fun=mean, geom="point", shape=18,
                   size=6, aes(color=rep_no))
  }
  return(plot)
}

###############################################################################################################


