library(survminer)
library(survival)
library(ggplot2)
library(dplyr)
library(extrafont)
library(readxl)
library(stringr) 

df <- read_excel('')

df[, sapply(df, is.numeric) & names(df) != "AGEGRP"][is.na(df[, sapply(df, is.numeric) & names(df) != "AGEGRP"])] <- 0

get_mutation_group <- function(row) {
  epigenetic_mutations <- c("ASXL1", "BCOR", "BCORL1", "DNMT3A", "EZH2", "IDH1", "IDH2", "TET2")
  #cytogenetics_mutations <- c("inv16_t16.16", "t8.21", "minus.5", "del.5q.", "minus.7", "minus.17", "del.7q.")
  cohesion_mutations <- c("RAD21", "SMC1A", "SMC3", "STAG2")
  splicing_mutations <- c("SF3B1", "SRSF2", "U2AF1", "ZRSR2")
  signaling_mutations <- c("CBL", "CSF3R", "FLT3I", "FLT3T", "JAK2", "KIT", "KRAS", "NOTCH1", "NRAS", "PTPN11")
  transcription_mutations <- c("CEBPA", "CEBPADM", "CEBPA.bZIP", "CEBPA.bZIP.inframe", "CUX1", "GATA2", "IKZF1", "PHF6", "RUNX1", "WT1")
  TP53_mutation <- c("TP53")
  NPM1_mutation <- c("NPM1")
  
  mutation_groups <- c()
  
  if (any(row[epigenetic_mutations] == 1)) {
    mutation_groups <- c(mutation_groups, "Epigenetic")
  }
  #if (any(row[cytogenetics_mutations] == 1)) {
  #  mutation_groups <- c(mutation_groups, "Cytogenetics")
  #}
  if (any(row[cohesion_mutations] == 1)) {
    mutation_groups <- c(mutation_groups, "Cohesin")
  }
  if (any(row[splicing_mutations] == 1)) {
    mutation_groups <- c(mutation_groups, "Splicing")
  }
  if (any(row[signaling_mutations] == 1)) {
    mutation_groups <- c(mutation_groups, "Signaling")
  }
  if (any(row[transcription_mutations] == 1)) {
    mutation_groups <- c(mutation_groups, "Transcription")
  }
  if (any(row[TP53_mutation] == 1)) {
    mutation_groups <- c(mutation_groups, "TP53")
  }
  if (any(row[NPM1_mutation] == 1)) {
    mutation_groups <- c(mutation_groups, "NPM1")
  }
  return(paste(mutation_groups, collapse = "/"))
}
# Apply the function to create the "Mutation group" column
df$Mutation_group <- apply(df, 1, get_mutation_group)

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

mutation_list <-c("Epigenetic", "Cohesin", "Splicing", "Signaling", "Transcription", "TP53", "NPM1") #"Cytogenetics", 
age_group_list <- c('infants','children','AYA','adults','seniors','elderly')

base_directory <- ""

exploded_df = lapply(mutation_list, 
                     function(x) {
                       df_filtered <- df %>% select(Pat, Mutation_group, AGEGRP, OSSTAT, OSTM)
                       df_filtered <- df_filtered %>% filter(grepl(x,Mutation_group))
                       df_filtered$Mutation_group2 = x
                       df_filtered
                     }) %>% bind_rows()

# Define a custom color palette for each mutation group
mutation_colors <- c("Epigenetic" = "#6358A7",
                     #"Cytogenetics" = "#33FF57",
                     "Cohesin" = "#610178",
                     "Splicing" = "#35B779",
                     "Signaling" = "#21918C",
                     "Transcription" = "#FDE725",
                     "TP53" = "#90D743",
                     "NPM1" = "#31688E")

for (groups in age_group_list) {
  fit_list <- list()
  legend_labels <- c()
  
  for (genetic_name in mutation_list) {
    # Filter rows based on the AGEGRP column and then the Mutation_group column
    df_AGEGRP_filtered <- exploded_df %>% filter(str_detect(AGEGRP, groups), Mutation_group2 == genetic_name)
    
    # If there are not enough patients, skip the analysis for this mutation group
    if (nrow(df_AGEGRP_filtered) == 0) {
      cat("Skipping", genetic_name, "due to insufficient data\n")
      next
    }
    
    # Create a survival object for all mutation groups in the current age group
    Y <- Surv(df_AGEGRP_filtered$OSTM, df_AGEGRP_filtered$OSSTAT)
    fit <- survfit(Y ~ Mutation_group2, data = df_AGEGRP_filtered)
    fit_list[[genetic_name]] <- fit
    
    # Add legend label for the mutation group
    legend_labels <- c(legend_labels, genetic_name)
  }
  
  # Dynamically generate legend labels based on levels in AGEGRP
  df_AGEGRP_filtered <- exploded_df %>% filter(str_detect(AGEGRP, groups))
  p <- ggsurvplot_combine(fit_list, data = df_AGEGRP_filtered,
                          legend.title = "Mutation groups: ", palette = mutation_colors, pval = TRUE,
                          legend.labs = legend_labels,
                          xlab = "Time in months",
                          xlim = c(0, 204), break.time.by = 12, risk.table = TRUE, ggtheme = theme_bw(),
                          tables.height = 0.35, conf.int = FALSE, title = paste(groups))
  
  # Save the plot
  filename <- paste0(base_directory, groups, ".tiff")
  ggsave(filename, p, width = 9, height = 6, units = "in", dpi = 300)
}
