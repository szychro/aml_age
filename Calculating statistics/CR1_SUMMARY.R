library(openxlsx)
library(survival)
library(readxl)

df <- read_excel('/Users/szymon/Desktop/Praca/Seattle project/Final_data/141223_target_sal_ohsu_modJE_20231227.xlsx')
wb <- createWorkbook()
endpoint <- 'CR1'

df$CEBPA.bZIP.inframe2 <- NA
for (i in 1:length(df$CEBPA)) {
  if (is.na(df$CEBPA[i])) {
    df$CEBPA.bZIP.inframe2[i] <- NA
  } else if (df$CEBPA[i] == 0) {
    df$CEBPA.bZIP.inframe2[i] <- 0
  } else if (df$CEBPA[i] == 1) {
    df$CEBPA.bZIP.inframe2[i] <- ifelse(df$CEBPA.bZIP.inframe[i] == 0, 0, 1)
  }
}

df <- subset(df, select = -c(CEBPA.bZIP.inframe))

names(df)[names(df) == "CEBPA.bZIP.inframe2"] <- "CEBPA.bZIP.inframe"

mutations <- c("ASXL1", "BCOR", "BCORL1", "CBL", "CEBPA", "CSF3R", "CUX1", "DNMT3A", "ETV6", "EZH2", "FLT3I",
               "GATA2", "IDH1", "IDH2", "IKZF1", "JAK2", "KIT", "KRAS", "NPM1", "NRAS", "PHF6", "PTPN11", "RAD21", "RUNX1",
               "SF3B1", "SMC3", "SRSF2", "STAG2", "TET2", "TP53", "U2AF1", "WT1","ZRSR2", "inv16_t16.16", "t8.21", "minus.5", "del.5q.",
               "minus.7", "minus.17", "t.v.11..v.q23.", "t.9.11..p21.23.q23.", "t.10.11.", "t.11.19..q23.p13.",
               "del.7q.", "del.9q.", "trisomy 8", "trisomy 21", "minus.Y", "minus.X", 'CEBPA.bZIP.inframe')

pvalue_mutations <- c()
unique_age_groups <- c("infants", "children", "AYA", "adults", "seniors", "elderly")

for (variable in mutations) {
  for (age_group in unique_age_groups) {
    
    age_filter <- df$AGEGRP == age_group
    
    filtered_df <- df[age_filter & !is.na(df[[endpoint]]), ]
    
    if (sum(!is.na(filtered_df[[variable]])) > 0) {
      
      number_patients <- sum(filtered_df[[variable]] == 1, na.rm = TRUE)
      
      y_age <- filtered_df[[endpoint]]
      x_age <- filtered_df[[variable]]

      if (sum(!is.na(x_age)) > 1) {
        model_age <- glm(y_age ~ x_age, family = binomial)
        model_summary_age <- summary(model_age)
        
        p_value_age <- coef(summary(model_age))[, "Pr(>|z|)"]
        pvalue_mutations <- c(pvalue_mutations, p_value_age[2])
        
      }
    }
  }
}

fdrs <- p.adjust(pvalue_mutations, method = "BH")
print(fdrs)

x <- 1  

for (variable in mutations) {
  # Create a worksheet for each mutation
  ws <- addWorksheet(wb, variable)
  headers <- c("Age Group", "OR", "Lower CI", "Upper CI", "p-value", "p_adj")
  writeData(wb, variable, headers, startRow = 1, startCol = 1)
  row_counter <- 1
  
  for (age_group in unique_age_groups) {
    
    age_filter <- df$AGEGRP == age_group
    filtered_df <- df[age_filter & !is.na(df[[endpoint]]), ]
    
    # Check if there are any non-NA values for the current variable
    if (sum(!is.na(filtered_df[[variable]])) > 0) {
      
      number_patients <- sum(filtered_df[[variable]] == 1, na.rm = TRUE)
      number_patients_0 <- sum(filtered_df[[variable]] == 0, na.rm = TRUE)
      sum_patients = number_patients + number_patients_0
      
      y_age <- filtered_df[[endpoint]]
      x_age <- filtered_df[[variable]]

      # Check if there are enough non-NA values to fit the Cox model
      if (sum(!is.na(x_age)) > 1) {
        model_age <- glm(y_age ~ x_age, family = binomial)
        model_summary_age <- summary(model_age)
        
        HR_age <- exp(coef(model_age))
        CI_age <- exp(confint(model_age, level = 0.95))
        p_value_age <- coef(summary(model_age))[, "Pr(>|z|)"]
        
        p_adj_value <- if (number_patients > 0) fdrs[x] else ""
        x <- x + 1  
        
        coef_table_age <- data.frame(
          Age_Group = age_group,
          Num_Pat = number_patients,
          Total_pat = sum_patients,
          OR = HR_age[2],
          Lower_CI = CI_age[2, 1],
          Upper_CI = CI_age[2, 2],
          p_value = p_value_age[2],
          p_adj = p_adj_value
        )
        
        # Write the coefficient table to the worksheet
        writeData(wb, variable, coef_table_age, startRow = row_counter, startCol = 1)
        style <- createStyle(fgFill = "#C7E0FC")
        addStyle(wb, sheet = variable, rows = row_counter, cols = c(1, 2, 3, 4, 5, 6, 7,8), style = style)
        
        row_counter <- row_counter + 2
      }
    }
  }
}

# Save the workbook
output_filename <- ''
saveWorkbook(wb, output_filename, overwrite = TRUE)
