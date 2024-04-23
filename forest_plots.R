#in filename should be made a file created via SUMMARY code
#change the color palette depending on the endpoint is loaded
#in OR_age the value should be change to either OR/HR depending on the endpoint used

library(ggplot2)
library(dplyr)
library(forestplot)
library(grid)
library(survival)
library(survminer)
library(forestploter)
library(viridisLite)
library(viridis) 
library(readODS)
library(formattable)
library(htmltools)

# Function to read all sheets from an Excel file
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if (!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#palette for cr1
color_palette <- c("#FF0000","#E2001C", "#DF001F", "#C60038","#AA0055" ,"#8D0071", "#71008D","#5500AA","#3800C6","#1C00E2", "#0000FF" )
#palette for ostm/efstm
#color_palette <- c("#0000FF", "#1C00E2", "#3800C6", "#5500AA",  "#71008D", "#8D0071", "#AA0055", "#C60038",  "#DF001F","#E2001C", "#FF0000") 

filename <- ""
mysheets <- read_excel_allsheets(filename)
forest_data <- data.frame()

cytogenetics <-  c("inv16_t16.16", "t8.21", "minus.5", "del.5q.",
                   "minus.7", "minus.17", "t.v.11..v.q23.", "t.9.11..p21.23.q23.", "t.10.11.", "t.11.19..q23.p13.",
                   "del.7q.", "del.9q.", "trisomy 8", "trisomy 21", "minus.Y", "minus.X")

molecular_genetics <- c("ASXL1", "BCOR", "BCORL1", "CBL", "CEBPA","CSF3R", "CUX1", "DNMT3A", "ETV6", "EZH2",
                        "FLT3I", "GATA2", "IDH1", "IDH2", "IKZF1","JAK2", "KIT", "KRAS", "NPM1", "NRAS",
                        "PHF6", "PTPN11", "RAD21", "RUNX1", "SF3B1","SMC3", "SRSF2", "STAG2", "TET2", "TP53",
                        "U2AF1", "WT1", "ZRSR2", "CEBPA.bZIP.inframe")
endpoint <- 'HR'
for (sheet_name in names(mysheets)) {
  if (sheet_name %in% molecular_genetics) {
    next  # Skip the iteration if sheet_name is in molecular genetics
  } else {
    sheet <- mysheets[[sheet_name]]
    number_patients <- sheet$Num_Pat[c(1, 3, 5, 7, 9, 11)]
    total_patients <- sheet$Total_pat[c(1, 3, 5, 7, 9, 11)]
    OR_age <- sheet$OR[c(1, 3, 5, 7, 9, 11)] # OR for CR1, HR for OSTM/EFSTM
    CI_age <- matrix(c(sheet$Lower_CI[c(1, 3, 5, 7, 9, 11)], sheet$Upper_CI[c(1, 3, 5, 7, 9, 11)]), ncol = 2)
    p_value_age <- sheet$p_value[c(1, 3, 5, 7, 9, 11)]
    p_adj <- sheet$p_adj[c(1, 3, 5, 7, 9, 11)]
    age_group_save <- sheet$Age_Group[c(1, 3, 5, 7, 9, 11)]
    
    row_data <- data.frame(
      feature = rep(sheet_name, 6),
      AGEGRP = age_group_save,
      Num_pat = as.double(number_patients),
      tot_pat = as.double(total_patients),
      HR = as.double(OR_age),
      CI_low = as.double(CI_age[, 1]),
      CI_high = as.double(CI_age[, 2]),
      pvalue = as.double(p_value_age),
      p_adjusted = as.double(p_adj)
    )
    
    forest_data <- rbind(forest_data, row_data)
  }
}
forest_data$feature <- gsub("t8.21", "t(8;21)", forest_data$feature)  # Change "t8.21" to "t(8;21)"
forest_data$feature <- gsub("minus.5", "-5", forest_data$feature)  
forest_data$feature <- gsub("minus.7", "-7", forest_data$feature)  
forest_data$feature <- gsub("minus.17", "-17", forest_data$feature)  
forest_data$feature <- gsub("del.5q.", "del(5q)", forest_data$feature) 
forest_data$feature <- gsub("t.v.11..v.q23.", "t(v;11q23.3)", forest_data$feature)
forest_data$feature <- gsub("t.9.11..p21.23.q23.", "t(9;11)(p21.3;q23.3)", forest_data$feature)  
forest_data$feature <- gsub("t.10.11.", "t(10;11)", forest_data$feature)  
forest_data$feature <- gsub("t.11.19..q23.p13.", "t(11;19)(q23.3;p13.3)", forest_data$feature) 
forest_data$feature <- gsub("del.7q.", "del(7q)", forest_data$feature)
forest_data$feature <- gsub("del.9q.", "del(9q)", forest_data$feature) 
forest_data$feature <- gsub("minus.X", "-X", forest_data$feature)  
forest_data$feature <- gsub("minus.Y", "-Y", forest_data$feature)  
forest_data$feature <- gsub("inv16_t16.16", "inv(16)(p13.1q22)", forest_data$feature)
forest_data$feature <- gsub("CEBPA.bZIP.inframe", "CEBPA-bZip-inf", forest_data$feature) 
forest_data$feature <- gsub("FLT3I", "FLT3-ITD", forest_data$feature) 

forest_data <- na.omit(forest_data)
forest_data$CI_all <- paste(sprintf('%.2f', forest_data$CI_low), " - ", sprintf('%.2f', forest_data$CI_high))
forest_data$pvalue <- ifelse(forest_data$pvalue <= 0.001, "< 0.001", sprintf("%.4f", round(forest_data$pvalue, 4)))
forest_data <- forest_data[order(forest_data$p_adjusted), ]

plot_data <- data.frame()
output_directory <- ""
age_groups = c('infants', 'children', 'AYA', 'adults', 'seniors', 'elderly')

for (variable in unique(forest_data$feature)) {
  filtered_data <- subset(forest_data, feature == variable)
  filtered_data <- filtered_data[complete.cases(filtered_data), ]
  
  for (i in 1:nrow(filtered_data)) {
    current_row <- filtered_data[i, ]
    row_data <- data.frame(
      age_group = current_row$feature,  
      AGEGRP = as.character(current_row$AGEGRP),
      Num_pat = current_row$Num_pat,
      Tot_pat = current_row$tot_pat,
      HR = current_row$HR,
      CI_low = current_row$CI_low,
      CI_high = current_row$CI_high,
      CI = current_row$CI_all,
      pvalue = current_row$pvalue,
      p_adj = current_row$p_adjusted
    )
    plot_data <- rbind(plot_data, row_data)
  } 
  
  if (variable == tail(unique(forest_data$feature), n = 1)) {
    plot_data$HR <- format(round(plot_data$HR, 2), nsmall = 2)
    plot_data$AGEGRP <- factor(plot_data$AGEGRP, levels = age_groups)
    plot_data <- plot_data[order(plot_data$AGEGRP, plot_data$HR), ]
    
    filtered_plot_data <- plot_data[plot_data$Num_pat >= 7 & plot_data$CI_high < 200 & plot_data$p_adj < 0.05, ]
    
    if (nrow(filtered_plot_data) > 0) {
      feature_label <- paste(filtered_plot_data$feature, " - ", filtered_plot_data$age_group)
      filtered_plot_data$p_adj <- ifelse(filtered_plot_data$p_adj <= 0.001, "< 0.001", sprintf("%.4f", round(filtered_plot_data$p_adj, 4)))
      
      tabletext <- cbind(
        c("Feature", filtered_plot_data$age_group),  
        c("Age group", as.character(filtered_plot_data$AGEGRP)),
        c("Patients with mutation", filtered_plot_data$Num_pat),
        c("Total patients", filtered_plot_data$Tot_pat), 
        c(" ",  rep("", nrow(filtered_plot_data))),
        c("OR", filtered_plot_data$HR), #update OR/HR based on the endpoint
        c("95% CI", filtered_plot_data$CI), 
        c("p-value", filtered_plot_data$pvalue),
        c("p-adjusted", filtered_plot_data$p_adj)
      )
      
      mean_values <- c(NA,  as.numeric(filtered_plot_data$HR))
      lower <- c(NA,  filtered_plot_data$CI_low)
      upper <- c(NA, filtered_plot_data$CI_high)
      
      # Determine box colors based on mean values
      box_colors <- ifelse(as.numeric(mean_values) < 0.25, color_palette[1],
                           ifelse(as.numeric(mean_values) >= 0.25 & as.numeric(mean_values) < 0.5, color_palette[2],
                                  ifelse(as.numeric(mean_values) >= 0.5 & as.numeric(mean_values) < 0.75, color_palette[3],
                                         ifelse(as.numeric(mean_values) >= 0.5 & as.numeric(mean_values) < 0.75, color_palette[4],
                                                ifelse(as.numeric(mean_values) >= 0.75 & as.numeric(mean_values) < 1, color_palette[5],
                                                       ifelse(as.numeric(mean_values) >= 1 & as.numeric(mean_values) < 2, color_palette[7],
                                                              ifelse(as.numeric(mean_values) >= 2 & as.numeric(mean_values) < 3, color_palette[8],
                                                                     ifelse(as.numeric(mean_values) >= 3 & as.numeric(mean_values) < 4, color_palette[9],
                                                                            ifelse(as.numeric(mean_values) >= 4 & as.numeric(mean_values) < 5, color_palette[10], color_palette[11])))))))))
      
      # Create forest plot
      p <- forestplot(
        labeltext = tabletext,
        clip = c(0, 8),
        graph.pos = 5,
        mean = mean_values,
        lower = lower,
        upper = upper,
        zero = 1,
        cex = 2,
        lineheight = "auto",
        boxsize = 0.04,
        colgap = unit(10, "mm"),
        xticks = c(0, 1, 2, 3, 4, 5, 6, 7,8),
        lwd.ci = 3.5,
        ci.vertices = TRUE,
        ci.vertices.height = 0.3
      ) |>
        fp_add_lines(h_2 = gpar(lwd = 1)) |> 
        fp_set_style(
          box = box_colors,
          line = box_colors,
          summary = gpar(fill = box_colors, clr = "black", zero = "gray50"),
          txt_gp = fpTxtGp(
            label = gpar(cex = 1.4),
            ticks = gpar(cex = 1.1),
            xlab = gpar(cex = 1.3),
            title = gpar(cex = 1.1)
          )
        ) |>
        fp_set_zebra_style("#EFEFEF")
      
      # Define output filename
      filename <- paste0(output_directory, ".tiff")
      
      # Save forest plot as TIFF
      tiff(filename, width = 18, height = 15, units = "in", res = 300) #45 for insigni
      print(p)
      dev.off()
      
      # Clear plot_data for the next iteration
      plot_data <- data.frame()
    }
  }
}

