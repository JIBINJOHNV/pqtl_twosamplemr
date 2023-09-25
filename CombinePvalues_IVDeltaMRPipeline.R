library(metap)
library(plyr)
library(glue)
library(ggplot2)
library(dplyr)


gwasnames <- c("PGC_ADHD2022_iPSYCH_deCODE","BIP_PGC3_noukb","ASD_PGC","PGC3_SCZ","Depression_iPSYCH_2023","PGC_AN2","PGC3_SCZ_NoUKB")
prefixes <- c("CisExposure_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta",
               "TransExposureNoMHCUnique_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta", 
              "TransExposureNoMHC_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta", 
              "TransExposure_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta")

pqtltype1 <- "Decode"
pqtltype2 <- "Biogen"
exclude_columns <- c("Gene_Symbol", "outcome")

for (gwasname in gwasnames) {
  for (prefix in prefixes) {
    ## Decode PQTL
    decode <- read.csv(glue("{pqtltype1}_{gwasname}_{prefix}.csv"))
    col_names <- c("exposure")
    new_col_names <- ifelse(columns_to_rename != "", paste0("Decode.PQTL_", columns_to_rename), columns_to_rename)
    colnames(decode)[colnames(decode) == "exposure"] <- new_col_names


    ## Biogen PQTL
    biogen <- read.csv(glue("{pqtltype2}_{gwasname}_{prefix}.csv"))
    col_names <- colnames(biogen)
    new_col_names <- ifelse(columns_to_rename != "", paste0("Biogen.PQTL_", columns_to_rename), columns_to_rename)
    colnames(biogen)[colnames(biogen) == "exposure"] <- new_col_names

    ## Merge biogen and decode data frames
    decode_biogen <- merge(decode, biogen, by = c("Gene_Symbol", "outcome"), all = TRUE)

    # Rest of your code for meta-analysis, reordering columns, and writing to CSV

    #MRIVWtest_Pvalue
    MRIVWtest_Pvalue_cols<-colnames(decode_biogen[, grep("IVWDelta_Pvalue$", colnames(decode_biogen))])
    MRIVWtest_miss<-decode_biogen[!complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]
    MRIVWtest_nomiss<-decode_biogen[complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]

    MRIVWtest_df <- MRIVWtest_nomiss[, grep("IVWDelta_Pvalue$", colnames(decode_biogen))]
    MRIVWtest_df <- MRIVWtest_df %>%mutate_at(MRIVWtest_Pvalue_cols, as.numeric)

    metap_MRIVWtest<-vector()

    for(indx in 1:nrow(MRIVWtest_df)) {
        teest<-unlist(as.vector(MRIVWtest_df[indx,]))
        metap_MRIVWtest<-c(metap_MRIVWtest,allmetap(teest,method="sumlog")$p[[1]] )}

    MRIVWtest_nomiss$metap_MRIVWtest<-metap_MRIVWtest
    MRIVWtest_miss$metap_MRIVWtest<-"NA"
    decode_biogen <- rbind(MRIVWtest_nomiss,MRIVWtest_miss)

    ##bland_altman_plot
    MRIVWtest_df[[colnames(MRIVWtest_df)[1]]] <- -log10(MRIVWtest_df[[colnames(MRIVWtest_df)[1]]])
    MRIVWtest_df[[colnames(MRIVWtest_df)[2]]] <- -log10(MRIVWtest_df[[colnames(MRIVWtest_df)[2]]])

    # Calculate the difference and mean between the two methods
    DecodePqtl_col <- colnames(MRIVWtest_df)[grep("Decode.PQTL_", colnames(MRIVWtest_df))]
    BiogenPqtl_col <- colnames(MRIVWtest_df)[grep("Biogen.PQTL_", colnames(MRIVWtest_df))]
    MRIVWtest_df$Difference <- MRIVWtest_df[[DecodePqtl_col]] - MRIVWtest_df[[BiogenPqtl_col]]
    MRIVWtest_df$Mean <- (MRIVWtest_df[[DecodePqtl_col]] + MRIVWtest_df[[BiogenPqtl_col]]) / 2

    # Calculate mean difference, standard deviation of differences, and limits of agreement
    mean_diff <- mean(MRIVWtest_df$Difference, na.rm = TRUE)
    std_diff <- sd(MRIVWtest_df$Difference, na.rm = TRUE)
    upper_limit <- mean_diff + 1.96 * std_diff
    lower_limit <- mean_diff - 1.96 * std_diff

    # Calculate coefficient of repeatability (CR)
    cr <- 1.96 * std_diff

    # Create a Bland-Altman plot
    bland_altman_plot <- ggplot(MRIVWtest_df, aes(x = Mean, y = Difference)) +
    geom_point(color = 'blue', size = 3) +
    geom_hline(yintercept = mean_diff, linetype = 'dashed', color = 'black') +
    geom_hline(yintercept = upper_limit, linetype = 'dashed', color = 'red') +
    geom_hline(yintercept = lower_limit, linetype = 'dashed', color = 'red') +
    labs(x = 'Mean of -log10(Decode_pval) and -log10(Biogen_pval)',
        y = 'Difference (-log10(Decode_pval) - -log10(Biogen_pval))',
        title = 'Bland-Altman Plot\nDifference and Limits of Agreement') +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
            plot.title = element_text(size = 16, hjust = 0.5),
            axis.title = element_text(size = 14, margin = margin(t = 20)),
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"),
            legend.position = "left",
            plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm"))  +
    annotate("text", x = -Inf, y = 1:5, label = c(
        paste('MeanDiff:', round(mean_diff, 4)),
        paste('StdDevDiff:', round(std_diff, 4)),
        paste('UpperLimit:', round(upper_limit, 4)),
        paste('LowerLimit:', round(lower_limit, 4)),
        paste('CR:', round(cr, 4))
    ), hjust = 0, vjust = 0.01, size = 4)

    # Save the Bland-Altman plot as a TIFF image
    outname2<-substr(prefix, 1, nchar(prefix) - 55)
    output_folder <- "bland_altman_plot"
    if (!file.exists(output_folder)) {dir.create(output_folder)}
    output_file <- file.path(output_folder, glue("{gwasname}_{outname2}_British_MRIVWtest_bland_altman_plot.tiff"))
    ggsave(output_file, plot = bland_altman_plot, dpi = 300, width = 10, height = 6)

    # Create a scatter plot using ggplot2
    MRIVWtest_correlation_plot <- ggplot(MRIVWtest_df, aes(x = MRIVWtest_df[[DecodePqtl_col]], y = MRIVWtest_df[[BiogenPqtl_col]])) +
    geom_point(color = 'blue') +
    labs(x = '-log10(Decode_MR_Pipeline_pval)',y = '-log10(Biogen_MR_Pipeline_pval)',title = 'Correlation Plot') +
    annotate("text", x =3 , y = 4,label = paste('Correlation:', round(cor(MRIVWtest_df[[DecodePqtl_col]], MRIVWtest_df[[BiogenPqtl_col]]), 4)))

    # Save the correlation plot as a TIFF image
    outname <- substr(prefix, 1, nchar(prefix) - 55)
    output_folder <- "correlation_plot"
    if (!file.exists(output_folder)) {dir.create(output_folder)}
    output_file <- file.path(output_folder, glue("{gwasname}_{outname2}_British_MRIVWtest_correlation_plot.tiff"))
    ggsave(output_file, plot = MRIVWtest_correlation_plot, dpi = 300, width = 10, height = 6)

       ##Reordering the columns
    first_cols <- colnames(decode_biogen[, grep("IVWDelta_Pvalue$|Biogen_pval$", colnames(decode_biogen))])
    second_cols <- colnames(decode_biogen[, grep("outcome$|exposure$", colnames(decode_biogen))])
    
    new_order <- c(
      c("Gene_Symbol"), first_cols, c( "metap_MRIVWtest"), second_cols,
      setdiff(colnames(decode_biogen), c("Gene_Symbol", first_cols, "metap_MRIVWtest", second_cols)))

    # Define the folder name
    output_folder <- "MetaP_AnalysisResults"
    if (!file.exists(output_folder)) {dir.create(output_folder)}
    output_file <- file.path(output_folder, glue("Biogen_Decode_pQTL_{gwasname}_{out_name}_MetapAnalysis.csv"))
    write.csv(decode_biogen_reordered_df, output_file, row.names = FALSE)

  }}
