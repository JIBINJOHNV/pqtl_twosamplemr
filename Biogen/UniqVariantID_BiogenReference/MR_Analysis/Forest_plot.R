.libPaths("/edgehpc/dept/human_genetics/jrobins_legacy/r_libs_4.1.0")
library(MendelianRandomization, lib.loc = '/home/jjohn1/modulefiles/mrpipeline')

suppressWarnings(suppressPackageStartupMessages({
    library(tidyr)
    library(glue)
    library(dplyr)
    library(plyr)
    library(data.table)
    library(magrittr)
    library(ggplot2)
}))


mr_forest_plot_custom <- function(singlesnp_results, exponentiate = FALSE) {
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("plyr", quietly = TRUE)
    res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"),
        function(d) {
            d <- plyr::mutate(d)
            if (sum(!grepl("All", d$SNP)) < 2) {
                return(blank_plot("Insufficient number of SNPs"))
            }
            levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
            levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
            am <- grep("All", d$SNP, value = TRUE)
            d$up <- d$b + 1.96 * d$se
            d$lo <- d$b - 1.96 * d$se
            d$tot <- 0.01
            d$tot[d$SNP %in% am] <- 1
            d$SNP <- as.character(d$SNP)
            nom <- d$SNP[!d$SNP %in% am]
            nom <- nom[order(d$b)]
            d <- rbind(d, d[nrow(d), ])
            d$SNP[nrow(d) - 1] <- ""
            d$b[nrow(d) - 1] <- NA
            d$up[nrow(d) - 1] <- NA
            d$lo[nrow(d) - 1] <- NA
            d$SNP <- factor(d$SNP, levels = unique(d$SNP))  # Keep original order
            xint <- 0
            if (exponentiate) {
                d$b <- exp(d$b)
                d$up <- exp(d$up)
                d$lo <- exp(d$lo)
                xint <- 1
            }

            ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) +
                ggplot2::geom_vline(xintercept = xint, linetype = "dotted") +
                ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo,
                  xmax = up, size = as.factor(tot), colour = as.factor(tot)),
                  height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) +
                ggplot2::geom_hline(yintercept = which(levels(d$SNP) %in% ""), colour = "grey") +
                ggplot2::scale_colour_manual(values = c("black", "red")) +
                ggplot2::scale_size_manual(values = c(0.3, 1)) +
                ggplot2::theme(legend.position = "none",
                axis.text.y = ggplot2::element_text(size = 8),
                axis.ticks.y = ggplot2::element_line(size = 0),
                axis.title.x = ggplot2::element_text(size = 8)) +
                ggplot2::labs(y = "", x = paste0("MR effect size for\n'",
                  d$exposure[1], "' on '", d$outcome[1], "'"))
        }
    )
    res
}



prefix="trans_exposure"

twosamplemr_single_file<-glue("{prefix}_TwoSampleMR_Analysis_SingleVariantWald_Test.csv")
twosamplemr_alltest_file<-glue("{prefix}_TwoSampleMR_Analysis_Multiple_MR_Test.csv")
biogen<-glue("{prefix}_MendelianPipelineTest.csv")
british_ivd_delta<-glue("{prefix}_MendelianRandomization_IVW_Delta_Test.csv")
british_all<-glue("{prefix}_MendelianRandomization_AllTest.csv")


##Two sample MR
two_samplemr_single_df <- fread(twosamplemr_single_file)
two_samplemr_single2_df<-two_samplemr_single_df[grepl("^rs", two_samplemr_single_df$SNP), ]
two_samplemr_single2_df <- two_samplemr_single2_df[, c("exposure", "outcome", "id.exposure", "id.outcome", "SNP", "b", "se")]

two_samplemr_alltest_df <- fread(twosamplemr_alltest_file)
two_samplemr_alltest_df<-two_samplemr_alltest_df[,c("exposure","outcome","id.exposure","id.outcome","method","b","se")]
colnames(two_samplemr_alltest_df)<-c("exposure", "outcome", "id.exposure", "id.outcome", "SNP", "b", "se")
two_samplemr_alltest_df$SNP <- gsub("Inverse variance weighted", "IVW", two_samplemr_alltest_df$SNP)
two_samplemr_alltest_df$SNP <- gsub("multiplicative random effects", "MRE", two_samplemr_alltest_df$SNP)
two_samplemr_alltest_df <- two_samplemr_alltest_df[!grepl("Unweighted", SNP)]
two_samplemr_alltest_df$SNP <- paste("TwoSampleMR-", two_samplemr_alltest_df$SNP, sep = "")


##biogen
biogen_df<-fread(biogen)
biogen_df$method[grepl("Wald", biogen_df$method)] <- biogen_df$snp[grepl("Wald", biogen_df$method)]
biogen_df<-biogen_df[,c("exposure","outcome","id.exposure","id.outcome","method","b","se")]
colnames(biogen_df)<-c("exposure", "outcome", "id.exposure", "id.outcome", "SNP", "b", "se")
biogen_df$SNP <- gsub("Inverse variance weighted", "IVW", biogen_df$SNP)
biogen_df$SNP <- gsub("Wald ratio", "WaldRatio", biogen_df$SNP)
biogen_df$SNP <- paste("Biogen-", biogen_df$SNP, sep = "")


#british_ivd_delta
british_ivd_delta_df<-fread(british_ivd_delta)
british_ivd_delta_df<-british_ivd_delta_df[,c("Model","Estimate","StdError","Exposure")]
colnames(british_ivd_delta_df) <- c("SNP", "b", "se", "exposure")
british_ivd_delta_df$SNP <- paste("British-", british_ivd_delta_df$SNP, sep = "")
two_samplemr_alltest_df_unique<-unique(two_samplemr_alltest_df[,c("exposure","outcome","id.exposure","id.outcome")])
british_ivd_delta_df <- merge(x = british_ivd_delta_df,y = two_samplemr_alltest_df_unique,by = c("exposure"),all.x = TRUE)

#british all methods
british_all_df<-fread(british_all)
british_all_df<-british_all_df[,c("exposure","outcome","Method","Estimate","Std Error")]
colnames(british_all_df)<-c("exposure","outcome","SNP", "b", "se")
british_all_df$SNP <- paste("BritishAll-", british_all_df$SNP, sep = "")
british_all_df$SNP <- gsub("Penalized", "P", british_all_df$SNP)

#allmr_methods<-bind_rows(two_samplemr_single2_df,biogen_df,british_ivd_delta_df,two_samplemr_alltest_df,british_all_df)
allmr_methods<-bind_rows(two_samplemr_single2_df,biogen_df,british_ivd_delta_df,two_samplemr_alltest_df)





Exposures<-unique(allmr_methods$exposure)
keywords<-C("AIF1","CCDC134","FES","PDCD5","PRDX6","TRDMT1","TNFRSF6B")
matching_exposures <- Exposures[grepl(paste(keywords, collapse = "|"), Exposures)]

matching_exposures<-c("AIF1:P55008:Oncology","CCDC134:Q9H6E4:Oncology_II","FES:P07332:Oncology","PDCD5:O14737:Neurology","PRDX6:P30041:Oncology",
                     "TRDMT1:O14717:Oncology_II","TNFRSF6B:O95407:Neurology")



for (Exposure in matching_exposures){

all_selected<-allmr_methods[grepl(Exposure, exposure)]
print(Exposure)
p2 <- mr_forest_plot_custom(all_selected)
modified_exposure <- gsub(":", "_", Exposure)
filename <- glue("{prefix}_{modified_exposure}_ForestPlot_TwosampleMR.png")
ggsave(p2[[1]], file = filename, width = 7, height = 7)

}
