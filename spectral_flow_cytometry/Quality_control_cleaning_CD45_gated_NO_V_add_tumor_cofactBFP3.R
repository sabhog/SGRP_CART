rm(list = ls())

library(flowWorkspace)
library(PeacoQC)
library(flowVS)
library(cowplot)
library(ggplot2)
library(tidyverse)

setwd("/Volumes/BTIT$/Sabrina/Tomas/Aurora_EGFRvIII_09_06_2023/Unmixed_Ndb_singlets_live_CD45_cleaned_NO_VEHICLE_add_tumor/cofactorBFP3/0_Quality_control_cleaning")

workingDir <- "/Volumes/BTIT$/Sabrina/Tomas/Aurora_EGFRvIII_09_06_2023/Unmixed_Ndb_singlets_live_CD45_cleaned_NO_VEHICLE_add_tumor/cofactorBFP3/0_Quality_control_cleaning"
prevDir <- "/Volumes/BTIT$/Sabrina/Tomas/Aurora_EGFRvIII_09_06_2023/Unmixed_Ndb_singlets_live_CD45_cleaned_NO_VEHICLE_add_tumor/0_Quality_control_cleaning"

data_folder <- "/Volumes/BTIT$/Sabrina/Tomas/Aurora_EGFRvIII_09_06_2023/data"

### Load data
panel <- read.csv(file.path(data_folder, "CD45_gated fcs-files/panel_mTagBFP2.csv"))

md <- read.csv(file.path(data_folder, "CD45_gated fcs-files/metadata.csv"))
md <- md[md$condition!="V",]

fcs.gated.files <- list.files(file.path(data_folder, "CD45_gated fcs-files"), pattern = '.fcs', full = F)
table(fcs.gated.files %in% md$file_name)
table(md$file_name %in% fcs.gated.files)

### create flowset
fs1 <- read.flowSet(md$full_file_name, transformation = FALSE, truncate_max_range = FALSE)

md <- dplyr::select(md, -full_file_name)

saveRDS(fs1, file.path(workingDir, "rds/flowSet_CD45_gated_NO_V.rds"))
# fs1 <- readRDS("rds/flowSet_CD45_gated_NO_V.rds")


### transform data (arcsinh)
features.keep <- panel$fcs_colname[!grepl("SC-|Zombie|AF", panel$fcs_colname)]

# remove FJComp-BFP-A because it doesn't run and doesn't make sense on CD45 gated cells
features.cofact <- panel$fcs_colname[!grepl("SC-|Zombie|AF|FJComp-BFP-A", panel$fcs_colname)]

#cofactors <- estParamFlowVS(fs1, channels=features.cofact)
## Warning messages:
##   1: In optimize(afun, b, maximum = TRUE) :
##   NA/Inf replaced by maximum positive value
# names(cofactors) <- features.cofact
# saveRDS(cofactors, file.path(workingDir, "rds/cofactors.rds"))
cofactors <- readRDS(file.path(prevDir, "rds/cofactors.rds"))

# ## visualise plots to find right cofactore for FJComp-BFP-A
x<-transFlowVS(fs1, channels="FJComp-BFP-A", 3)
densityplot(~`FJComp-BFP-A`, x)
# CytoExploreR::cyto_plot(x,
#                         parent = "root",
#                         channels = c("FJComp-BFP-A", "FSC-H"))
ggcyto::ggcyto(x, aes(x = `FJComp-BFP-A`, y =  `FSC-H`)) + geom_hex(bins = 128)


# add cofactor for FJComp-BFP-A because the function couldn't calculate it
cofactors2 <- c(cofactors[1:4], 3, cofactors[5:19])
names(cofactors2) <- c(names(cofactors)[1:4], "FJComp-BFP-A", names(cofactors)[5:19])
saveRDS(cofactors2, file.path(workingDir, "rds/cofactors_add_tumor_BFPcofact3.rds"))
# cofactors2 <- readRDS("rds/cofactors_add_tumor_BFPcofact3.rds")


fs1.VS <- transFlowVS(fs1, channels=names(cofactors2), cofactors2)
saveRDS(fs1.VS, file.path(workingDir, "rds/flowSet_CD45_gated_NO_V_arcsinh_transformed_add_tumor_BFPcofact3.rds"))

# fs1.VS <- readRDS("rds/flowSet_CD45_gated_NO_V_arcsinh_transformed_add_tumor_BFPcofact3.rds")

pdf(file.path(workingDir, "plots/arcsinh_transformation_add_tumor.pdf"))
densityplot(~`FJComp-APC-A`, fs1.VS, main="Siglec-H - Transfromed channels")
densityplot(~`FJComp-APC-Fire 810-A`, fs1.VS, main="P2ry12 - Transfromed channels")
densityplot(~`FJComp-Alexa Fluor 488-A`, fs1.VS, main="Axl - Transfromed channels")
densityplot(~`FJComp-Alexa Fluor 700-A`, fs1.VS, main="CD206 - Transfromed channels")
densityplot(~`FJComp-BFP-A`, fs1.VS, main="mTagBFP2 - Transfromed channels")
densityplot(~`FJComp-BUV496-A`, fs1.VS, main="CD4 - Transfromed channels")
densityplot(~`FJComp-BUV563-A`, fs1.VS, main="Ly6G - Transfromed channels")
densityplot(~`FJComp-BUV737-A`, fs1.VS, main="CD45 - Transfromed channels")
densityplot(~`FJComp-BV421-A`, fs1.VS, main="XCR1 - Transfromed channels")
densityplot(~`FJComp-BV510-A`, fs1.VS, main="CD11b - Transfromed channels")
densityplot(~`FJComp-BV605-A`, fs1.VS, main="F4/80 - Transfromed channels")
densityplot(~`FJComp-BV650-A`, fs1.VS, main="CD8a - Transfromed channels")
densityplot(~`FJComp-BV711-A`, fs1.VS, main="MerTK - Transfromed channels")
densityplot(~`FJComp-BV785-A`, fs1.VS, main="Ly6C - Transfromed channels")
densityplot(~`FJComp-NFBlue 610-70S-A`, fs1.VS, main="MHCII - Transfromed channels")
densityplot(~`FJComp-PE-A`, fs1.VS, main="TNFa - Transfromed channels")
densityplot(~`FJComp-PE-Cy5-A`, fs1.VS, main="CD86 - Transfromed channels")
densityplot(~`FJComp-PE-Cy7-A`, fs1.VS, main="CD11c - Transfromed channels")
densityplot(~`FJComp-PerCP-A`, fs1.VS, main="EGFRvIII - Transfromed channels")
densityplot(~`FJComp-Zombie NIR-A`, fs1.VS, main="live_dead - Transfromed channels")
densityplot(~`FJComp-mCherry-A`, fs1.VS, main="Tcells - Transfromed channels")
dev.off()

pdf(file.path(workingDir, "plots/comparison_transformation_add_tumor.pdf"))
lapply(names(cofactors2), function(x) {
  a1 <-ggplot(fs1, aes(x = get(x), group = name)) + 
    geom_density(alpha = 0.2) +
    xlim(c(-1e4,5e4))
  a2 <-ggplot(fs1.VS, aes(x = get(x), group = name)) + 
    geom_density(alpha = 0.2)+ 
    scale_x_continuous(name = paste0(x))
  plot_grid(a1, a2, labels = c('no transformation', 'asinh'))
})
dev.off()

pdf(file.path(workingDir, "plots/CD11b_scatter.pdf"), width=20, height=20)
CytoExploreR::cyto_plot(fs1,
                        parent = "root",
                        channels = c("FJComp-BV510-A", "FSC-H"),
                        title = "no transformation")
dev.off()

pdf(file.path(workingDir, "plots/CD11b_density.pdf"), width=10)
CytoExploreR::cyto_plot(fs1.VS,
                        parent = "root",
                        channels = c("FJComp-BV510-A"),
                        title = "arcsinh",
                        legend=T)
dev.off()



### QC

for (i in 1:length(fs1.VS)){
  peacoqc_res <- PeacoQC(
    ff=fs1.VS[[i]], 
    channels=features.cofact, 
    determine_good_cells="all", 
    save_fcs=TRUE, 
    plot=TRUE,
    output_directory = "plots",
    IT_limit = 0.55,
    MAD=6,
    suffix_fcs= paste0(substr(fs1.VS@phenoData@data$name[i], 1, nchar(fs1.VS@phenoData@data$name[i])-4)))
}
# Warning messages:
#   1: In FindIncreasingDecreasingChannels(breaks, ff, channels, plot,  :
#   There seems to be an increasing or decreasing trend in a channel for file1e241f7080d0 . 
#   Please inspect this in the overview figure.
#   2: In PeacoQC(ff = fs1.VS[[i]], channels = features.keep, determine_good_cells = "all",  :
#   There are not enough bins for a robust isolation tree analysis.

setwd("/Volumes/BTIT$/Sabrina/Tomas/Aurora_EGFRvIII_09_06_2023/Unmixed_Ndb_singlets_live_CD45_cleaned_NO_VEHICLE_add_tumor/cofactorBFP3/0_Quality_control_cleaning/plots/PeacoQC_results/fcs_files")

file.rename(list.files(".", full=TRUE),
            gsub(".*CD45_pregated_export_", "CD45_pregated_export_", list.files(".", full=TRUE)))

setwd("/Volumes/BTIT$/Sabrina/Tomas/Aurora_EGFRvIII_09_06_2023/Unmixed_Ndb_singlets_live_CD45_cleaned_NO_VEHICLE_add_tumor/cofactorBFP3/0_Quality_control_cleaning")


location <- "plots/PeacoQC_results/PeacoQC_report.txt"

# Make heatmap overview of the quality control run
pdf(file.path(workingDir, "plots/PeacoQC_results/heatmap_overview.pdf"))
PeacoQCHeatmap(report_location=location, show_values = FALSE,
               show_row_names = FALSE)
dev.off()

# PlotPeacoQC(fs1.VS[[1]], features.keep, display_peaks=TRUE, prefix = "PeacoQC_peaks_")
# PlotPeacoQC(fs1.VS[[1]], features.keep, display_peaks=TRUE, prefix = "PeacoQC_nopeaks_")

a <- lapply(c(1:10), function(x) {a<-fs1.VS[[x]]@description$GUID})
x <- cbind(fs1.VS@phenoData@data, Filename=unlist(a))

report <- read.delim(file.path(workingDir, "plots/PeacoQC_results/PeacoQC_report.txt"))
report <- report%>%
  left_join(x, by = "Filename")%>%
  relocate(name)


write.csv(report, file.path(workingDir, "plots/PeacoQC_results/PeacoQC_report.csv"), row.names = F)
