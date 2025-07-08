# devtools::install_github(repo = 'cole-trapnell-lab/monocle3', ref = "d4a9a35")
# Specific commit used in performing all analyses

# Input path to github directory
setwd("/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-6th-TPM1-DSPKI/")
path_to_github ="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-6th-TPM1-DSPKI/code/sci-plex-master/"




# Set directory for sciPlex github 
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggpubr)
  library(Matrix)
  library(monocle3)
})


DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))


# Read in the CDS
cds <- readRDS("../../../output/combo/cds_umi300_750_nuclei44705.RDS")
rownames(cds) <- substr(rownames(rowData(cds)),21,nchar(rownames(rowData(cds))))
rowData(cds)$id <- rownames(cds)

##### go to line 212


############## Barnyard for sciPlex Screen  ############## 

# Get human and mouse genes
human.genes_names = rownames(rowData(cds))[grepl("^ENSG", rowData(cds)$id)]
mouse.genes_names = rownames(rowData(cds))[grepl("^ENSMUSG", rowData(cds)$id)]
colData(cds)$n_human_umi = Matrix::colSums(exprs(cds)[human.genes_names,])
colData(cds)$n_mouse_umi = Matrix::colSums(exprs(cds)[mouse.genes_names,])

pdf("../../../output/docs/supplemental_figure4_barnyard.pdf", width = 1.5, height = 1.5)
#tiff("../../../output/docs/supplemental_figure4_barnyard.tiff", units="in", width=1.5, height=1.5, res=300)
colData(cds) %>%
  as.data.frame() %>%
  filter(combo %in% c("Plate1_D6","Plate1_E6","Plate1_A7","Plate1_B7","Plate1_C7")) %>% 
  ggplot()  + 
  geom_point(aes(x = n_human_umi, y = n_mouse_umi, color = cell_type), 
             size = .5, stroke = 0) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "top",
        text = element_text(size = 3),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0) +
  xlab("Human UMIs") +
  ylab("Mouse UMIs") +
  
  scale_y_continuous(breaks = c(0,1000,2000),
                     limits = c(0,2000)) +
  scale_x_continuous(breaks = c(0,1000,2000),
                     limits = c(0,2000))

#ggsave("../../../output/docs/supplemental_figure4_barnyard.png",width = 1.5, height = 1.5)
dev.off()





############## Possible Hashes in sciPlex Screen  ############## 

possible_barcodes <- data.frame(
  group = c("Legal", "Illegal"),
  value = c(96*52, 52*768)
)

pdf("../../../data/GSM4150378-sciPlex3-large-screen/supplemental_figure4_possible_hashSpace.pdf", 
    width = .5, height = .5)
ggplot(possible_barcodes, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, color = "black", stat = "identity")+
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual("Legal",values = 
                      c("Illegal" ="#67795D", 
                        "Legal" = "#E48899")) +
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0) 
dev.off()




`%notin%` <- Negate(`%in%`)
# Make a plot with the legal hash pairings observed 
pdf("../../../output/docs/supplemental_figure4_observed_hashSpace.pdf", 
    width = 2, height = 1)
colData(cds) %>%
  as.data.frame() %>%
  filter(combo %notin% c("Plate1_A5","Plate1_A6","Plate1_A7","Plate1_A8","Plate1_A9","Plate1_A10","Plate1_A11","Plate1_A12")) %>%
  mutate(tmp = is.na(drug_dose)) %>%
  ggplot() +
  geom_bar(aes(x = tmp, fill = tmp ),width = 1, color = "black", position = "stack")+
  monocle:::monocle_theme_opts() +
  scale_fill_manual("Legal",values = 
                      c("TRUE" ="#67795D", 
                        "FALSE" = "#E48899")) +
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1))  +
  ylab("Cells") +
  xlab("")  +
  scale_x_discrete(labels = c("TRUE" = "Illegal" ,
                              "FALSE" = "Legal")) + 
  coord_flip()
dev.off()


############## Hash Enrichment Ratios and Number of Hash Reads in sciPlex Screen  ############## 


pdf("../../../output/docs/supplemental_figure4_histogram_hash_umis_P.pdf", 
    width = 1.5, height = 1.5)
colData(cds) %>%
  as.data.frame() %>%
  filter(!is.na(hash_umis_p)) %>%
  ggplot() + 
  geom_histogram(aes(x = log10(hash_umis_p))) + 
  monocle:::monocle_theme_opts() +
  geom_vline(xintercept = log10(5), color = "#FF0000", size = .5) +
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0) +
  xlab("log10(Plate Hash UMIs)") +
  ylab("Cells")
dev.off()


pdf("../../../output/docs/supplemental_figure4_histogram_hash_umis_W_umi100_200_nuclei30003.pdf", 
    width = 1.5, height = 1.5)
colData(cds) %>%
  as.data.frame() %>%
  filter(!is.na(hash_umis_well)) %>%
  ggplot() + 
  geom_histogram(bins = 30, aes(x = log10(hash_umis_well))) + 
  monocle:::monocle_theme_opts() +
  geom_vline(xintercept = log10(5), color = "#FF0000", size = .5) +
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0) +
  xlab("log10(Well Hash UMIs)") +
  ylab("Cells")
dev.off()


pdf("../../../data/GSM4150378-sciPlex3-large-screen/supplemental_figure4_Enrichment_Ratio_histogram_P.pdf", 
    width = 1.5, height = 1.5)
colData(cds) %>%
  as.data.frame() %>%
  filter(!is.na(top_to_second_best_ratio_P)) %>%
  ggplot() + 
  geom_histogram(aes(x = log10(top_to_second_best_ratio_P))) + 
  monocle:::monocle_theme_opts() +
  geom_vline(xintercept = log10(5), color = "#FF0000", size = .5) +
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0) +
  xlab("log10(Enrichment Ratio)") +
  ylab("Cells")
dev.off()

pdf("../../../output/docs/supplemental_figure4_Enrichment_Ratio_histogram_W.pdf", 
    width = 1.5, height = 1.5)
colData(cds) %>%
  as.data.frame() %>%
  filter(!is.na(top_to_second_best_ratio_W)) %>%
  ggplot() + 
  geom_histogram(bins=30,aes(x = log10(top_to_second_best_ratio_W))) + 
  monocle:::monocle_theme_opts() +
  geom_vline(xintercept = log10(5), color = "#FF0000", size = .5) +
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0) +
  xlab("log10(Enrichment Ratio)") +
  ylab("Cells") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5),
                     limits = c(0,5))
dev.off()



# Write out the cells that pass our enrichment ratio cutoffs
`%notin%` <- Negate(`%in%`)
validly_labeled_cells = 
  colData(cds) %>%
  as.data.frame() %>%
  #filter(!is.na(top_oligo_P), 
  #       !is.na(top_oligo_W)) %>%
  filter(hash_umis_well >= 5,
  #       hash_umis_well >= 5,
  #       top_to_second_best_ratio_P >= 5,
          top_to_second_best_ratio_W >= 5) %>%
  filter(combo %notin% c("Plate1_D6","Plate1_E6","Plate1_F6","Plate1_G6","Plate1_H6","Plate1_A7","Plate1_B7","Plate1_C7")) %>%
  pull(cell)

`%notin%` <- Negate(`%in%`)
validly_labeled_cells = 
  colData(cds) %>%
  as.data.frame() %>%
  #filter(!is.na(top_oligo_P), 
  #      !is.na(top_oligo_W)) %>%
  filter(hash_umis_well >= 5,
         #       hash_umis_well >= 5,
         #       top_to_second_best_ratio_P >= 5,
         top_to_second_best_ratio_W >= 5) %>%
  pull(cell)

##summary
#umi100_200_nuclei40006 now has 29471 valid nuclei
write.table(x = validly_labeled_cells, 
            file =  "../../../output/combo/valid_cells_umi100_200_nuclei105694.tsv",
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")


# Isolate singly labeled cells and drop all mouse genes 

human.genes = grepl("^ENSG", rowData(cds)$id)
cds = cds[human.genes,validly_labeled_cells]

print(paste(dim(cds)[1]," genes and ",dim(cds)[2], " cells after filtering"))
#"20633  genes and  105694  cells after filtering"

saveRDS(cds,"../../../output/combo/cds_validnuclei_umi300_750_nuclei39240.RDS")
############## UMIs per cell type and correlation between replicates  ############## 


pdf("../../../output/docs/umi100_200_nuclei30003/Notebook2_supplemental_figure5_Boxplot_cell_line.pdf", width = 1.5, height = 1.5)
`%notin%` <- Negate(`%in%`)
colData(cds) %>% 
as.data.frame() %>%
#filter(time_point == 24) %>%
#(combo %notin% c("Plate2_A1","Plate2_A2","Plate2_A3","Plate2_A4","Plate2_A5","Plate2_A6","Plate2_A7","Plate2_A8","Plate2_A9","Plate2_A10","Plate2_A11","Plate2_A12")) %>%
filter(cell_line %in% c("patient","patient.corrected")) %>%
ggplot() +
geom_boxplot(aes(x = cell_line, y = log10(n.umi), fill =  cell_line), 
              size = .25,outlier.size = 1, outlier.stroke = 0) +
#scale_fill_manual(values = c("#438FCD","#D6C939","#DB3C6A","darkred")) +
scale_fill_brewer(palette="Paired")+
monocle3:::monocle_theme_opts()  + 
theme(legend.position = "none", 
      text = element_text(size = 6),  
      legend.key.width = unit(0.2,"line"), 
      legend.key.height = unit(0.2,"line"),
      legend.text.align = 0) +
xlab("Cell Type") +
ylab("log10(RNA UMIs)") 
dev.off()

# Make a correlation plot between genes for different replicates
cells_24hrs = 
  colData(cds) %>% 
  as.data.frame() %>%
  filter(time_point == 24) %>%
  pull(cell)

cds_24hrs = cds[,cells_24hrs]
cds_24hrs = estimate_size_factors(cds_24hrs)
cds = estimate_size_factors(cds)  ## tried this line of code for 2nd experiment and changed "cds_24hrs" to "cds" for the rest of the following codes to make this figure.

counts_per_gene = 
  colData(cds) %>% 
  as.data.frame() %>%
  group_by(cell_line, replicate)%>%
  nest() %>%
  mutate(counts_per_gene = purrr::map(data, .f = function(pdata_subset, cds) {
    cds_subset = cds[,pdata_subset$cell]
    norm_exprs_mat = Matrix::t(Matrix::t(counts(cds_subset))/size_factors(cds_subset))
    totals = Matrix::rowSums(norm_exprs_mat)
    tibble(gene_id = rownames(exprs(cds_subset)),
           totals = totals )
  },cds))



counts_per_gene = 
  counts_per_gene %>%
  unnest(counts_per_gene)

counts_per_gene %>%
  spread(key = replicate, value = totals) %>% 
  ggplot() +
  geom_point(aes(x = log10(rep1 + 1), y = log10(rep2 + 1))
             , stroke = 0, size = .25, color = "grey80") +
  geom_abline(intercept = 0, slope = 1, size = .25) +
  geom_smooth(aes(x = log10(rep1 + 1), y = log10(rep2 + 1)),
              method = "lm", color = "red", size = .25 ) +
  stat_cor(aes(x = log10(rep1 + 1), y = log10(rep2 + 1)),
               size=2, label.y.npc=0.01, color="grey31", method = "pearson") +
  monocle:::monocle_theme_opts() +
  facet_wrap(~cell_type) +
  theme(legend.position = "none",
        text = element_text(size = 6),
        strip.text.x = element_text(size = 6)) +
  ggsave("supplemental_figure5_per_gene_replicate_correlation.png",
         height = 1.5, width = 4, dpi = 1200)

counts_per_gene_drug_dose = 
  colData(cds_24hrs) %>% 
  as.data.frame() %>%
  group_by(cell_type, replicate, catalog_number ,dose)%>%
  nest() %>%
  mutate(counts_per_gene = purrr::map(data, .f = function(pdata_subset, cds) {
    cds_subset = cds_24hrs[,pdata_subset$cell]
    norm_exprs_mat = Matrix::t(Matrix::t(counts(cds_subset))/size_factors(cds_subset))
    totals = Matrix::rowSums(norm_exprs_mat)
    tibble(gene_id = rownames(exprs(cds_subset)),
           totals = totals )
  },cds_24hrs)) %>%
  dplyr::select(-data) %>%
  unnest()

counts_per_gene_drug_dose =
  counts_per_gene_drug_dose %>% 
  spread(key = replicate, value = totals) 

per_drug_dose_correlations =
  counts_per_gene_drug_dose %>%
  drop_na() %>%
  group_by(cell_type, catalog_number,dose) %>%
  nest() %>% 
  mutate(correlation = purrr::map(data, .f = function(count_mat_long){
    cor(log10(count_mat_long$rep1 + 1) , log10(count_mat_long$rep2 + 1)) 
  })) %>%
  dplyr::select(-data) %>%
  unnest()

ggplot(per_drug_dose_correlations) +
  geom_violin(aes(x = as.factor(cell_type), y = correlation, fill =  cell_type), 
size = .25,outlier.size = 1, outlier.stroke = 0) +
  scale_fill_manual(values = c("#438FCD","#D6C939","#DB3C6A")) +
  monocle_theme_opts()  +
  theme(legend.position = "none", 
        text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        legend.text.align = 0) +
  xlab("Cell Type") +
  ylab("Replicate Correlation") +
  ggsave("supplemental_figure5_well_level_correlations between replicates.pdf",height = 3, width =3)


cor_matrix_A549 =
  counts_per_gene %>%
  filter(cell_type == "A549") %>% 
  spread(key = replicate, value = totals,fill = 0) 

cor_a549 = cor(x = log10(cor_matrix_A549$rep1 + 1), y = log10(cor_matrix_A549$rep2 + 1), method = "pearson"  )  

cor_matrix_K562 =
  counts_per_gene %>%
  filter(cell_type == "K562") %>% 
  spread(key = replicate, value = totals,fill = 0) 

cor_k562 = cor(x = log10(cor_matrix_K562$rep1 + 1), y = log10(cor_matrix_K562$rep2 + 1), method = "pearson"  )  

cor_matrix_MCF7 =
  counts_per_gene %>%
  filter(cell_type == "MCF7") %>% 
  spread(key = replicate, value = totals,fill = 0) 

cor_mcf7 = cor(x = log10(cor_matrix_MCF7$rep1 + 1), y = log10(cor_matrix_MCF7$rep2 + 1), method = "pearson"  )  

############## Vehicle Cells Recovered per replicate   ############## 

##"applied some changes to the following code."~Neda
pdf("../../../output/docs/Notebook2_supplemental_figure5_Boxplot_cellsRecovered.pdf", width = 4.5, height = 4)
`%notin%` <- Negate(`%in%`)
colData(cds) %>%
as.data.frame() %>%
#filter(vehicle) %>%
#group_by(Combo,replicate, cell_type) %>%
filter(combo %notin% c("Plate2_A1","Plate2_A2","Plate2_A3","Plate2_A4","Plate2_A5","Plate2_A6","Plate2_A7","Plate2_A8","Plate2_A9","Plate2_A10","Plate2_A11","Plate2_A12")) %>%
group_by(condition, combo) %>%
summarise(count = n()) %>%
ggplot(aes(x = condition, y = count, fill = condition)) +
geom_boxplot(size = 0.25, outlier.size = .25, outlier.stroke = 0) +
ylim(0,NA) +
theme(legend.position = "none", 
      text = element_text(size = 6), 
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
      axis.text.y = element_text(size = 6),
      legend.key.width = unit(0.2,"line"), 
      legend.key.height = unit(0.2,"line")) +
scale_fill_brewer(palette="Paired")+
#scale_fill_manual(values = c("patient" = "firebrick2",
#                             "patient_4F" = "navy",
#                             "patient_GFP" = "firebrick2",
#                             "patient_corrected" = "navy",
#                             "patient_corrected_4F" = "firebrick2",
#                             "patient_corrected_GFP" = "navy",
#                             "WT" = "firebrick2",
#                             "WT_4F" = "navy",
#                             "WT_GFP" = "firebrick2",
#                             "WT_mutation" = "navy",
#                             "WT_mutation_4F" = "firebrick2",
#                             "WT_mutation_GFP" = "navy"))+
xlab("Condition")+
ylab("Cells recovered\nper condition") +
#ylab("Cells recovered\nper vehicle well") +
monocle:::monocle_theme_opts()
#monocle:::monocle_theme_opts() +
#ggsave("supplemental_figure5_cells_recovered_by_vehicle_well.pdf", width = 4.5, height = 1.25)
dev.off()
  

colData(cds) %>%
  as.data.frame() %>%
  filter(combo %notin% c("Plate2_A1","Plate2_A2","Plate2_A3","Plate2_A4","Plate2_A5","Plate2_A6","Plate2_A7","Plate2_A8","Plate2_A9","Plate2_A10","Plate2_A11","Plate2_A12")) %>%
  group_by(condition, combo) %>%
  summarise(count = n()) %>%
  group_by(condition) %>%
  distinct() %>%
  summarise(median_n = median(count))

#### the following code for vehicle is the original one for the one above.
colData(cds_24hrs) %>%
  as.data.frame() %>%
  filter(vehicle) %>%
  group_by(Combo,replicate, cell_type) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = cell_type, y = count, fill = replicate)) +
  geom_boxplot(size = 0.25, outlier.size = .25, outlier.stroke = 0) +
  ylim(0,NA) +
  theme(legend.position = "right", 
        text = element_text(size = 6), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line")) +
  scale_fill_manual(values = c("rep1" = "firebrick2",
                                "rep2" = "navy")) +
  xlab("Replicate")+
  ylab("Cells recovered\nper vehicle well") +
  monocle:::monocle_theme_opts() +
  ggsave("supplemental_figure5_cells_recovered_by_vehicle_well.pdf", width = 4.5, height = 1.25)
  
  
colData(cds_24hrs) %>%
  as.data.frame() %>%
  filter(!vehicle) %>%
  group_by(Combo,replicate, cell_type) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = cell_type, y = count, fill = replicate)) +
  geom_boxplot(size = 0.25, outlier.size = .25, outlier.stroke = 0) +
  theme_void() +
  theme(legend.position = "none", 
        text = element_text(size = 6), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line"),
        panel.grid.major.x = element_line(colour = "grey50",linetype = 3)) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(breaks = c(0,100, 200,300,400),
                     limits = c(0,350))+
  xlab("")+
  ylab("Cells recovered") +
  coord_flip() +
  ggsave("figure3_cells_recovered_per_condition.pdf", width = 1, height = 1.25)

colData(cds_24hrs) %>%
  as.data.frame() %>%
  filter(!vehicle) %>%
  group_by(Combo,replicate, cell_type) %>%
  summarise(count = n()) %>%
  group_by(replicate, cell_type) %>%
  distinct() %>%
  summarise(median_n = median(count))
  
