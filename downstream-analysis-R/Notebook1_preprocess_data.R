path_to_github ="~/master/"

# Set directory for sciPlex github 
bin_directory = paste0(path_to_github,
                       "bin",
                       sep = "")

suppressPackageStartupMessages({  #Disable messages upon loading a package 
  library(dplyr)
  library(monocle3)
})

#Wrapping an array-like object (typically an on-disk object) in a DelayedArray object allows one to perform common array operations on it without loading the object in memory. In order to reduce memory usage and optimize performance, operations on the object are either delayed or executed using a block processing mechanism.
# First, let's see block processing in action:
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)

setwd(paste0(path_to_github,
             "large_screen",
             sep = ""))


# Read in the CDS
cds = readRDS("../../../output/combo/cds.RDS")
dim(cds)
#[1]43214 1298096
hash = readRDS("../../../output/combo/hash_df_all.RDS")
dim(hash)
#[1]117084      6
hash <- hash %>% filter(hash$qval<0.05)
dim(hash)
#[1] 115330      6
index <- match(hash$Cell,colnames(cds))
cds <- cds[,index]
dim(cds)
#[1] 43214 115330
################## Append Cell Path Metadata To CDS ################## 
#Pull out CDS column data to append relevant metadata
library(stringr)
coldata_cds = colData(cds) %>% as.data.frame()
coldata_cds$RTwell <- paste0("RT_BC_",str_split(coldata_cds$Cell,pattern = "_", simplify = T)[,5])
coldata_cds$Ligwell <- paste0("Lig_BC_",str_split(coldata_cds$Cell,pattern = "_", simplify = T)[,8])
coldata_cds$PCR_well <- substr(coldata_cds$Cell,1,7)
coldata_cds$PCR_plate <- rep(1,time=length(coldata_cds$Cell))
coldata_cds = 
  coldata_cds %>%
  dplyr::rename(cell = Cell,
                rt_well = RTwell,
                lig_well = Ligwell,
                pcr_well = PCR_well,
                pcr_plate = PCR_plate) %>%  
                dplyr::mutate(combo = hash$top_oligo)%>%
                dplyr::mutate(hash_umis_well = hash$hash_umis) %>%
                dplyr::mutate(top_to_second_best_ratio_W = hash$top_to_second_best_ratio)

RT_data <- read.delim("../process_from_raw/sci-RNA-seq-pipeline-scripts/RT_384_with_plate_3lvl",stringsAsFactors = FALSE, sep="", header = T, row.names = 1)
lig_data <- read.delim("../process_from_raw/sci-RNA-seq-pipeline-scripts/ligation_384_with_plate_3lvl",stringsAsFactors = FALSE, sep="", header = T, row.names = 1)

#### <no need to run the following codes> this is part of the original code
#rt_plate = 
#  data.frame(rt_well = seq(1,384),
#             rt_plate = seq(0,383) %/% 96  + 1)
#lig_plate = 
#  data.frame(lig_well = seq(1,384),
#             lig_plate = seq(0,383) %/% 96  + 1)

#coldata_cds = left_join(coldata_cds, rt_plate, by = "rt_well")
#coldata_cds = left_join(coldata_cds, lig_plate, by = "lig_well")
#coldata_cds$Combo = paste0(coldata_cds$top_oligo_W,coldata_cds$top_oligo_P) 
####</end>

################## Append Hash Encoding Metadata To CDS ################## 

#### <no need to run the following codes> this is part of the original code
#hash_meta_data = read.table("../large_screen/bin/sciChem3_screen2_hashMeta_data.tsv",header = F, sep = "\t")
#colnames(hash_meta_data) = 
#  c("combo","plate_oligo","cell_type","replicate","time_point","drug_dose")
####</end>

hash_meta_data = read.delim("../large_screen/bin/hash_meta_data_3lvl_Enzo_Porrello",header = F, sep = "\t",stringsAsFactors = FALSE)
colnames(hash_meta_data) = 
  c("combo","plate_oligo","cell_line","replicate","density","seed","media")
coldata_cds = left_join(coldata_cds, hash_meta_data, by = "combo")

#### <no need to run the following codes> this is part of the original code
#hash_meta_data$Combo =  paste0(hash_meta_data$well_oligo,hash_meta_data$plate_oligo)
#coldata_cds = left_join(coldata_cds, hash_meta_data, by = "Combo")
####</end>

coldata_cds = left_join(coldata_cds, hash_meta_data, by = "combo")

coldata_cds$catalog_number = stringr::str_split_fixed(coldata_cds$drug_dose, pattern = "_",n = 2)[,1]

coldata_cds$vehicle = grepl(pattern = "S000",x = coldata_cds$catalog_number)

coldata_cds$dose_pattern = stringr::str_split_fixed(coldata_cds$drug_dose, pattern = "_",n = 2)[,2]

dose_keys = data.frame(dose_pattern = c("1","2","3","4"),
           dose_character = c("10000","1000","100","10"))

coldata_cds = left_join(coldata_cds,dose_keys, by = "dose_pattern")

# Set the dose of vehicle cells to 0
coldata_cds$dose_character = ifelse(coldata_cds$vehicle, "0", as.character(coldata_cds$dose_character))
coldata_cds$dose = as.numeric(coldata_cds$dose_character)

coldata_cds$catalog_number = ifelse(coldata_cds$vehicle, "S0000", as.character(coldata_cds$catalog_number))
coldata_cds$treatment = as.character(coldata_cds$catalog_number)

################## Append Drug Pathway Metadata To CDS ################## 

pathway_annotations = read.table("../large_screen/bin/Supplementary_Table_3.txt",
                                 header = T, sep = "\t")

pathway_annotations = 
  pathway_annotations %>%
  dplyr::select(pathway_level_1,
                pathway_level_2,
                catalog_number) %>%
  dplyr::mutate(catalog_number  = as.character(catalog_number)) %>%
  dplyr::distinct()


pathway_annotations$catalog_number[is.na(pathway_annotations$catalog_number)] <- "S0000"

coldata_cds = left_join(coldata_cds,pathway_annotations, by ="catalog_number" )

################## Append Drug Metadata To CDS ################## 

drug_annotations = read.table("../large_screen/bin/drugProperties_short.tsv",
                              header = T, 
                              sep = "\t")

drug_annotations =
  drug_annotations %>%
  dplyr::rename(catalog_number = Catalog.Number,
                product_name = Product.Name,
                pathway = Pathway,
                target = Target)


coldata_cds = left_join(coldata_cds,drug_annotations, by = "catalog_number")


################## Write and save new CDS ################## 

rownames(coldata_cds) = coldata_cds$cell
mat = counts(cds)
gene_metadata = rowData(cds)

cds = new_cell_data_set(mat,
                        cell_metadata =coldata_cds,
                        gene_metadata = gene_metadata)

saveRDS(object = cds, 
        file = "../../../output/combo/cds_umi300_750_nuclei18113.RDS")

