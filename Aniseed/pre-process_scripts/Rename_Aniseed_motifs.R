#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "optparse",
                        "tidyr",
                        "UniprotR")
for (lib in required.libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}


####################
## Read arguments ##
####################
message("; Reading arguments from command-line")
option_list = list(
  make_option( c("-i", "--input_motif_folder"), type = "character", default = NULL, help = "Folder containing the Aniseed motifs. One file per motif.", metavar = "character"),
  make_option( c("-o", "--output_motif_folder"), type = "character", default = NULL, help = "Output folder where the Logos and PFMs folders will be created", metavar = "character"),
  make_option( c("-a", "--annotation_table"), type = "character", default = NULL, help = "Annotation table from Aniseed", metavar = "character"),
  make_option( c("-s", "--species"), type = "character", default = NULL, help = "Species", metavar = "character"),
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

input.motif.folder    <- opt$input_motif_folder
output.motif.folder   <- opt$output_motif_folder
motif.source          <- opt$annotation_table
species               <- opt$species



input.motif.folder       <- "/run/user/280010/gvfs/sftp:host=biotin3.hpc.uio.no,user=jamondra/storage/scratch/JASPAR_2022/Aniseed/data/Best_round/pfm_best_round"
output.motif.folder        <- "/run/user/280010/gvfs/sftp:host=biotin3.hpc.uio.no,user=jamondra/storage/scratch/JASPAR_2022/Aniseed/data/Best_round"
annotation.table   <- "/run/user/280010/gvfs/sftp:host=biotin3.hpc.uio.no,user=jamondra/storage/scratch/JASPAR_2022/Aniseed/data/Best_round/SELEX_project_TableS1_S3_Nitta_et_al_2019.csv"
organism           <- "Ciona_robusta"

## Create output directory
pfm.output.dir <- file.path(output.motif.folder, "PFMs")
dir.create(pfm.output.dir, showWarnings = F, recursive = T)


#############################################
## Read and parse Aniseed annotation table ##
#############################################
message("; Reading Aniseed annotation table")
aniseed.tab <- fread(annotation.table, header = T)

aniseed.tab <- aniseed.tab[,c(2,3,4,5,6,7,11:15)]
  
colnames(aniseed.tab) <- c("TF_fam",
                           "TF_class",
                           "Gene_name",
                           "Alt_gene_name",
                           "DBD",
                           "SELEX_clone",
                           "InterproID_1",
                           "InterproID_2",
                           "InterproID_3",
                           "InterproID_4",
                           "InterproID_5")

## PArse the Aniseed table
##
## Keep entries with SELEX_clone ID
## Transform the interproID colums to one entry per line
aniseed.tab <- aniseed.tab %>% 
                dplyr::filter(SELEX_clone != "") %>% 
                mutate(InterproID = paste(InterproID_1, InterproID_2, InterproID_3, InterproID_4, InterproID_5, sep = ",")) %>% 
                transform(InterproID = gsub(x = InterproID, pattern = ",{2,}", replacement = "")) %>% 
                select(-c(InterproID_1, InterproID_2, InterproID_3, InterproID_4, InterproID_5)) %>% 
                separate_rows(InterproID, sep = ",")


## New file with all motifs
all.motifs.file <- file.path(motif.folder, paste0("All_", organism, "_motifs_concatenated.tf"))

TF.name.counter <- list()
tmp.files.list <- NULL

## Iterate over each motif file
for (f in list.files(input.motif.folder)) {
  
  ## Get Selex probe ID
  selex.id <- gsub(f, pattern = "^([A-Za-z0-9]+)_\\d+.+$", replacement = "\\1")
  tf <- as.vector(subset(aniseed.tab, SELEX_clone == selex.id)$Gene_name)
  
  ## Count the number of times this TF appears in the file list
  if (is.null(TF.name.counter[[tf]])) {
    TF.name.counter[[tf]] <- 1
  } else {
    TF.name.counter[[tf]] <- TF.name.counter[[tf]] + 1
  }
  

  message("; Processing ", tf, " motif")
  
  ## Set new motif file name
  new.tmp.file <- file.path(pfm.output.dir, paste0(tf, "_pfm.jaspar"))
  tmp.files.list <- append(tmp.files.list, new.tmp.file)
  
  ## Add a header to the motifs
  add.header.cmd <- paste0( '(echo ">',tf,'" && cat ', file.path(input.motif.folder, f), ') > ', new.tmp.file )
  message("; ", add.header.cmd)
  system(add.header.cmd)
  
  
  ## RSAT convert-matrix command
  ## From cis-bp to transfac format
  ## Add manually the TF name
  convert.matrix.cmd <- paste0("convert-matrix -i ", file.path(input.motif.folder, f), " -from cis-bp -to tf -multiply 1000 -attr name ", tf," -attr accession ", tf, " > ", new.file)
  message("; ", convert.matrix.cmd)
  system(convert.matrix.cmd)
  
  ## Concat all motifs in a single file
  cat.all.motifs.cmd <- paste0("cat ", new.file, " >> ", all.motifs.file)
  system(cat.all.motifs.cmd)
  
  ## Remove the txt (cis-bp format) file
  rm.cisbp.file <- paste0("rm -rf ", f)
  system(rm.cisbp.file)
}


## Generate logos
logos.cmd <- paste("convert-matrix -i ", all.motifs.file, " -from tf -to tf -return logo -logo_dir ", logo.folder)
message("; ", logos.cmd)
system(logos.cmd)

