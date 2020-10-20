#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "optparse",
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
  make_option( c("-m", "--motif_folder"), type = "character", default = NULL, help = "Folder containing the C2H2 motifs in cis-bp format. One file per motif.", metavar = "character"),
  make_option( c("-l", "--logo_folder"), type = "character", default = NULL, help = "Output folder for the motif logos", metavar = "character"),
  make_option( c("-t", "--motif_source"), type = "character", default = NULL, help = "Experimental type where the motif are derived from: ChIP-seq or ChIP-exo", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

motif.folder   <- opt$motif_folder
logo.folder    <- opt$logo_folder
motif.source   <- opt$motif_source



## ChIP-seq motifs
# motif.folder <- "/storage/mathelierarea/processed/jamondra/Projects/JASPAR/JASPAR_2022/Zinc_fingers_C2H2_ChIP-seq/PFMs"
# logo.folder <- "/storage/mathelierarea/processed/jamondra/Projects/JASPAR/JASPAR_2022/Zinc_fingers_C2H2_ChIP-seq/Logos"
# motif.source <- "ChIP-seq"
# 
# ## ChIP-exo motifs
# motif.folder <- "/storage/mathelierarea/processed/jamondra/Projects/JASPAR/JASPAR_2022/Zinc_fingers_C2H2_ChIP-exo/PFMs"
# logo.folder <- "/storage/mathelierarea/processed/jamondra/Projects/JASPAR/JASPAR_2022/Zinc_fingers_C2H2_ChIP-exo/Logos"
# motif.source <- "ChIP-exo"


## New file with all motifs
all.motifs.file <- file.path(motif.folder, paste0("All_", motif.source, "_motifs_concatenated.tf"))


## Initialize dataframe
uniprot.name.dict <- data.frame(Uniprot = list.files(motif.folder),
                                TF_name = rep(NA, times = length(list.files(motif.folder))),
                                File_ori = file.path(motif.folder, list.files(motif.folder)))


## Get uniprot ID from motif file name
uniprot.name.dict$Uniprot <- gsub(uniprot.name.dict$Uniprot, pattern = ".pfm.txt", replacement = "")

## Map uniprot ID to TF name
message("; Mapping uniprot ID")
uniprot.call <- ConvertID(uniprot.name.dict$Uniprot, ID_from = "ACC+ID", ID_to = "GENENAME")
uniprot.name.dict$TF_name <- as.vector(uniprot.call[,2])

## Export the Uniprot - TF name mapping table
message("; Export uniprot ID <-> TF name dictionary")
fwrite(uniprot.name.dict, file = file.path(motif.folder, "Uniprot_TF_name_dict.txt"), sep = "\t")


## Iterate over each motif file
for (i in seq_along(uniprot.name.dict$File_ori)) {
  
  ## Get file and Tf name
  f <- as.vector(uniprot.name.dict$File_ori[i])
  tf <- as.vector(uniprot.name.dict$TF_name[i])
  
  message("; Processing ", tf, " motif")
  
  ## Set new motif file name
  new.file <- file.path(motif.folder, paste0(tf, "_pfm.tf"))
  
  ## RSAT onvert-matrix command
  ## From cis-bp to transfac format
  ## Add manually the TF name
  convert.matrix.cmd <- paste0("convert-matrix -i ", f, " -from cis-bp -to jaspar -multiply 1000 -attr name ", tf," -attr accession ", tf, " > ", new.file)
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
logos.cmd <- paste("convert-matrix -i ", all.motifs.file, " -from jaspar -to jaspar -return logo -logo_dir ", logo.folder)
message("; ", logos.cmd)
system(logos.cmd)

