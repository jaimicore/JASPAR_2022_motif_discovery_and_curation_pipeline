#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "UniprotR")
for (lib in required.libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}


## ChIP-seq motifs
motif.folder <- "/storage/mathelierarea/processed/jamondra/Projects/JASPAR/JASPAR_2022/Zinc_fingers_C2H2_ChIP-seq/PFMs"
logo.folder <- "/storage/mathelierarea/processed/jamondra/Projects/JASPAR/JASPAR_2022/Zinc_fingers_C2H2_ChIP-seq/Logos"

## ChIP-exo motifs
motif.folder <- "/storage/mathelierarea/processed/jamondra/Projects/JASPAR/JASPAR_2022/Zinc_fingers_C2H2_ChIP-exo/PFMs"
logo.folder <- "/storage/mathelierarea/processed/jamondra/Projects/JASPAR/JASPAR_2022/Zinc_fingers_C2H2_ChIP-exo/Logos"


all.motifs.file <- file.path(motif.folder, "All_motifs_concatenated.tf")


uniprot.name.dict <- data.frame(Uniprot = list.files(motif.folder),
                                TF_name = rep(NA, times = length(list.files(motif.folder))),
                                File_ori = file.path(motif.folder, list.files(motif.folder)))


uniprot.name.dict$Uniprot <- gsub(uniprot.name.dict$Uniprot, pattern = ".pfm.txt", replacement = "")


uniprot.call <- ConvertID(uniprot.name.dict$Uniprot, ID_from = "ACC+ID", ID_to = "GENENAME")
uniprot.name.dict$TF_name <- as.vector(uniprot.call[,2])

## Export the Uniprot - TF name mapping table
fwrite(uniprot.name.dict, file = file.path(motif.folder, "Uniprot_TF_name_dict.txt"), sep = "\t")


for (i in seq_along(uniprot.name.dict$File_ori)) {
  
  ## Get file and Tf name
  f <- as.vector(uniprot.name.dict$File_ori[i])
  tf <- as.vector(uniprot.name.dict$TF_name[i])
  
  ## Set new motif file name
  new.file <- file.path(motif.folder, paste0(tf, "_pfm.tf"))
  
  ## RSAT onvert-matrix command
  ## From cis-bp to transfac format
  ## Add manually the TF name
  convert.matrix.cmd <- paste0("convert-matrix -i ", f, " -from cis-bp -to tf -multiply 1000 -attr name ", tf," -attr accession ", tf, " > ", new.file)
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

