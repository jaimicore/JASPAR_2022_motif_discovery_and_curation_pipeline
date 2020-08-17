#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "optparse")
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
  make_option( c("-i", "--peaks_in_folder"), type = "character", default = NULL, help = "Folder containing the C2H2 motifs in cis-bp format. One file per motif.", metavar = "character"),
  make_option( c("-o", "--peaks_out_folder"), type = "character", default = NULL, help = "Output folder for the motif logos", metavar = "character"),
  make_option( c("-s", "--suffix"), type = "character", default = NULL, help = "Experimental type where the motif are derived from: ChIP-seq or ChIP-exo", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

in.folder    <- opt$peaks_in_folder
out.folder   <- opt$peaks_out_folder
suffix       <- opt$suffix


## Debug
# in.folder <- "/home/jamondra/Downloads/Yeast_peaks/data/yep-peaks"
# out.folder <- "/home/jamondra/Downloads/Yeast_peaks/data/Yeast_SacCer3"
# suffix <- "_peaks.narrowPeak"


## The ChIP-seq does not have an ID or an associated condition
## We will use a generic ID: ChExMix
id <- "ChExMix"


## Conversion chromosome IDs
## From decimal to roman numbers
## S cereviseae has 16 chromosomes + M
chrom <- c(1:16, "M")
sacCer3.chrom <- data.frame(Old = chrom,
                            New = as.character(as.roman(chrom)))
sacCer3.chrom$Old <- paste0("chr", sacCer3.chrom$Old)
sacCer3.chrom$New <- paste0("chr", sacCer3.chrom$New)


## TF <-> Experiment ID table 
TF.exp.tab <- data.frame()
TF.exp.tab.file <- file.path(out.folder, "ChExMix_SacCer3_TF_exp_ID.txt")
dataset.counter <- 0

## Iterate over each motif file
for (f in list.files(in.folder)) {
  
  dataset.counter <- dataset.counter + 1
  
  ## Get protein name
  tf <- unlist(strsplit(f,split = "\\."))[1]
  dataset.id <- paste(id, tf, sep = "_")
  
  message("; Processing ", tf, " peaks - Dataset ", dataset.counter)
  
  ## Set new motif file name
  new.bed.file <- file.path(out.folder, dataset.id, paste0(dataset.id, suffix))
  
  ## Rename the chromosomes
  #message("; Renaming the chromosomes using roman numbers")
  old.bed <- fread(file.path(in.folder, f))
  
  if (nrow(old.bed) <= 1) {
    message("; Empty file: ", tf, " - Skipped")
  } else {
    
    new.bed <- merge(old.bed, sacCer3.chrom, by.x = "V1", by.y = "Old") %>% 
      select(New, V2, V3, V4)
    
    ## Create directory for separated file
    new.dir <- dirname(new.bed.file)
    #message("; Creating output directory: ", new.dir)
    dir.create(new.dir, recursive = T, showWarnings = F)
    
    ## TF-ExpID table
    tf.id.info <- data.frame(A = dataset.id,
                             B = dataset.id,
                             C = tf,
                             D = nrow(new.bed))
    TF.exp.tab <- rbind(TF.exp.tab, tf.id.info)
    
    ## Export the new BED file in its own folder
    #message("; Exporting parsed BED file: ", new.bed.file)
    fwrite(new.bed, file = new.bed.file, sep = "\t", row.names = F, col.names = F)
    
  }
}

## Export TF-ExperimentID table
TF.exp.tab <- TF.exp.tab %>% 
                arrange(desc(D))
message("; Exporting TF - Experiment ID table: ", TF.exp.tab.file)
fwrite(TF.exp.tab, file = TF.exp.tab.file, sep = "\t", row.names = F, col.names = F)
