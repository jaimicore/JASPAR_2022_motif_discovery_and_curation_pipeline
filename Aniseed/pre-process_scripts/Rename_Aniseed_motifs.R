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


# /run/user/280010/gvfs/sftp:host=biotin3.hpc.uio.no,user=jamondra
input.motif.folder    <- "/storage/scratch/JASPAR_2022/Aniseed/data/Best_round/pfm_best_round"
output.motif.folder   <- "/storage/scratch/JASPAR_2022/Aniseed/data/Best_round"
annotation.table      <- "/storage/scratch/JASPAR_2022/Aniseed/data/Best_round/SELEX_project_TableS1_S3_Nitta_et_al_2019.csv"
organism              <- "Ciona_robusta"

## Create output directory
pfm.output.dir <- file.path(output.motif.folder, "PFMs")
logo.output.dir <- file.path(output.motif.folder, "Logos")
dir.create(pfm.output.dir, showWarnings = F, recursive = T)
dir.create(logo.output.dir, showWarnings = F, recursive = T)


#############################################
## Read and parse Aniseed annotation table ##
#############################################
message("; Reading Aniseed annotation table")
aniseed.tab <- fread(annotation.table, header = T)

aniseed.tab <- aniseed.tab[,c(2:7)]

colnames(aniseed.tab) <- c("TF_fam",
                           "TF_class",
                           "Gene_name",
                           "Alt_gene_name",
                           "DBD",
                           "SELEX_clone")


#############################
## Parse the Aniseed table
##
## Keep entries with SELEX_clone ID
## Transform the interproID colums to one entry per line
aniseed.tab <- aniseed.tab %>% 
  dplyr::filter(SELEX_clone != "")
aniseed.tab$SELEX_clone <- substr(aniseed.tab$SELEX_clone, start = 1, stop = 5)

## New file with all motifs
all.motifs.file <- file.path(pfm.output.dir, paste0("All_", organism, "_motifs_concatenated.tf"))

TF.name.counter <- list()

## Iterate over each motif file
for (f in list.files(input.motif.folder)) {
  
  ## Get Selex probe ID
  selex.id <- gsub(f, pattern = "^([A-Za-z0-9]+)_\\d+.+$", replacement = "\\1")
  selex.id <- substr(selex.id, start = 1, stop = 5)

  tf <- as.vector(subset(aniseed.tab, SELEX_clone == selex.id)$Gene_name)
  tf <- unique(tf)
  
  ## Renamed problematic names
  ##
  ## Tbx2/3  -> Tbx2_3
  tf <- gsub(tf, pattern = "/", replacement = "-")
  
  ## Process the motif or skip
  if (length(tf) == 0) {
    message("; SELEX probe: ", selex.id, " not found in annotation table. Skip.")
  } else {
    
    message("; Processing ", tf, " motif - SELEX probe: ", selex.id)

    ## Count the number of times this TF appears in the file list
    if (is.null(TF.name.counter[[tf]])) {
      TF.name.counter[[tf]] <- NULL
      TF.name.counter[[tf]] <- 1
    } else {
      TF.name.counter[[tf]] <- TF.name.counter[[tf]] + 1
    }
    

    ## Set new motif file name
    new.tmp.file <- file.path(pfm.output.dir, paste0(tf, "_pfm_", TF.name.counter[[tf]],".jaspar"))
    tmp.files.list <- append(tmp.files.list, new.tmp.file)
    
    ## Add a header to the motifs
    add.header.cmd <- paste0( '(echo ">',tf,'"_v', TF.name.counter[[tf]], ' && cat ', file.path(input.motif.folder, f), ') > ', new.tmp.file )
    # message("; ", add.header.cmd)
    system(add.header.cmd)
    
    
    ## RSAT convert-matrix command
    ## From cis-bp to transfac format
    ## Add manually the TF name
    new.tf.file <- file.path(pfm.output.dir, paste0(tf, "_pfm_", TF.name.counter[[tf]], ".tf"))
    convert.matrix.cmd <- paste0("convert-matrix -i ", new.tmp.file, " -from jaspar -to tf -multiply 1000  > ", new.tf.file)
    # message("; ", convert.matrix.cmd)
    system(convert.matrix.cmd)
    
    ## Concat all motifs in a single file
    cat.all.motifs.cmd <- paste0("cat ", new.tf.file, " >> ", all.motifs.file)
    system(cat.all.motifs.cmd)
    
    ## Remove the jaspar temporal file
    rm.tmp.file <- paste0("rm -rf ", new.tmp.file)
    system(rm.tmp.file)
  }
}


## Generate logos
logos.cmd <- paste("convert-matrix -i ", all.motifs.file, " -from tf -to tf -return logo -logo_dir ", logo.output.dir)
message("; ", logos.cmd)
system(logos.cmd)




