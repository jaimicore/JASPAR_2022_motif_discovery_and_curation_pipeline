#####################
## Load R packages ##
#####################
required.libraries <- c("optparse",
                        "universalmotif",
                        "data.table")

for (lib in required.libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = T))
}


####################
## Read arguments ##
####################
option_list = list(
  
  make_option(c("-m", "--input_motif"), type = "character", default = NULL,
              help = "Input motif file. (Mandatory)", metavar = "character"),
  
  make_option(c("-n", "--format"), type = "character", default = NULL,
              help = "Input motif format. (Mandatory)", metavar = "character"),
  
  make_option(c("-b", "--both"), type = "numeric", default = NULL,
              help = "Trim b nucleotides in both sides", metavar = "numeric"),
  
  make_option(c("-l", "--left"), type = "numeric", default = NULL,
              help = "Trim l nucleotides in left side", metavar = "numeric"),
  
  make_option(c("-r", "--right"), type = "numeric", default = NULL,
              help = "Trim r nucleotides in right side", metavar = "numeric"),
  
  make_option(c("-f", "--from"), type = "numeric", default = NULL,
              help = "Keep nucleotides from f position", metavar = "numeric"),
  
  make_option(c("-t", "--to"), type = "numeric", default = NULL,
              help = "Keep nucleotides until t position", metavar = "numeric"),
  
  make_option(c("-o", "--output_directory"), type = "character", default = NULL,
              help = "Output directory to export the trimmed motifs.", metavar = "character")
  
);

message("; Reading arguments from the command line.")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


########################
## Set variable names ##
########################
motif.file    <- opt$input_motif
motif.format  <- tolower(opt$format)
trim.b        <- opt$both
trim.l        <- opt$left
trim.r        <- opt$right
from          <- opt$from
to            <- opt$to
out.dir       <- opt$output_directory

## Debug:
# setwd("JASPAR_2022")
# motif.file <- "ZFP42_rounded.jaspar"
# motif.file <- "ZFP42.jaspar"
# motif.file <- "MA1651.1.jaspar"

# motif.format <- "jaspar"
# from <- 1
# to <- 2
# trim.b <- NULL
# trim.l <- NULL
# trim.r <- NULL


############################
## Check input parameters ##
############################
b.flag <- 0
l.flag <- 0
r.flag <- 0
t.flag <- 0
f.flag <- 0



if (!is.null(trim.b)) {
  if (trim.b >= 0) {
    b.flag <- 1
  } else {
    message("; --both must be a positive integer")
  }
}


if (!is.null(trim.l)) {
  
  if (trim.l >= 0) {
    l.flag <- 1
  } else {
    message("; --left must be a positive integer")
  }
  
}


if (!is.null(trim.r)) {
  
  if (trim.r >= 0) {
    r.flag <- 1
  } else {
    message("; --right must be a positive integer")
  }
}


if (!is.null(from)) {
  
  if (from >= 0) {
    f.flag <- 1
  } else {
    message("; --from must be a positive integer")
  }
}


if (!is.null(to)) {
  
  if (to >= 0) {
    t.flag <- 1
  } else {
    message("; --to must be a positive integer")
  }
}
  

# print(b.flag)
# print(l.flag)
# print(r.flag)
# print(f.flag)
# print(t.flag)



if (b.flag & sum(c(l.flag, r.flag) > 0) ) {
  stop("; --both option is not compatible with --left or --right")
}

if (b.flag & sum(c(t.flag, f.flag) > 0) ) {
  stop("; --both option is not compatible with --from or --to")
}

if (sum(c(t.flag, f.flag) > 0 & sum(c(l.flag, r.flag) > 0)) ) {
  stop("; --left and --right options are not compatible with --from or --to")
}


supported.formats <- c("jaspar")
if (!tolower(motif.format) %in% supported.formats) {
  message("; Input format ", motif.format, " not recognized.")
  stop("Accepted motif formats: ", paste(supported.formats, collapse = ", "))
}


#####################################
## Read motif according its format ##
#####################################
message("; Reading input motif file: ", motif.file)
if (motif.format == "jaspar") {
  
  
  dt.motif <- fread(motif.file)
  
  ## JASPAR header name
  relevant.header.pos <- unlist(sapply(colnames(dt.motif), grep, pattern = "^V\\d+$", invert = T, value = T))
  
  if (grepl(relevant.header.pos, pattern = "\t")) {
    relevant.header.pos <- strsplit(relevant.header.pos, split = "\t")
  }
  
  header.name <- paste(relevant.header.pos, collapse = " ")
  
  annotations <- sapply(relevant.header.pos, gsub, pattern = "^>", replacement = "")
  motif.id <- annotations[1]
  
  if(length(annotations) > 1){
    motif.name <- annotations[2]
  } else {
    motif.name <- ""
  }
  
  
  dt.counts <- dt.motif[, 2:(ncol(dt.motif)-1)]
  
  if(all(dt.counts[,1] == "[")){
    dt.counts <- dt.counts[,-1]
  }
  
  dt.counts <- data.frame(apply(dt.counts, 2, gsub, pattern = "\\[", replacement = ""))
  
  dt.counts.rounded <- round(apply(dt.counts,2,as.numeric))
  
  
  motif.pfm <- universalmotif::create_motif(dt.counts.rounded, alphabet = "ACTG", type = "PFM")
  
  
  # motif.pfm  <- read_jaspar(motif.file)s
  # motif.id   <-  motif.pfm["name"]
  # motif.name <-  motif.pfm["altname"]

}



##############################
## Define missing variables ##
##############################

motif.size <- ncol(motif.pfm["motif"])

if (is.null(out.dir)) {
  out.dir <- dirname(motif.file)
} else {
  dir.create(out.dir, showWarnings = F, recursive = T)
  message("; Output file: ", out.dir)
}


if (b.flag) {
  
  from <- trim.b
  to   <- motif.size - trim.b
  
} else if (l.flag | r.flag) {
  
  from <- ifelse(is.null(trim.l), yes = 1, no = trim.l)
  to   <- ifelse(is.null(trim.r), yes = motif.size, no = motif.size - trim.r)
  
} else if (f.flag | t.flag) {
  
  from <- ifelse(is.null(from), yes = 1, no = from)
  to   <- ifelse(is.null(to), yes = motif.size, no = to)
  
}


if (to > motif.size) {
  to <- motif.size
  message("; Right side adjusted to motif size: ", motif.size)
}



###################
## Core function ##
###################
trim_motif_friseur <- function(m = NULL,
                               from = 0,
                               to   = 0) {
  
    new.motif <- m["motif"][,from:to]
    return(new.motif)
}


## Trim the input motif + Rename it
trimmed.motif            <- trim_motif_friseur(m = motif.pfm, from, to)
trimmed.motif            <- convert_type(trimmed.motif, "PCM")
trimmed.motif["name"]    <- motif.id 
trimmed.motif["altname"] <- motif.name


####################################################
## Write trimmed motif according its input format ##
####################################################
new.motif.file <- file.path(out.dir, paste0(basename(motif.file), ".trimmed"))
message("; Exporting trimmed motif: ", new.motif.file)

if (motif.format == "jaspar") {
  
  if (file.exists(new.motif.file)) {
    file.remove(new.motif.file)
  }
  
  write_jaspar(motifs = trimmed.motif, file = new.motif.file)
}


#######################
## End of the Script ##
#######################
message("; Motif trimming done.")
