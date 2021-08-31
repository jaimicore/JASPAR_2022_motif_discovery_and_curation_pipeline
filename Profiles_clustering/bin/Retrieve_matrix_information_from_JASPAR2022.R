#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "jsonlite",
                        "optparse",
                        "purrr",
                        "tidyverse")

for (lib in required.libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = T))
}

## How to run: 
## Rscript R-scripts/Retrieve_matrix_information_from_JASPAR2022.R -t Vertebrates -c CORE -o metadata_script_output_testing -v nonredundant

####################
## Read arguments ##
####################
message("; Reading arguments from command-line.")

option_list = list(
  make_option( c("-t", "--taxon"), type = "character", default = "Vertebrates", help = "Taxon available in JASPAR: Vertebrates | Plants | Fungi | Insects | Nematodes | Urochordata", metavar = "character"),
  make_option( c("-o", "--output_directory"), type = "character", default = ".", help = "Output folder", metavar = "character"),
  make_option( c("-c", "--collection"), type = "character", default = "CORE", help = "Collections available in JASPAR: CORE | UNVALIDATED", metavar = "character"),
  make_option( c("-v", "--version"), type = "character", default = "redundant", help = "By default all versions will be returned (option: redundant). For latest versions of profiles: nonredundant", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

taxon        <- opt$taxon
out.folder   <- opt$output_directory
collection   <- opt$collection
version      <- opt$version

# ## Debug:
# taxon      <- "Nematodes"
# out.folder <- "/run/user/316574/gvfs/sftp:host=biotin2.hpc.uio.no/storage/mathelierarea/processed/ieva/projects/JASPAR_2022/metadata_script_output_testing"
# ## out.folder <- "/run/user/280010/gvfs/sftp:host=biotin2.hpc.uio.no,user=jamondra/storage/mathelierarea/processed/jamondra/Projects/JASPAR/JASPAR_2022/Annotation_table/Example"
# collection <- "UNVALIDATED"
# version    <- "rnonredundant"


###############
## Functions ##
###############
# profileID <- "UN0105.1"
get_info_from_api <- function(profileID) {
  
  indiv.jaspar.url <- paste0("http://testjaspar.uio.no/api/v1/matrix/", profileID, "/?format=json")
  indiv.mat.info <- fromJSON(indiv.jaspar.url)
  indiv.mat.info
} 


###############################
## Creating output directory ##
###############################

message("; Creating output directory.")
dir.create(out.folder, showWarnings = F, recursive = T)

#####################################
## Retrieving information from API ##
#####################################

message("; Retrieving ", collection, " collection information from Jaspar API for ", taxon, ".")

## Getting the page size:
if (version == "nonredundant") {
  initial.jaspar.url <- paste0("http://testjaspar.uio.no/api/v1/matrix/?collection=", collection, "&tax_group=", taxon, "&version=latest")
} else {
  initial.jaspar.url <- paste0("http://testjaspar.uio.no/api/v1/matrix/?collection=", collection, "&tax_group=", taxon)
}


## Split the total number of entries in groups of X elements (1000 in this case)
group.size     <- 1000
initial_result <- fromJSON(initial.jaspar.url)
nb_matrices    <- initial_result$count
nb_pages       <- ceiling(nb_matrices/group.size)
nb_pages       <- ifelse(nb_pages == 0, yes = 1, no = nb_pages)  ## When there are less than 1000 entries in a taxon, the ceiling function will be 0 


## Final table:
complete_profiles_tab <- vector(mode = "list", length = nb_pages)

## Iterating through pages:
for (i in 1:nb_pages) {
  message("; Quering page number: ", i)
  
  # ## Requesting all matrices:
  if (version == "nonredundant") {
    jaspar.url <- paste0("http://testjaspar.uio.no/api/v1/matrix/?page=", i, "&page_size=", group.size,"&collection=", collection, "&tax_group=", taxon, "&version=latest&?format=json")
  } else {
    jaspar.url <- paste0("http://testjaspar.uio.no/api/v1/matrix/?page=", i, "&page_size=", group.size,"&collection=", collection, "&tax_group=", taxon, "&?format=json")
  }
  # jaspar.url <- paste0("http://testjaspar.uio.no/api/v1/matrix/?page=", i, "&page_size=1000&collection=CORE&?tax_group=", taxon, "&version=latest&?format=json")
  print(jaspar.url)
  result <- fromJSON(jaspar.url)
  
  # print(result$results$matrix_id)
  ## Getting all information:
  
  # all.profiles.info <- map(result$results$matrix_id, get_info_from_api)
  
  ## If the size of the oputput vector is known, initilize it before calling the iteration function
  ## Allocate memory on the fly is computationally expensive and makes the program slower
  all.profiles.info <- vector(mode = "list", length = length(result$results$matrix_id))
  all.profiles.info <- lapply(result$results$matrix_id, get_info_from_api)
  
  ## Fields:
  # names(all.profiles.info)
  #
  # [1] "pubmed_ids"    "description"   "family"        "pfm"           "tax_group"     "matrix_id"     "sequence_logo" "remap_tf_name"
  # [9] "pazar_tf_ids"  "versions_url"  "collection"    "base_id"       "class"         "tffm"          "tfe_ids"       "name"        
  # [17] "tfbs_shape_id" "uniprot_ids"   "sites_url"     "species"       "alias"         "version"       "unibind"       "type"        
  # [25] "symbol"  
  
  ## Mapping all the info:
  all.profiles.info.subset <- map(all.profiles.info, `[`, c("base_id", "version",
                                                            "name", "uniprot_ids",
                                                            "tax_group",
                                                            "class", "family",
                                                            "species",
                                                            "collection", "type",
                                                            "source", "pubmed_ids", "comment", "tax_group"))
  ### Jaspar curation table format:
  ## PWM
  ## current_BASE_ID
  ## current_VERSION
  ## TF NAME
  ## Uniprot
  ## TAX_ID
  ## class
  ## family
  ## TFBSshape ID
  ## Data
  ## Source
  ## Validation
  ## Comment
  ## Addition or Upgrade or Non-validated (A or U or N )
  
  ## Species is a nested list with two entries: name and tax_id.
  ## They must be processed separately
  species.df <- data.frame( species = do.call(rbind, lapply(map(all.profiles.info.subset, c("species", "name")), paste, collapse = "::") ))
  tax.id.df <- data.frame( tax_id = do.call(rbind, lapply(map(all.profiles.info.subset, c("species", "tax_id")), paste, collapse = "::") ))
  
  ## Family/Class/Uniprot_ids may contain two or more entries (e.g., dimers), therefore, they must be processed separately
  family.df <- data.frame( family = do.call(rbind, lapply(map(all.profiles.info.subset, "family"), paste, collapse = "::") ))
  class.df <- data.frame( class = do.call(rbind, lapply(map(all.profiles.info.subset, "class"), paste, collapse = "::") ))
  tax_group.df <- data.frame( tax_group = do.call(rbind, lapply(map(all.profiles.info.subset, "tax_group"), paste, collapse = "::") ))
  uniprot.df <- data.frame( uniprot_ids = do.call(rbind, lapply(map(all.profiles.info.subset, "uniprot_ids"), paste, collapse = "::") ))
  
  ## some profiles may contain two pubmed ids (this is a mistake when we curated the database). To avoid problems, we concatenate them
  ## But this problem must be fixed for future releases
  pubmed.df <- data.frame( pubmed_ids = do.call(rbind, lapply(map(all.profiles.info.subset, "pubmed_ids"), paste, collapse = "::") ))
  
  ## Conver list to data.frame
  all.profiles.info.subset <- map(all.profiles.info.subset, `[`, c("base_id", "version", "name", "collection", "type", "source", "comment"))
  oldw <- getOption("warn")
  options(warn = -1)
  all.profiles.info.tab <- rbindlist(all.profiles.info.subset, fill = TRUE)
  options(warn = oldw)
  # colnames(all.profiles.info.tab) <- c("base_id", "version", "name", "collection", "type", "source", "comment")
  
  ## Concat all the data.frames
  all.profiles.info.tab.clean <-
    cbind(all.profiles.info.tab, species.df, tax.id.df, class.df, family.df, uniprot.df, pubmed.df, tax_group.df) %>%
    dplyr::select(base_id, version,
                  name, uniprot_ids,
                  tax_id,
                  class, family,
                  collection, type,
                  source, pubmed_ids, comment, tax_group)
  
  
  ## Adding to the final list:
  complete_profiles_tab[[i]] <- all.profiles.info.tab.clean
}

complete_profiles_tab <- do.call("rbind", complete_profiles_tab) %>% 
  dplyr::filter(tax_group == tolower(taxon))


message("; Exporting table")
fwrite(complete_profiles_tab,
       sep = "\t",
       file = file.path(out.folder, paste0("JASPAR_2022_metadata_", version, "_", collection, "_", taxon,"_table.tab")))

###################
## End of script ##
###################