#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "jsonlite",
                        "optparse",
                        "purrr")

for (lib in required.libraries) {
  suppressPackageStartupMessages(library(lib, character.only=TRUE, quietly = T))
}

## How to run: 
## 
## Rscript R-scripts/Retrieve_matrix_information_from_JASPAR2022.R -t Vertebrates -o .


####################
## Read arguments ##
####################
message("; Reading arguments from command-line")
option_list = list(
  make_option( c("-t", "--taxon"), type = "character", default = "Vertebrates", help = "Taxon available in JASPAR: Vertebrates | Plants | Fungi | Insects | Nematodes | Urochordata", metavar = "character"),
  make_option( c("-o", "--output_directory"), type = "character", default = ".", help = "Output folder", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

taxon        <- opt$taxon
out.folder   <- opt$output_directory

# ## Debug:
# taxon <- "Nematodes"
# out.folder <- "/run/user/316574/gvfs/sftp:host=biotin2.hpc.uio.no/storage/mathelierarea/processed/ieva/projects/JASPAR_2022"
###############################
## Creating output directory ##
###############################

dir.create(out.folder, showWarnings = F, recursive = T)

#####################################
## Retrieving information from API ##
#####################################

message("; Retrieving information from Jaspar API.")

## Getting the page size:
initial.jaspar.url <- file.path("http://testjaspar.uio.no/api/v1/matrix/?collection=CORE&?tax_group=", taxon, "&version=latest")
# initial.jaspar.url <- paste0("http://testjaspar.uio.no/api/v1/matrix/?collection=CORE&?tax_group=", taxon, "&version=latest")
initial_result <- fromJSON(initial.jaspar.url)
nb_matrices <- initial_result$count
nb_pages <- ceiling(nb_matrices/1000)

## Final table:
complete_profiles_tab <- vector(mode = "list", length = nb_pages)

## Iterating through pages:
for (i in 1:nb_pages) {
  message("; Quering page number: ", i)
  
  ## Requesting all matrices:
  jaspar.url <- paste0("http://testjaspar.uio.no/api/v1/matrix/?page=", i, "&page_size=1000&collection=CORE&?tax_group=", taxon, "&version=latest&?format=json")
  result <- fromJSON(jaspar.url)
  
  #print(result$results$matrix_id)
  ## Getting all information:
  all.profiles.info <- sapply(result$results$matrix_id, function(id){
    
    #print(id)
    indiv.jaspar.url <- paste0("http://testjaspar.uio.no/api/v1/matrix/", id, "/?format=json")
    indiv.mat.info <- fromJSON(indiv.jaspar.url)
    
    indiv.mat.info
    
  })
  
  ## Fields:
  # names(all.profiles.info)
  #
  # [1] "pubmed_ids"    "description"   "family"        "pfm"           "tax_group"     "matrix_id"     "sequence_logo" "remap_tf_name"
  # [9] "pazar_tf_ids"  "versions_url"  "collection"    "base_id"       "class"         "tffm"          "tfe_ids"       "name"        
  # [17] "tfbs_shape_id" "uniprot_ids"   "sites_url"     "species"       "alias"         "version"       "unibind"       "type"        
  # [25] "symbol"  
  
  ## Selecting relevant info:
  message("; Parsing table.")
  
  ## Mapping all the info:
  all.profiles.info.subset <- map(all.profiles.info, `[`, c("base_id", "version",
                                                            "name", "uniprot_ids",
                                                            "tax_group",
                                                            "class", "family",
                                                            "species",
                                                            "tfbs_shape_id", "type",
                                                            "source", "pubmed_ids", "comment"))
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
  all.profiles.info.subset <- map(all.profiles.info.subset, `[`, c("base_id", "version", "name", "tax_group", "tfbs_shape_id", "type", "source", "comment"))
  all.profiles.info.tab <- rbindlist(all.profiles.info.subset, fill = TRUE)
  
  ## Concat all the data.frames
  all.profiles.info.tab.clean <-
    cbind(all.profiles.info.tab, species.df, tax.id.df, class.df, family.df, uniprot.df, pubmed.df) %>%
    dplyr::select(base_id, version,
                  name, uniprot_ids,
                  tax_id,
                  class, family,
                  tfbs_shape_id, type,
                  source, pubmed_ids, comment, tax_group)
  
  
  ## Adding to the final list:
  complete_profiles_tab[[i]] <- all.profiles.info.tab.clean
}

complete_profiles_tab <- do.call("rbind", complete_profiles_tab)

message("; Exporting table")
fwrite(complete_profiles_tab, sep = "\t", file = file.path(out.folder, paste0("JASPAR_2022_info_", taxon,"_table.tab")))

###################
## End of script ##
###################