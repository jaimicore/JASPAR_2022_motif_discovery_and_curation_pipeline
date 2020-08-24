#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "jsonlite",
                        "optparse",
                        "purrr")
for (lib in required.libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}


## How to run: 
## 
## Rscript Retrieve_matrix_information_from_JASPAR2020.R -t Vertebrates -o .



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

dir.create(out.folder, showWarnings = F, recursive = T)


# setwd(".")
# taxon <- "Vertebrates"
message("; Retrieving information from API")
jaspar.url <- paste0("http://jaspar.genereg.net/api/v1/matrix/?page_size=800&collection=CORE&tax_group=", taxon, "&version=latest&format=json")
result <- fromJSON(jaspar.url)

all.profiles.info <- sapply(result$results$matrix_id, function(id){
  
  indiv.jaspar.url <- paste0("http://jaspar.genereg.net//api/v1/matrix/", id,"/format=json")
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


message("; Parsing table")
## Vertebrates: 746 profiles
all.profiles.info.subset <- map(all.profiles.info, `[`, c("name", "matrix_id", "class", "family", "tax_group", "species", "type", "pubmed_ids", "uniprot_ids"))

## Species is a nested list with two entries: name and tax_id.
## They must be processed separately
species.df <- data.frame( species = do.call(rbind, lapply(map(all.profiles.info.subset, c("species", "name")), paste, collapse = "::") ))
tax.id.df <- data.frame( tax_id = do.call(rbind, lapply(map(all.profiles.info.subset, c("species", "tax_id")), paste, collapse = "::") ))

## Family/Class/Uniprot_ids may contain two or more entries (e.g., dimers), therefore, they must be processed separately
family.df <- data.frame( family = do.call(rbind, lapply(map(all.profiles.info.subset, "family"), paste, collapse = "::") ))
class.df <- data.frame( class = do.call(rbind, lapply(map(all.profiles.info.subset, "class"), paste, collapse = "::") ))
uniprot.df <- data.frame( uniprot_ids = do.call(rbind, lapply(map(all.profiles.info.subset, "uniprot_ids"), paste, collapse = "::") ))

## some profiles may contain two pubmed ids (this is a mistake when we curated the database). To avoid problems, we concatenate them
## But this problem must be fixed for future releases
pubmed.df <- data.frame( pubmed_ids = do.call(rbind, lapply(map(all.profiles.info.subset, "pubmed_ids"), paste, collapse = "::") ))

## Conver list to data.frame
all.profiles.info.subset <- map(all.profiles.info.subset, `[`, c("name", "matrix_id", "tax_group", "type"))
all.profiles.info.tab <- rbindlist(all.profiles.info.subset)


## Concat all the data.frames
all.profiles.info.tab.clean <-
  cbind(all.profiles.info.tab, species.df, tax.id.df, class.df, family.df, uniprot.df, pubmed.df) %>%
  dplyr::select(name, matrix_id, class, family, tax_group, species, tax_id, type, uniprot_ids, pubmed_ids)

message("; Exporting table")
fwrite(all.profiles.info.tab.clean, sep = "\t", file = file.path(out.folder, paste0("Jaspar_2020_info_", taxon,"_table.tab")))
