#############################
## Load required libraries ##
#############################
required.packages = c("dplyr",
                      "data.table",
                      "optparse",
                      "purrr",
                      "tidyr")


for (lib in required.packages) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}


####################
## Read arguments ##
####################
option_list = list(
  
  make_option(c("-m", "--metadata_table"), type = "character", default = "NULL", 
              help = "JASPAR CORE collection metadata (one table for taxon). (Mandatory) ", metavar = "character"),
  
  make_option(c("-c", "--cluster_assignment"), type = "character", default = NULL, 
              help = "Cluster table from RSAT matrix-clistering (Mandatory)", metavar = "character"),
  
  make_option(c("-s", "--suffix"), type = "character", default = "example_collection", 
              help = "Suffix to be added to the exported table file names", metavar = "character"),
  
  make_option(c("-o", "--output_directory"), type = "character", default = 0,
              help = "Output directory where the resulting tables will be exported. (Mandatory) ", metavar = "character")
);
message("; Reading arguments from command-line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


## Set variable names
metadata.tab.file   <- opt$metadata_table
cluster.assign.file <- opt$cluster_assignment
out.dir             <- opt$output_directory
file.suffix         <- opt$suffix


###########
## Debug ##
###########
# metadata.tab.file   <- "/home/jamondra/Downloads/JASPAR_2022_metadata_nonredundant_CORE_vertebrates_table_parsed.tab"
# cluster.assign.file <- "/home/jamondra/Downloads/clusters.tab"
# out.dir             <- "/home/jamondra/Downloads/Annotate_FBPs"
# file.suffix         <- "JASPAR_2022_nonredundant_CORE_vertebrates"


###########################
## Create output folders ##
###########################
dir.create(out.dir, recursive = T)


#########################
## Read metadata table ##
#########################
message("; Reading metadata table: ", metadata.tab.file)
metadata <- fread(metadata.tab.file) %>% 
              mutate(ID = paste0(base_id, ".", version))


####################
## Read FBP table ##
####################
message("; Reading cluster assignment table: ", cluster.assign.file)
cluster.assignment <- fread(cluster.assign.file, header = F) %>% 
                        rename(cluster  = V1,
                               ID       = V2) %>% 
                        separate_rows(ID, sep = ",") %>% 
                        mutate(ID      = gsub(ID, pattern = ".+_m\\d+_", replacement = ""),
                               cluster = as.numeric(gsub(cluster, pattern = "cluster_", replacement = ""))) %>% 
                        mutate(ID = gsub(ID, pattern = "_", replacement = "\\."))
cluster.assignment <- data.table(cluster.assignment)
  
  
###############################################
## Combine metadata with cluster information ##
###############################################
  
## Merge tables to obtain the TF class/family
message("; Combining metadata with cluster assignation")
motif.names <- merge.data.table(cluster.assignment, metadata, by = "ID") %>% 
                  select(ID, name, cluster, class, family)
  
  
###################################################
## Separate singleton and non-singleton clusters ##
###################################################
message("; Separating singleton and non-singleton clusters")
nb.motifs.per.cluster <- table(motif.names$cluster)
singleton.clusters    <- as.numeric(names(nb.motifs.per.cluster[which(nb.motifs.per.cluster == 1)])) ## Clusters of size 1
nb.motifs.per.cluster <- nb.motifs.per.cluster[which(nb.motifs.per.cluster > 1)] ## Remove clusters of size 1
cluster.order         <- as.numeric(names(sort(nb.motifs.per.cluster, decreasing = T)))
  
message("; Singleton clusters: ", length(singleton.clusters))
message("; Non-Singleton clusters: ", length(cluster.order))
  
  
#################################################
## Count the number of classes on each cluster ##
#################################################
message("; Counting the number of TF classes on each cluster")
classes.per.cluster.list <- vector(mode = "list", length = length(cluster.order))
classes.per.cluster.list <- lapply(cluster.order, function(cl){
    
  cluster.tf.classes <- subset(motif.names, cluster == cl) %>% 
                          separate_rows(name, class, sep = "::") %>% 
                          distinct()
    
  cluster.tf.classes <- cluster.tf.classes %>% 
                          group_by(class) %>% 
                            mutate(Names = paste(name, collapse = ","),
                                   ID    = paste(ID, collapse = ","),
                                   N = n()) %>% 
                            mutate(N_rel   = N/nrow(cluster.tf.classes)) %>%
                            ungroup() %>% 
                            mutate(Max_N   = max(N_rel)) %>% 
                            # dplyr::filter(N_rel == Max_N) %>% 
                            select(cluster, N, N_rel, class, Names, ID) %>% 
                            unique() %>% 
                            arrange(desc(cluster))
    
    
    cluster.tf.classes
    
})
classes.per.cluster.list <- rbindlist(classes.per.cluster.list)  
  
  
#####################################################
## Select the most common TF class on each cluster ##
#####################################################
message("; Selecting the most common class on each cluster")
cluster.fbp <- classes.per.cluster.list %>% 
                group_by(cluster) %>% 
                filter(N_rel == max(N_rel)) %>% 
                arrange(cluster) %>% 
                select(cluster, N_rel, class) %>% 
                mutate(class = paste(class, collapse = "&")) %>% 
                unique() %>% 
                select(cluster, class) %>% 
                rename(fbp = class) %>% 
                data.table()

  
singleton.fbp <- motif.names %>% 
                  dplyr::filter(cluster %in% singleton.clusters) %>% 
                  select(cluster, name) %>% 
                  rename(fbp = name)
  
cluster.fbp <- rbind(cluster.fbp, singleton.fbp)
  
# length(unique(cluster.fbp$cluster))
# length(cluster.fbp$cluster)
  

#################################
## Assign a number to each FBP ##
#################################
message("; Assigning names to FBPs")
## Use rowid function to assign a number to each entry in a vector
## The first time a name appears will have assigned the number 1, the second time number 2 and so on
cluster.fbp$Row_ID <- rowid(cluster.fbp$fbp)
  
cluster.fbp <- cluster.fbp %>% 
                mutate(FBP = paste(fbp, Row_ID, sep = "::")) %>% 
                select(cluster, FBP)
  
cluster.fbp <- merge(cluster.fbp, motif.names, by = "cluster") %>% 
                group_by(cluster) %>% 
                mutate(TFs = paste(name, collapse = ","),
                       ID  = paste(ID, collapse = ","),
                       N   = n()) %>% 
                select(cluster, FBP, N, TFs, ID) %>% 
                distinct()
  
cluster.fbp <- cluster.fbp %>% 
                mutate(FBP = ifelse(N == 1,
                                    yes = gsub(FBP, pattern = "::1", replacement = ""),
                                    no  = FBP))
  


#############################
## Export annotation table ##
#############################
cluster.fbp.file <- file.path(out.dir, paste0("FBP_names_", file.suffix, ".txt"))
fwrite(cluster.fbp, file = cluster.fbp.file, sep = "\t")
  
classes.per.cluster.list.file <- file.path(out.dir, paste0("Classes_per_cluster_", file.suffix, ".txt"))
fwrite(classes.per.cluster.list, file = classes.per.cluster.list.file, sep = "\t")

