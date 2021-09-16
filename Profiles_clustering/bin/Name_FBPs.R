#############################
## Load required libraries ##
#############################

## List of packages to install from CRAN
required.packages = c("dplyr",
                      "data.table",
                      "purrr",
                      "tidyr")


for (lib in required.packages) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}



data.dir        <- "/home/jamondra/Downloads/JASPAR_2022_plots/Archetype_names/data"
results.dir     <- "/home/jamondra/Downloads/JASPAR_2022_plots/Archetype_names/results"


taxa <- c("vertebrates", "plants", "insects", "fungi", "urochordates", "nematodes")
for (taxon in taxa) {
  
  metadata.tab    <- file.path(data.dir, taxon, paste0("JASPAR_2022_metadata_nonredundant_CORE_", taxon, "_table_parsed.tab")) 
  archetype.files <- file.path(data.dir, taxon, "archetypes", list.files(file.path(data.dir, taxon, "archetypes"), pattern = "motifs"))
  
  
  ###########################
  ## Create output folders ##
  ###########################
  tab.results.dir <- file.path(results.dir, taxon)
  dir.create(tab.results.dir, recursive = T)
  
  
  #########################
  ## Read metadata table ##
  #########################
  metadata <- fread(metadata.tab) %>% 
    mutate(ID = paste0(base_id, ".", version))
  
  
  ##################################################
  ## Read archetype files and merge with metadata ##
  ##################################################
  motif.names <- rbindlist(map(archetype.files, fread)) %>% 
    select(id2, name2, cluster) %>% 
    rename(ID = id2)
  
  ## Merge tables to obtain the TF class/family
  motif.names <- merge.data.table(motif.names, metadata, by = "ID") %>% 
                  select(ID, name, cluster, class, family)
  
  
  ###################################################
  ## Separate singleton and non-singleton clusters ##
  ###################################################
  nb.motifs.per.cluster <- table(motif.names$cluster)
  singleton.clusters    <- as.numeric(names(nb.motifs.per.cluster[which(nb.motifs.per.cluster == 1)])) ## Clusters of size 1
  nb.motifs.per.cluster <- nb.motifs.per.cluster[which(nb.motifs.per.cluster > 1)] ## Remove clusters of size 1
  cluster.order         <- as.numeric(names(sort(nb.motifs.per.cluster, decreasing = T)))
  
  message("; Singleton clusters: ", length(singleton.clusters))
  message("; Non-Singleton clusters: ", length(cluster.order))
  
  
  #################################################
  ## Count the number of classes on each cluster ##
  #################################################
  classes.per.cluster.list <- vector(mode = "list", length = length(cluster.order))
  classes.per.cluster.list <- lapply(cluster.order, function(cl){
    
    cluster.tf.classes <- subset(motif.names, cluster == cl) %>% 
      separate_rows(name, class, sep = "::") %>% 
      unique()
    
    cluster.tf.classes <- cluster.tf.classes %>% 
      group_by(class) %>% 
      mutate(Names = paste(name, collapse = ","),
             N = n()) %>% 
      mutate(N_rel   = N/nrow(cluster.tf.classes)) %>%
      ungroup() %>% 
      mutate(Max_N   = max(N_rel)) %>% 
      # dplyr::filter(N_rel == Max_N) %>% 
      select(cluster, N, N_rel, class, Names) %>% 
      unique() %>% 
      arrange(desc(cluster))
    
    
    cluster.tf.classes
    
  })
  classes.per.cluster.list <- rbindlist(classes.per.cluster.list)  
  
  
  #####################################################
  ## Select the most common TF class on each cluster ##
  #####################################################
  cluster.fbp <- classes.per.cluster.list %>% 
    group_by(cluster) %>% 
    filter(N_rel == max(N_rel)) %>% 
    arrange(cluster) %>% 
    select(cluster, N_rel, class) %>% 
    mutate(class = paste(class, collapse = "&")) %>% 
    distinct() %>% 
    select(cluster, class) %>% 
    rename(fbp = class)
  
  singleton.fbp <- motif.names %>% 
    dplyr::filter(cluster %in% singleton.clusters) %>% 
    select(cluster, name) %>% 
    rename(fbp = name)
  
  cluster.fbp <- rbind(cluster.fbp, singleton.fbp)
  
  
  length(unique(cluster.fbp$cluster))
  length(cluster.fbp$cluster)
  
  #################################
  ## Assign a number to each FBP ##
  #################################
  cluster.fbp$Row_ID <- rowid(cluster.fbp$fbp)
  
  cluster.fbp <- cluster.fbp %>% 
    mutate(FBP = paste(fbp, Row_ID, sep = "::")) %>% 
    select(cluster, FBP)
  
  cluster.fbp <- merge(cluster.fbp, motif.names, by = "cluster") %>% 
    group_by(cluster) %>% 
    mutate(TFs = paste(name, collapse = ","),
           N   = n()) %>% 
    select(cluster, FBP, N, TFs) %>% 
    distinct()
  
  cluster.fbp <- cluster.fbp %>% 
    mutate(FBP = ifelse(N == 1,
                        yes = gsub(FBP, pattern = "::1", replacement = ""),
                        no  = FBP))
  
  cluster.fbp.file <- file.path(tab.results.dir, paste0("FBP_names_", taxon, ".txt"))
  fwrite(cluster.fbp, file = cluster.fbp.file, sep = "\t")
  
  classes.per.cluster.list.file <- file.path(tab.results.dir, paste0("Classes_per_cluster_", taxon, ".txt"))
  fwrite(classes.per.cluster.list, file = classes.per.cluster.list.file, sep = "\t")
}
