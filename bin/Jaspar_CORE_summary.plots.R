#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "ggplot2",
                        "optparse",
                        "plotly",
                        "RColorBrewer",
                        "reshape2",
                        "rcartocolor",
                        "tidyr")

for (lib in required.libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}


####################
## Read arguments ##
####################
option_list = list(
  
  make_option(c("-i", "--input_table"), type = "character", default = "NULL", 
              help = "Table with the number of motifs on each JASPAR release. (Mandatory) ", metavar = "character"),
  
  make_option(c("-o", "--output_directory"), type = "character", default = NULL, 
              help = "Output directory to export the results (Mandatory)", metavar = "character"),

  make_option(c("-p", "--plotly_server"), type = "numeric", default = 0,
              help = "Export the figures in a plotly server. [Default \"%default\"] ", metavar = "number"),
);
message("; Reading arguments from command-line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


## Set variable names
results.dir               <- opt$output_directory
motifs.per.taxon.tab.file <- opt$input_table
plotly.export             <- as.numeric(opt$plotly_server)
# motifs.per.taxon.tab.file <- "/home/jamondra/Downloads/Motifs_per_taxon_per_release.txt"
# results.dir               <- "/home/jamondra/Downloads/JASPAR_2022_plots"

dir.create(results.dir, showWarnings = F, recursive = T)

##########################################
## Read JASPAR motifs per release table ##
##########################################
motifs.per.taxon.tab <- fread(motifs.per.taxon.tab.file)

## Add column with the total number of motifs
motifs.per.taxon.tab$All_taxa <- rowSums(motifs.per.taxon.tab[,c("Vertebrates", "Plants", "Insects", "Nematodes", "Fungi", "Urochordates", "Cnidaria", "Protozoa", "Trematodes")],na.rm = T)

## Get the year of the latest release
release.year <- max(motifs.per.taxon.tab$Year, na.rm = T)
releases     <- as.vector(na.omit(motifs.per.taxon.tab$Year))

## Convert the matrix in a data.frame
nb.motfs.per.release <- data.frame(melt(motifs.per.taxon.tab[,c(-2)], id.vars = c("Collection", "Year")))
nb.motfs.per.release <- nb.motfs.per.release %>% 
                          rename(Taxon     = variable,
                                 Nb_motifs = value)
                          
nb.taxa <- length(as.vector(unique(nb.motfs.per.release$Taxon)))

## Add the year column
nb.motfs.per.release$Year <- rep(motifs.per.taxon.tab$Year, nb.taxa)
nb.motfs.per.release$Year <- ordered(nb.motfs.per.release$Year, levels = as.vector(unique(nb.motfs.per.release$Year)))


##########################
## Create color palette ##
##########################

## Determine order of taxon (descending, depending on the number of motifs)
taxon.order.CORE        <- subset(nb.motfs.per.release, Year == release.year & Collection == "CORE" & Taxon != "All_taxa") %>%
                            drop_na(Year, Nb_motifs) %>% 
                            arrange(desc(Nb_motifs))
taxon.order.CORE        <- as.vector(taxon.order.CORE$Taxon)


taxon.order.UNVALIDATED <- subset(nb.motfs.per.release, Year == release.year & Collection == "UNVALIDATED" & Taxon != "All_taxa") %>%
  drop_na(Year, Nb_motifs) %>% 
  arrange(desc(Nb_motifs))
taxon.order.UNVALIDATED <- as.vector(taxon.order.UNVALIDATED$Taxon)

taxon.order <- unique(c(taxon.order.CORE, taxon.order.UNVALIDATED))

tax.cols     <- carto_pal(length(taxon.order), "Bold")
cols         <- c( "#666666", tax.cols)
names(cols)  <- c("All_taxa", taxon.order)

## Order factors by number of motifs
nb.motfs.per.release$Taxon <- factor(nb.motfs.per.release$Taxon, levels = as.vector(c("All_taxa",
                                                                                      taxon.order)))




######################################################
## Generate plots: CORE and UNVALIDATED collections ##
######################################################

for (collection in c("CORE", "UNVALIDATED")) {
  
  message("; Plots for JASPAR ", collection, " collection")
  
  ################
  ## Line chart ##
  ################
  nb.motfs.per.release.gg <- nb.motfs.per.release %>% 
                                dplyr::filter(Collection == collection) %>% 
                                drop_na(Year, Nb_motifs) %>% 
                             ggplot(aes(x = as.factor(Year), y = Nb_motifs, color = Taxon, group = Taxon)) +
                                geom_line(size = 2) +
                                geom_point(size = 4, shape = 21, fill = "white") + 
                                scale_colour_manual(values = cols) +
                                theme_classic() +
                                labs(title = paste0("JASPAR ", collection, " data growth per release"), y = "# Profiles", x = "Year of release")
                            
  ggsave(plot     = nb.motfs.per.release.gg,
         filename = file.path(results.dir, paste0("Jaspar_", collection, "_growth_lines.pdf")),
         width    = 9.5,
         height   = 7.5)
  
  nb.motfs.per.release.gg <- ggplotly(nb.motfs.per.release.gg,
                                      tooltip = c("y", "x", "group"))
  
  htmlwidgets::saveWidget(nb.motfs.per.release.gg, file.path(results.dir, paste0("Jaspar_", collection, "growth_lines.html")))
  message("; Line chart ready")
  
  
  
  ################
  ## Donut plot ##
  ################
  jaspar.donut <- nb.motfs.per.release %>% 
    dplyr::filter(Collection == collection) %>% 
    drop_na(Year, Nb_motifs) %>% 
    dplyr::filter(Year == release.year & Taxon != "All_taxa") %>% 
    ggplot(aes(x = 2, y = Nb_motifs, fill = Taxon)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = cols) +
    theme_void()
  
  ggsave(plot     = jaspar.donut,
         filename = file.path(results.dir, paste0("Jaspar_", collection, "_", release.year,"_pie-donut.pdf")),
         width    = 9.5,
         height   = 7.5)
  
  ## Generate the donut plot
  jaspar.donut <- 
    nb.motfs.per.release %>% 
    dplyr::filter(Collection == collection) %>% 
    drop_na(Year, Nb_motifs) %>% 
    dplyr::filter(Year == release.year & Taxon != "All_taxa") %>% 
    group_by(Taxon) %>% 
    plot_ly(labels = ~Taxon, values = ~Nb_motifs) %>%
    add_pie(hole = 0.6) %>%
    layout(title = paste0("JASPAR ", release.year, " ", collection),  showlegend = TRUE,
           xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE),
           yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE), colorway = tax.cols)
  
  htmlwidgets::saveWidget(jaspar.donut, file.path(results.dir, paste0("Jaspar_", collection, "_", release.year,"_pie-donut.html")))
  message("; Donut chart ready")
  
  
  ##############
  ## Bar plot ##
  ##############
  jaspar.bars <- nb.motfs.per.release %>% 
                    dplyr::filter(Taxon != "All_taxa") %>% 
                    dplyr::filter(Collection == collection) %>% 
                    drop_na(Year, Nb_motifs) %>% 
                 ggplot(aes(x=Year, y=Nb_motifs, fill=Taxon)) +
                    geom_bar(stat = "identity", position = "stack") +
                    scale_fill_manual(values = cols) +
                    theme_classic() +
                    theme(text = element_text(size = 15),
                          axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
                          axis.text.y = element_text(hjust = 1, size = 20)) +
                    labs(title = paste0("JASPAR ", collection, " data growth"), y = "# Profiles", x = "Year of release") +
                    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(plot     = jaspar.bars,
         filename = file.path(results.dir, paste0("Jaspar_", collection, "_", release.year,"_growth_barplot.pdf")),
         width    = 9.5,
         height   = 7.5)
  
  jaspar.bars <- ggplotly(jaspar.bars,
                          tooltip = c("y", "x", "fill"))
  
  
  htmlwidgets::saveWidget(jaspar.bars, file.path(results.dir, paste0("Jaspar_", collection, "_", release.year,"_growth_barplot.html")))
  message("; Bar chart ready")
}

















#########################################################
## If required, upload the plots at the plotly website ##
#########################################################
if (plotly.export) {
  
  ## Set environment variable for plotly API key
  Sys.setenv("plotly_username" = "jaimicore")
  Sys.setenv("plotly_api_key"  = "your_api_key")  ## https://plot.ly/settings/api
  
  # Create a shareable link to the charts
  # Set up API credentials: https://plot.ly/r/getting-started
  api_create(nb.motfs.per.release.gg, filename = "JASPAR_CORE_growth_line_plot")
  api_create(jaspar.donut, filename = paste0("Jaspar_", release.year, "_pie-donut"))
  api_create(jaspar.bars, filename = paste0("Jaspar_", release.year, "_bars"))
}