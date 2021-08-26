#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "ggplot2",
                        "optparse",
                        "plotly",
                        "RColorBrewer",
                        "reshape2")

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
# results.dir <- "/home/jamondra/Downloads/JASPAR_2020_plots"


##########################################
## Read JASPAR motifs per release table ##
##########################################
motifs.per.taxon.tab <- fread(motifs.per.taxon.tab.file)

## Add column with the total number of motifs
motifs.per.taxon.tab$All_taxa <- rowSums(motifs.per.taxon.tab[,c("Vertebrates", "Plants", "Insects", "Nematodes", "Fungi", "Urochordates")],na.rm = T)

## Get the year of the latest release
release.year <- max(motifs.per.taxon.tab$Year)

## Convert the matrix in a data.frame
nb.motfs.per.release <- data.frame(melt(motifs.per.taxon.tab[,c(3:9)]))
colnames(nb.motfs.per.release) <- c("Taxon", "Nb_motifs")

## Add the year column
nb.motfs.per.release$Year <- rep(motifs.per.taxon.tab$Year, 7)
nb.motfs.per.release$Year <- ordered(nb.motfs.per.release$Year, levels = as.vector(unique(nb.motfs.per.release$Year)))

## Order factors by number of motifs
nb.motfs.per.release$Taxon <- factor(nb.motfs.per.release$Taxon, levels = as.vector(c("All_taxa", "Vertebrates", "Plants", "Fungi", "Insects", "Nematodes", "Urochordates")))

## Define color palette
taxon.colors <- colorRampPalette(brewer.pal(7, "Dark2"), space = "Lab")(length(unique(nb.motfs.per.release$Taxon)))
cols <- c( "#666666", taxon.colors)
names(cols) <- c("All_taxa", "Vertebrates", "Plants", "Fungi", "Insects", "Nematodes", "Urochordates")


################
## Line chart ##
################
nb.motfs.per.release.gg <- 
  ggplot(nb.motfs.per.release, aes(x = Year, y = Nb_motifs, color = Taxon, group = Taxon)) +
  geom_line(size = 2) +
  geom_point(size = 4, shape = 21, fill = "white") + 
  # scale_color_brewer(palette = "Dark2") +
  scale_colour_manual(values = cols) +
  theme_classic() +
  labs(title = "JASPAR CORE data growth per release", y = "# Profiles", x = "Year of release")

ggsave(plot     = nb.motfs.per.release.gg,
       filename = file.path(results.dir, "Jaspar_growth_lines.pdf"),
       width = 9.5, height = 7.5)

nb.motfs.per.release.gg <- ggplotly(nb.motfs.per.release.gg,
                                    tooltip = c("y", "x", "group"))

htmlwidgets::saveWidget(nb.motfs.per.release.gg, file.path(results.dir, "Jaspar_growth_lines.html"))




################
## Donut plot ##
################
jaspar.donut <- nb.motfs.per.release %>% 
  dplyr::filter(Year == release.year & Taxon != "All_taxa") %>% 
  ggplot(aes(x = 2, y = Nb_motifs, fill = Taxon)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  xlim(0.5, 2.5) +
  scale_fill_brewer(palette = "Dark2") +
  theme_void()

ggsave(plot     = jaspar.donut,
       filename = file.path(results.dir, paste0("Jaspar_", release.year,"_pie-donut.pdf")),
       width    = 9.5,
       height   = 7.5)

## Generate the donut plot
jaspar.donut <- 
  nb.motfs.per.release %>% 
  dplyr::filter(Year == release.year & Taxon != "All_taxa") %>% 
  group_by(Taxon) %>% 
  plot_ly(labels = ~Taxon, values = ~Nb_motifs) %>%
  add_pie(hole = 0.6) %>%
  layout(title = paste0("JASPAR ", release.year, " CORE profile stats"),  showlegend = TRUE,
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE), colorway = taxon.colors)

htmlwidgets::saveWidget(jaspar.donut, file.path(results.dir, "Jaspar_growth_lines.html"))



##############
## Bar plot ##
##############
jaspar.bars <- nb.motfs.per.release %>% 
  dplyr::filter(Taxon != "All_taxa") %>% 
  ggplot(aes(x=Year, y=Nb_motifs, fill=Taxon)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
          axis.text.y = element_text(hjust = 1, size = 20)) +
    labs(title = "JASPAR CORE data growth", y = "# Profiles", x = "Year of release") +
    theme(plot.title = element_text(hjust = 0.5))


ggsave(plot     = jaspar.bars,
       filename = file.path(results.dir, paste0("Jaspar_", release.year,"_growth_barplot.pdf")),
       width    = 9.5,
       height   = 7.5)



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