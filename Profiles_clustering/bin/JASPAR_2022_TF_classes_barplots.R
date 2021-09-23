#############################
## Load required libraries ##
#############################
required.packages = c("dplyr",
                      "data.table",
                      "ggplot2",
                      "optparse",
                      "plotly",
                      "rcartocolor",
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

  make_option(c("-s", "--suffix"), type = "character", default = "example_collection", 
              help = "Suffix to be added to the exported table file names", metavar = "character"),
  
  make_option(c("-o", "--output_directory"), type = "character", default = 0,
              help = "Output directory where the resulting tables will be exported. (Mandatory) ", metavar = "character"),
  
  make_option(c("-c", "--collection"), type = "character", default = NULL,
              help = "JASPAR motif collection. [CORE | UNVALIDATED]", metavar = "character")
);
message("; Reading arguments from command-line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


## Set variable names
metadata.tab.file   <- opt$metadata_table
out.dir             <- opt$output_directory
file.suffix         <- opt$suffix
collection.type     <- opt$collection


#################################################
## Set environment variable for plotly API key ##
#################################################
# Sys.setenv("plotly_username" = "jaimicore")
# Sys.setenv("plotly_api_key"  = "your_api_key")  ## https://plot.ly/settings/api


###########
## Debug ##
###########
# out.dir             <- "/home/jamondra/Downloads/TF_Classes_barplots"

# metadata.tab.file   <- "/home/jamondra/Downloads/JASPAR_2022_metadata_nonredundant_CORE_urochordates_table_parsed.tab"
# collection.type     <- "C"
# file.suffix         <- "JASPAR_2022_nonredundant_CORE_nematodes"

# metadata.tab.file   <- "/home/jamondra/Downloads/JASPAR_2022_metadata_nonredundant_UNVALIDATED_vertebrates_table_parsed.tab"
# collection.type     <- "U"
# file.suffix         <- "JASPAR_2022_nonredundant_UNVALIDATED_vertebrates"



###########################
## Create output folders ##
###########################
dir.create(out.dir, recursive = T)


#########################
## Read metadata table ##
#########################
message("; Reading metadata table: ", metadata.tab.file)
metadata <- fread(metadata.tab.file, header = T) %>% 
              mutate(ID = paste0(base_id, ".", version)) %>%
              mutate(collection = ifelse(grepl(pattern = "^UN", x = ID), yes = "UNVALIDATED", no = "CORE"),
                     class      = ifelse(class == "", yes = "Unknown", no = class)) %>% 
              select(name, class, ID, collection) %>% 
              separate_rows(name, class, sep = "::") %>% 
              distinct() %>% 
              dplyr::filter(name != "name")


##################################
## TF class - Color assignation ##
##################################
message("; Assigning colours to TF classes")

## Get the TF class names ordered by number of motifs
TF.class.order <- metadata %>% 
                    group_by(class) %>% 
                    tally() %>% 
                    arrange(desc(n)) %>% 
                    dplyr::filter(class != "Unknown")
TF.class.order <- TF.class.order$class

TF.known.classes    <- TF.class.order
TF.known.classes.nb <- length(TF.known.classes)

## Grey color for Unknown class
unknown.class.color <- "#888888"

nb.classes.palette <- 12  ## We want to use the first 12 colors of the Safe palette
nb.seed.colors     <- ifelse(TF.known.classes.nb < nb.classes.palette,
                             yes = TF.known.classes.nb,
                             no = nb.classes.palette)

## Generate a carto palette (remove gray value)
carto.pal.classes  <- carto_pal(nb.seed.colors, "Safe")
carto.pal.classes  <- carto.pal.classes[which(carto.pal.classes != "#888888")]


## Expand the color palette and add the gray color at the end
class.colors        <- c(colorRampPalette(carto.pal.classes, space = "Lab")(TF.known.classes.nb),
                         unknown.class.color)
names(class.colors) <- c(TF.known.classes, "Unknown")


df.class.colour <- data.frame(colour   = class.colors,
                              class    = names(class.colors),
                              class_nb = seq_len(TF.known.classes.nb + 1))
class.colors    <- rev(class.colors)


#################################################
## Prepare metadata table to be read by ggplot ##
#################################################
metadata.w.colors       <- merge(metadata, df.class.colour, by = "class") %>% 
                              arrange(class_nb)
metadata.w.colors$class <- factor(metadata.w.colors$class, levels = rev(c(TF.class.order, "Unknown")))


## Adapt to collection type
metadata.w.colors$collection <- factor(metadata.w.colors$collection, levels = rev(c("CORE", "UNVALIDATED")))


## This dataframe contains the labels to be displayed on each bar
## Adapt to collection type
label.df <- metadata.w.colors %>% 
              group_by(class) %>% 
              mutate(Total_Class = n()) %>%
              group_by(class, collection) %>% 
              mutate(COREc = sum(collection %in% "CORE"),
                     UNVc  = sum(collection %in% "UNVALIDATED")) %>% 
              select(class, COREc, UNVc) %>% 
              distinct() %>% 
              ungroup() %>% 
              group_by(class) %>% 
              mutate(CORE        = max(COREc),
                     UNVALIDATED = max(UNVc)) %>% 
              select(class, CORE, UNVALIDATED) %>% 
              distinct()
          

if (collection.type == "UNVALIDATED") {
  
  label.df <- label.df %>% 
                mutate(lab         = paste0(CORE, "/", UNVALIDATED),
                       Total_Class = sum(c(CORE, UNVALIDATED))) %>% 
                data.table()
  
} else if (collection.type == "CORE") {
  
  label.df <- label.df %>% 
                mutate(lab         = CORE,
                       Total_Class = sum(c(CORE, UNVALIDATED))) %>% 
                data.table()
}

  
  

## Use this variable to have enough space to insert the labels
max.y <- max(label.df$Total_Class) + 25

## The alpha parameter and plot title change depending on the collection type
alpha.values  <- NULL
y.axis.lab    <- NULL
barplot.title <- NULL
if (collection.type == "UNVALIDATED") {
  
  alpha.values  <- c(0.4, 1)
  y.axis.lab    <- "Number of motifs (CORE/UNVALIDATED)"
  barplot.title <- gsub(file.suffix, pattern = "_UNVALIDATED$", replacement = "_CORE+UNVALIDATED")
  
} else if (collection.type == "CORE") {
  
  alpha.values  <- 1
  y.axis.lab    <- "Number of motifs (CORE)"
  barplot.title <- file.suffix
}

## Barplot code
TF.class.barplot <- ggplot(metadata.w.colors, aes(alpha = collection, x = class, fill = class)) +
                      geom_bar() + 
                      coord_flip() +
                      scale_fill_manual(values = class.colors) +
                      theme_classic() +
                      labs(x = "TF class", y = y.axis.lab, title = barplot.title) +
                      scale_alpha_manual(values = alpha.values) +
                      theme(text        = element_text(size = 15),
                            axis.text.y = element_text(angle = 0, hjust = 1, size = 10),
                            axis.text.x = element_text(hjust = 0.5, size = 15),
                            legend.position = "none",
                            legend.title    = element_blank(),
                            legend.text     = element_text(size = 9),
                            legend.box      = "vertical",
                            plot.title      = element_text(hjust = 0.5),
                            panel.grid.major.x = element_line(color = "#969696",
                                                              size = 0.25,
                                                              linetype = 2)) +
                      geom_text(data = label.df, aes(x = class, y = Total_Class, label = lab), vjust = 0.5, hjust = -0.2, inherit.aes = F, size = 3.5) +
                      guides(fill = guide_legend(reverse = T))  +
                      scale_y_continuous(limits = c(0, max.y), expand = c(0,2), breaks = seq(0, max.y, by = 50)[-1])


## Convert ggplot to plotly
TF.class.barplotly <- ggplotly(TF.class.barplot,
                               tooltip = c("alpha", "x", "y")) %>% 
                      config(displayModeBar = F)


###################################################
## Export barplot: interactive and static format ##
###################################################
TF.class.barplotly.html <- file.path(out.dir, paste0("TF_classes_barplot_", file.suffix, ".html"))
htmlwidgets::saveWidget(widget = TF.class.barplotly,
                        file   = TF.class.barplotly.html, selfcontained = F)
message("; Interactive barplot created: ", TF.class.barplotly.html)


# Create a shareable link to the charts
# Set up API credentials: https://plot.ly/r/getting-started
# api_create(TF.class.barplotly, filename = paste0("TF_classes_barplot_", file.suffix, ".html"), )
# message("; Interactive barplot created online, available through plotly account")


TF.class.barplot.img <- file.path(out.dir, paste0("TF_classes_barplot_", file.suffix, ".jpeg"))
ggsave(plot     = TF.class.barplot,
       filename = TF.class.barplot.img,
       width    = 12,
       height   = 8)
message("; JPEG barplot created: ", TF.class.barplot.img)