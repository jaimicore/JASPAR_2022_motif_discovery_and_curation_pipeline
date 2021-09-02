#!/usr/bin/Rscript

# Import command line parser
suppressMessages(library(optparse))
suppressMessages(library(dplyr))

###############################################################################################
# Parse
###############################################################################################
# Declare parsing options
option_list = list( make_option( c("-a", "--annotation"), type="character",  default=NULL, 
                                 help=paste0("Mandatory option. Name of the file containing the annotations to be added to the radial tree.\n",
                                             "\t\tThe file must have the following column content without column names and should be tab separated:\n",
                                             "\t\t\t1\tName of the matrix collection\t:\tThe name should be identical to the matrix collection name used in the RSAT::matrix-clustering command.\n",
                                             "\t\t\t2\tID of the matrix\t:\tThe matrix ID should be identical to the one used in the collection as input in RSAT::matrix-clustering command.\n",
                                             "\t\t\t3\tColor code\t:\tHexadecimal color code associated to each matrix class.The color must be always provided betwen double quotes.\n",
                                             "\t\t\t4\tColor code ID\t:\tA unique integer number associated to each matrix class.\n",
                                             "\t\tPlease ensure in all moment that no special characters from regular expressions are included in the name."), metavar="character"),
                    make_option( c("-i", "--input"),      type="character",  default=NULL, 
                                 help=paste0("Mandatory option. Name of the output directory produced by RSAT::matrix-clustering containing the radial tree.\n",
                                             "\t\tPlease ensure in all moment that no special characters from regular expressions are included in the name."), metavar="character") );

# Parse options
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# Mandatory options # 1
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied to --input directory_matrix-clustering_radial_tree", call.=FALSE)
}
# Mandatory options # 2
if (is.null(opt$annotation)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied to --annotation file.tab", call.=FALSE)
}

# Print command
cmd_annotate_matrix_clustering <- paste0("Rscript annotate_matrix-clustering.R --input=",opt[["input"]],
                                          " --annotation=",opt[["annotation"]])
print(cmd_annotate_matrix_clustering)

###############################################################################################
# Import all Libraries
###############################################################################################
suppressMessages(library(jsonlite))

###############################################################################################
#                                            MAIN                                             #
###############################################################################################


#########################################################
# Parse the JSON containing the radial tree 

# Strip absolute path from directory
prefix.files <- list.files(opt[["input"]], pattern = "SUMMARY.html")
prefix.files <- gsub(prefix.files, pattern = "_SUMMARY.html", replacement = "")

# Create name of output tree subdirectory
#tree_subdir_name          <- paste0(opt[["input"]],"/", tree_dir_name,"_trees" )
tree_subdir_name          <- paste0(opt[["input"]],"/", prefix.files,"_trees" )
# Create output name for matrix_order.txt
matrix_order.file         <- paste0(tree_subdir_name,"/matrix_order.txt")
# Create output name for matrix_cluster_color.txt
matrix_cluster_color.file <- paste0(tree_subdir_name,"/matrix_cluster_color.txt")
# Create filename for JSON file
tree_json.file            <- paste0(opt[["input"]],"/",prefix.files,"_trees/parsed_tree_cluster_1.json")
# Test that file exists
if(!file.exists(tree_json.file)){
  print(paste0("; ERROR: Unable to find file ",tree_json.file, " !"))
  stop("Unable to proceed")
}

# Create cmd for parsing json and extracting the matrix order
cmd_grep_json_order <- paste0("grep 'label' ",tree_json.file,
                        " | awk '{print $2}'| perl -pe 's/\\\"|\\,//g' > ", 
                        matrix_order.file )

# Create cmd for parsing json and extracting the cluster colors
cmd_grep_clust_color <- paste0("grep -A 11 'label' ", tree_json.file,
                               " | grep 'branch_color' | awk '{print $3}' |",
                               "perl -pe 's/,//g' > ",
                               matrix_cluster_color.file)

# Execute it!
print(cmd_grep_json_order)
system(cmd_grep_json_order)

print(cmd_grep_clust_color)
system(cmd_grep_clust_color)

#########################################################
# Read the matrix order and clusters colors
matrix_order.df  <- read.table(file  = matrix_order.file,         header = FALSE, sep = "\t")
print(paste0("; INFO Successful read of file containing matrix information ",matrix_order.file))
cluster_color.df <- read.table(file  = matrix_cluster_color.file, header = FALSE, sep = "\t")
print(paste0("; INFO Successful read of file containing matrix cluster color information ",matrix_cluster_color.file))

# Combine columns
matrix_order.df <- cbind(matrix_order.df, cluster_color.df)
# Name columns
colnames(matrix_order.df) <- c("matrix_name", "cluster_nb")
# Add number of position
matrix_order.df$id_motif <- 1:dim(matrix_order.df)[1]

#########################################################
# Merge  annotation matrix with matrix order
annotation.df           <- read.table(file  = opt[["annotation"]], header = FALSE, sep = "\t")
print(paste0("; INFO Successful read of file containing matrix annotation information ",opt[["annotation"]]))


# Name columns
colnames(annotation.df) <- c("collection_name", "matrix_name", "class", "class_nb")

# Unquote columns
annotation.df$collection_name <- gsub("\"","",as.character(annotation.df$collection_name))
annotation.df$matrix_name     <- gsub("\"","",as.character(annotation.df$matrix_name))
# annotation.df$collection_name <- gsub("\"","",as.character(annotation.df$collection_name))

print("; INFO Merging color code with  matrix IDs sorted by their position in the tree.")

#NOTE WSG I need to add inside the apply a sanity check: replace all regex special characters with "\\specialchar"
# Merge data frames
annotation.merge.df <- apply(annotation.df, 1 ,function(matrix_element){
  # Create regex 
  regex_string    <- paste0("^",matrix_element["collection_name"],"_m\\d+_",matrix_element["matrix_name"],"$")
  # NOTE WSG: spcial character from regex must be forbidden in 
  # the file names. Need to add a sanity check
  # Safe check for special characters in regex
  #regex_string <- sub("\\+", "\\\\+", regex_string) 
  # Find complete motif name
  match_matrix.df <- matrix_order.df[(grepl(regex_string,as.character(matrix_order.df$matrix_name))),]
  # Transpose match 
  matrix_element <- as.data.frame(t(as.data.frame(matrix_element)))
  # Add annotation columns to matched matrix
  match_matrix.df <- cbind(match_matrix.df,matrix_element)

  return(match_matrix.df)
}) %>% (plyr::rbind.fill)

print("; INFO Calculating degrees intervals for each annotation layer.")

# Sort the data frame by id_motif
annotation.merge.df       <- annotation.merge.df[order(annotation.merge.df$id_motif),]
annotation.merge.df$start <- (0:(dim(annotation.merge.df)[1] - 1)) * (360/dim(annotation.merge.df)[1])
annotation.merge.df$end   <- (1:(dim(annotation.merge.df)[1])) * (360/dim(annotation.merge.df)[1])

#########################################################
# Convert to JSON
print("; INFO Converting data.frame to JSON format.")

annotation.json <- toJSON(annotation.merge.df)
#########################################################
# Write JSON file to matrix-clustering directory folder
print("; INFO Writing JSON file")
annotation_json_file <- paste0(tree_subdir_name,"/annotation_matrix.json")
write_json(annotation.json,annotation_json_file)

#########################################################
# Modify radial tree html result
print("; INFO Creating annotated radial tree html.")
# Create python script for radial tree
cmd_annotate_htmltree <- paste0("python annotate-html-radialtree.py ",
                                "-i ", file.path(opt[["input"]], prefix.files))

#Execute it!
print(cmd_annotate_htmltree)
system(cmd_annotate_htmltree)

print("Succesful annotation!\n\n")

            

