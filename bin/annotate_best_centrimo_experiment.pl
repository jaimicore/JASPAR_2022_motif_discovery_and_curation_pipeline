#!/usr/bin/env perl
#################################################################
## Authors:                                                    ##
##    - Robin van der Lee                                      ##
##        - robinvanderlee AT gmail DOT com                    ##
##                                                             ##
##    - Jaime A Castro-Mondragon                               ##
##        - j DOT a DOT c DOT mondragon AT ncmm DOT uio DOT no ##
#################################################################

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use File::Spec;


### MAKE SURE THE EXPERIMENT_MAP HAS:
my $EXP_ID_FIELD = 0;
my $EXP_TF_FIELD = 2;



### PARSE OPTIONS
my $tf_jaspar_map = ""; # tf_jaspar_map = config["TF_Experiment_map"],
my $best_exp = ""; # best_exp = os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best")
my $output = ""; # os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best.TF_associated")

my $USAGE = "Usage: $0 --best BEST_EXP_FILE --map TF_JASPAR_MAP_FILE --output OUTPUT_FILE\n";
if( scalar @ARGV != 6 ){
	die $USAGE;
}
GetOptions( 'best=s' => \$best_exp,
			'map=s' => \$tf_jaspar_map,
			'output=s' => \$output) or die $USAGE;

# print $tf_jaspar_map . "\n";
# print $best_exp . "\n";
# print $output . "\n";



############
## STEP 1 ##
############
my $motif = "";
my $motif_folder  = "";

my $centrimo_best_file_name = "";
my $centrimo_best_file_name__basename = "";
my $centrimo_best_pvalue  = "";
my $PDF_motif_info = "";

open(BEST_EXP_FH, "<$best_exp") or die "File $best_exp does not exist";
while(<BEST_EXP_FH>){
	chomp;
	
	my @F = split /\t/;
	$centrimo_best_file_name = $F[0];
	$centrimo_best_pvalue = "$F[1]";
	
	$centrimo_best_file_name__basename = basename($centrimo_best_file_name);
	my $centrimo_best_file_name__dirname  = dirname($centrimo_best_file_name);
	
	# Get first part of centrimo_best_file_name and rename it, e.g.
	# From : dm_Stat92E-GFP_Stat92E_E12-24_m1
	# To   : dm_Stat92E-GFP_Stat92E_E12-24_peak-motifs_m1
	#
	# ENCSR000BST_GATA3_MCF-7_m1.501bp.fa.sites.centrimo
	$centrimo_best_file_name__basename =~ /^(\S+)(\.\d+bp\.fa\.sites\.centrimo)$/; # get the file extension, tricky because some worm ID contain a period, e.g. ce_OP532_Y116A8C.19_L1_WA_m5.501bp.fa.sites.centrimo	$1 will contain "ce_OP532_Y116A8C.19_L1_WA_m5", $2 will contain ".501bp.fa.sites.centrimo"
	$motif = $1;
	my $TF_exp_name = $1;
	my $extension = $2;
	$motif =~ s/_m(\d+)$/_peak-motifs_m$1/g;

	## PDF with selected motif
	$PDF_motif_info = $centrimo_best_file_name__dirname;
	$TF_exp_name =~ s/_m\d+$//gi;
	$PDF_motif_info = $PDF_motif_info."/selected_motif/".$TF_exp_name.$extension.".best.TF_associated.pdf";

	# rename the path name, e.g.
	# from /home/rvdlee/JASPAR/ModERN/results/ModERN_fly/output/dm_Stat92E-GFP_Stat92E_E12-24/central_enrichment
	# to /home/rvdlee/JASPAR/ModERN/results/ModERN_fly/output/dm_Stat92E-GFP_Stat92E_E12-24/motifs/jaspar/logos
	$motif_folder = $centrimo_best_file_name__dirname;
	$motif_folder =~ s/central_enrichment/motifs\/jaspar\/logos/gi;

	#print $motif . "\n";
	#print $motif_folder . "\n";
}
close(BEST_EXP_FH);


############
## STEP 2 ##
############

## Get all the parts required for the output

# find the information corresponding to the best_exp in the tf_jaspar_map
my $exp_id = "";
my $exp_tf = "";
open(TF_JASPAR_MAP_FH, "<$tf_jaspar_map");
while(<TF_JASPAR_MAP_FH>){
	chomp;
	my @F = split /\t/;

	# the experiment ID is in the first field of the map file
	my $exp_id_current = $F[$EXP_ID_FIELD];

	# if the current line of the map file corresponds to the experiment ID in the centrimo best file, then we've found the right information, and we can save it.
	if(begins_with($centrimo_best_file_name__basename, $exp_id_current)){
		$exp_id = $exp_id_current;

		## The TF name must be in the third field [0,1,2]
		$exp_tf = $F[$EXP_TF_FIELD];
	}
}
close(TF_JASPAR_MAP_FH);

my $logo_png_filepath = File::Spec->catfile($motif_folder, $motif . "_logo.png");

## create the output file
# open(OUTFILE_FG, ">$output");
# print OUTFILE_FG $exp_id . "\t" . $exp_tf . "\t" . $exp_tf . "\t" . $centrimo_best_file_name . "\t" . $centrimo_best_pvalue . "\t" . $logo_png_filepath . "\t" . $PDF_motif_info . "\n";
# close(OUTFILE_FG);

################################################################
## Print table with ready to curate
##
## FIELDS:
## 1) PWM path file
## 2) current_BASE_ID
## 3) current_VERSION
## 4) TF NAME
## 5) Uniprot
## 6) TAX_ID
## 7) class
## 8) family
## 9) TFBSshapeID
## 10) Data
## 11) Source
## 12) Validation
## 13) Comment	Addition or Upgrade or Non-validated (A or U or N )
## 14) BED
## 15) FASTA
## 16) Centrality pval
## 17) Centrality plot
## 18) Motif logo path
## 19) Experiment ID
## 20) Centrality p-value
## 21) PDF motif info

$centrimo_best_file_name__basename =~ /^(\S+)(\.\d+bp\.fa\.sites\.centrimo)$/; # get the file extension, tricky because some worm ID contain a period, e.g. ce_OP532_Y116A8C.19_L1_WA_m5.501bp.fa.sites.centrimo	$1 will contain "ce_OP532_Y116A8C.19_L1_WA_m5", $2 will contain ".501bp.fa.sites.centrimo"
$motif = $1;
my $TF_exp_name = $1;
$TF_exp_name = $exp_tf;
$motif =~ s/_m(\d+)$/_peak-motifs_m$1/g;

my $centrimo_best_file_pdf = $centrimo_best_file_name.".pdf";

$motif_folder = dirname($centrimo_best_file_name);
$motif_folder =~ s/central_enrichment/motifs\/jaspar\/pfm/gi;
my $pwm_filepath = File::Spec->catfile($motif_folder, $motif.".jaspar");

$motif_folder = dirname($centrimo_best_file_name);
$motif_folder =~ s/central_enrichment/matrix_sites/gi;
my $bed_filepath = File::Spec->catfile($motif_folder, $motif.".tf.sites.bed");
my $fasta_filepath = File::Spec->catfile($motif_folder, $motif.".tf.sites.fasta");


open(OUTFILE_FG, ">$output");
print OUTFILE_FG $pwm_filepath . "\t\t\t" . $exp_tf . "\t\t\t\t\t\t" . "ChIP-seq" . "\t\t\t\t" . $bed_filepath . "\t" . $fasta_filepath . "\t" . $centrimo_best_file_name . "\t" . $centrimo_best_file_pdf . "\t" . $logo_png_filepath . "\t" . $exp_id . "\t" . $centrimo_best_pvalue . "\t" . $PDF_motif_info . "\n"; 
close(OUTFILE_FG);

#################
### FUNCTIONS ###
#################

# check if the string in the first agrument starts with the string in the second argument
sub begins_with
{
    return substr($_[0], 0, length($_[1])) eq $_[1];
}
