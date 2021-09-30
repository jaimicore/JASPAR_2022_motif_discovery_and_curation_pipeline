################################################################
## Import functions
import os
import itertools
import subprocess
from snakemake.utils import R

## git push origin HEAD:master



#####################################################################
## Select config file according to the given command line argument ##
## snakemake --config analysis_id=REMAP2022_Human                  ##
#####################################################################
if config['analysis_id'] == "ReMap2020_Human":
    configfile: "config_files/config_Hsapiens_hg38.yaml"
elif config['analysis_id'] == "ChExMix_sacCer3":
    configfile: "config_files/config_yeast_sacCer3_chexmix.yaml"
elif config['analysis_id'] == "ReMap2020_Athaliana":
    configfile: "config_files/config_Athaliana_TAR10_REMAP.yaml"
else :
    print("; ERROR: analysis name not found. Supported: ReMap2020_Human | ReMap2020_Athaliana | ChExMix_sacCer3")


##########################
## Initialize variables ##
##########################

## Transcription Factor names
(FOLDERS, TF_NAMES) = glob_wildcards(os.path.join(config["data_folder"], "{Folders}", "{TF}_peaks.narrowPeak"))


## Peak summit extension (both sides)
PEAK_LENGTH = "101 501".split()
PEAK_EXTENSIONS =  {
    '101' : '50',
    '501' : '250'
}


## This is done to launch the rule 'Concat_annotated_experiments'
## It expects as input all logo files, so the rule is launched only
## once the previous rule is completed for all datasets
MOST_ENRICHED_MOTIF_ASSOC_LOGO = expand(os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best.TF_associated.tex"), TF = TF_NAMES)

## JASPAR 2020 annotation table
## Output from rule JASPAR_annotation_table
JASPAR_ANN_TAB = os.path.join(config["out_dir"], "Jaspar_2020_info_" + config["taxon"] + "_table.tab")

LOGPVAL = [config['central_pvalue']]

################################################################
## Rules
rule all:
    input:
        expand(os.path.join(config["out_dir"], "{TF}", "peak-motifs", "results", "discovered_motifs", "{TF}_motifs_discovered.tf"), TF = TF_NAMES), \
        expand(os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best.TF_associated"), TF = TF_NAMES), \
        JASPAR_ANN_TAB, \
        expand(os.path.join(config["curation_dir"], "Renamed_log_pval_{logpval}.txt"), logpval = LOGPVAL), \
	expand(os.path.join(config["curation_dir"], "Selected_motifs_to_curate_log10_pval_{logpval}.tab"), logpval = LOGPVAL), \
        expand(os.path.join(config["curation_dir"], "Selected_motifs_to_curate_log10_pval_{logpval}.pdf"), logpval = LOGPVAL)

        

################################################################
################################################################

#####################################################
## Prepare BED file with summit coordinates        ##
## This rule may vary according to the data source ##
#####################################################

#############################
## REMAP 2022 Homo sapiens ##
#############################
if config['analysis_id'] == "ReMap2020_Human":
    
    rule extract_peak_summits:
        """
        Extract peak summits: chromosome, start, end
        """
        input:
            os.path.join(config["data_folder"], "{TF}", "{TF}_peaks.narrowPeak")
        output:
            os.path.join(config["out_dir"], "{TF}", "peak_summits", "{TF}_peak_summits.bed")
        message:
            "; Peak summits - TF : {wildcards.TF}"
        priority:
            100
        shell:
            """
            awk '{{ print $1"\\t"($2+$10)"\\t"($2+$10+1)}}' {input} > {output}
            """
            
#####################################
## REMAP 2022 Arabidopsis thaliana ##
#####################################
if config['analysis_id'] == "ReMap2020_Athaliana":
   
   rule extract_peak_summits:
        """
        Extract peak summits: chromosome, start, end
        """
        input:
            os.path.join(config["data_folder"], "{TF}", "{TF}_peaks.narrowPeak")
        output:
            os.path.join(config["out_dir"], "{TF}", "peak_summits", "{TF}_peak_summits.bed")
        message:
            "; Peak summits - TF : {wildcards.TF}"
        priority:
            100
        shell:
            """
            awk '{{ print "chr"$1"\\t"($2+$10)"\\t"($2+$10+1)}}' {input} > {output}
            """

            
#############
## ChExMix ##
#############       
elif config['analysis_id'] == "ChExMix_sacCer3":

        rule extract_peak_summits:
            """
            Extract peak summits: chromosome, start, end
            """
            input:
                os.path.join(config["data_folder"], "{TF}", "{TF}_peaks.narrowPeak")
            output:
                os.path.join(config["out_dir"], "{TF}", "peak_summits", "{TF}_peak_summits.bed")
            message:
                "; Peak summits - TF : {wildcards.TF}"
            priority:
                100
            shell:
                """
                awk '{{ print $1"\\t"$2"\\t"$3}}' {input} > {output}
                """


################################################################
################################################################

rule extend_ENCODE_peak_files:
    """
    Extend the peak summits to both sides.
    """
    input:
        os.path.join(config["out_dir"], "{TF}", "peak_summits", "{TF}_peak_summits.bed")
    output:
        os.path.join(config["out_dir"], "{TF}", "extended_peaks", "{TF}.{length}bp.bed")
    message:
        "; Extending peak summits - TF : {wildcards.TF} - Final length: {wildcards.length}"
    params:
        genome = config["genome_size"],
        slop = lambda wildcards: PEAK_EXTENSIONS[wildcards.length]
    priority:
        99
    shell:
        """
        bedtools slop -b {params.slop} -i {input} -g {params.genome} > {output}
        """


rule get_fasta_ENCODE_peak_files:
    """
    Get the fasta sequences for the extended peaks
    """
    input:
        os.path.join(config["out_dir"], "{TF}", "extended_peaks", "{TF}.{length}bp.bed")
    output:
        os.path.join(config["out_dir"], "{TF}", "fasta", "{TF}.{length}bp.fa")
    message:
        "; Retrieving FASTA sequences from extended peaks - TF : {wildcards.TF} - Length: {wildcards.length}"
    params:
        genome = config["genome_fasta"]
    priority:
        98
    shell:
        """
        bedtools getfasta -fi {params.genome} -bed {input} -fo {output}
        """


###################
## RSAT analysis ##
###################
rule RSAT_peakmotifs_per_exp:
    """
    Run RSAT peak-motifs (motif discovery with different algorithms) on every peakset.
    The set of N discovered motifs is stored in a single transfac file.
    """
    input:
        os.path.join(config["out_dir"], "{TF}", "fasta", "{TF}.101bp.fa")
    output:
        os.path.join(config["out_dir"], "{TF}", "peak-motifs", "results", "discovered_motifs", "{TF}_motifs_discovered.tf")
    message:
        "; Running RSAT peak-motifs - TF : {wildcards.TF} "
    priority:
        97
    params:
        RSAT = config["RSAT"],
        disco = config["peakmotifs_disco"],
        nb_motifs = config["peakmotifs_disco_nb_motifs"],
        min_oligo = config["peakmotifs_disco_min_oligo"],
        max_oligo = config["peakmotifs_disco_max_oligo"],
        ci = config["peakmotifs_class_interval"],
        cisbp = os.path.join(config["RSAT"], "public_html/motif_databases/cisBP/cisBP_Homo_sapiens_2014-10.tf"),
        jaspar_motifs = os.path.join(config["RSAT"], "public_html/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.tf"),
        #jaspar_motifs = os.path.join(config["RSAT"], "public_html/motif_databases/JASPAR/Jaspar_2020/nonredundant/JASPAR2020_CORE_vertebrates_non-redundant_pfms.tf"),
        #task = "purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan",
        task = "purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis",
        prefix = "{TF}",
        peakmo_outdir = os.path.join(config["out_dir"], "{TF}", "peak-motifs")
    shell:
        """
        {params.RSAT}/perl-scripts/peak-motifs -v 2 \
        -r_plot \
        -title {wildcards.TF} \
        -i {input} \
        -markov auto \
        -disco {params.disco} \
        -nmotifs {params.nb_motifs} \
        -minol {params.min_oligo} \
        -maxol {params.max_oligo} \
        -no_merge_lengths \
        -ci {params.ci} \
        -noov \
        -2str \
        -origin center \
        -motif_db jaspar_vertebrates tf {params.jaspar_motifs} \
        -scan_markov 1 \
        -task {params.task} \
        -prefix {params.prefix} \
        -img_format png \
        -outdir {params.peakmo_outdir} ;
        """



##############################
## Matrix format conversion ##
##############################
rule Motif_number_to_ID:
    """
    Map table to associate motif numbers (m1, m2, ...) with motif names from peak-motifs (e.g., oligos_6nt_mkv3_m1)
    """
    input:
        os.path.join(config["out_dir"], "{TF}", "peak-motifs", "results", "discovered_motifs", "{TF}_motifs_discovered.tf")
    output:
       os.path.join(config["out_dir"], "{TF}", "peak-motifs", "results", "discovered_motifs", "{TF}_motifs_map.tab")
    message:
        "; Associating motif numbers with motif IDs : {wildcards.TF} "
    priority:
        96
    shell:
        """
        grep '^AC' {input} | cut -d' ' -f3 | perl -pe '$_ = "m$.\t$_"' > {output}

        """


     
checkpoint  RSAT_PSSM_to_JASPAR_format:
    """
    Convert the RSAT discovered matrices from transfac (tf) format to JASPAR format.
    The transfac file is split in several files (one per motif) adding as suffix the motif number (e.g., 1,2,3,...).
    """
    input:
        os.path.join(config["out_dir"], "{TF}", "peak-motifs", "results", "discovered_motifs", "{TF}_motifs_discovered.tf")
    output:
       directory(os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "pfm"))
    message:
        "; Converting RSAT motifs (tf format) to JASPAR format - TF : {wildcards.TF} "
    priority:
        95
    params:
        RSAT = config["RSAT"],
        return_fields = "counts",
        prefix = "peak-motifs",
        prefix_motif = os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "pfm", "{TF}"),
        out_dir = os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "pfm")
    shell:
        """
        mkdir -p {params.out_dir} ;
        {params.RSAT}/perl-scripts/convert-matrix -v 2 \
        -from tf -to jaspar \
        -i {input} \
        -return {params.return_fields} \
        -split \
        -prefix {params.prefix} \
        -o {params.prefix_motif} ;
        {params.RSAT}/perl-scripts/convert-matrix -v 2 \
        -from tf -to tab \
        -i {input} \
        -return {params.return_fields} \
        -split \
        -prefix {params.prefix} \
        -o {params.prefix_motif}
        """



rule JASPAR_PSSM_to_PWM:
    """
    Convert the JASPAR PSSMs in PWMs.
    This rule is executed for each discovered motif.
    """
    input:
        motif = os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "pfm", "{TF}_peak-motifs_m{n}.jaspar")
    output:
        os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "pwm", "{TF}_peak-motifs_m{n}.jaspar.pssm")
    message:
        "; Generating PWM from JASPAR matrices - TF : {wildcards.TF} "
    priority:
        94
    params:
        scripts_bin = config["bin"]
    shell:
        """
        perl {params.scripts_bin}/PCM_to_PWM.pl \
        -f {input.motif} > {output}
        """


rule Generate_matrix_logo:
    """
    Generate motif logos using RSAT convert-matrix
    This rule is executed for each discovered motif.
    """
    input:
        os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "pfm", "{TF}_peak-motifs_m{n}.jaspar")
    output:
        os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "logos", "{TF}_peak-motifs_m{n}_logo.png")
    message:
        "; Generating PWM from JASPAR matrices - TF : {wildcards.TF} "
    params:
        logo_dir = os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "logos"),
        logo_name = "{TF}_peak-motifs_m{n}",
        RSAT = config["RSAT"]
    priority:
        93
    shell:
        """
        {params.RSAT}/perl-scripts/convert-matrix -v 2 \
        -i {input} \
        -from jaspar -to jaspar \
        -return logo \
        -logo_dir {params.logo_dir} \
        -logo_no_title \
        -prefix {params.logo_name}
        """


############################################################
## Matrix scan: find sites to build the discovered motifs ##
############################################################
rule find_RSAT_matrix_sites:
    """
    Find the TFBSs used to built the matrices.
    This rule is executed for each discovered motif.
    """
    input:
        logos = os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "logos", "{TF}_peak-motifs_m{n}_logo.png"), \
        sequences = os.path.join(config["out_dir"], "{TF}", "peak-motifs", "data", "sequences", "{TF}_test.fasta"), \
        matrix = os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "pfm", "{TF}_peak-motifs_m{n}.tab"), \
        bg_file = os.path.join(config["out_dir"], "{TF}", "peak-motifs", "results", "composition", "{TF}_test_inclusive-1str-ovlp_2nt.txt")
    output:
        os.path.join(config["out_dir"], "{TF}", "matrix_sites", "{TF}_peak-motifs_m{n}.tf.sites")
    message:
        "; Scanning PSSM on 101bp peaks - TF : {wildcards.TF} - Matrix number: {wildcards.n}"
    params:
        RSAT = config["RSAT"]
    priority:
        92
    shell:
        """
        {params.RSAT}/bin/matrix-scan-quick \
        -i {input.sequences} \
        -m {input.matrix} \
        -bgfile {input.bg_file} \
        -t 5 \
        -return sites > {output}
        """

#/lsc/rsat/perl-scripts/convert-matrix -v 2         -from tf -to jaspar         -i /storage/scratch/JASPAR_2022/REMAP_2020_Athaliana/results/GSE60141_WRKY30_Col-0-leaves-tnt-colamp/peak-motifs/results/discovered_motifs/GSE60141_WR\KY30_Col-0-leaves-tnt-colamp_motifs_discovered.tf         -return counts         -split         -prefix peak-motifs         -o /storage/scratch/JASPAR_2022/REMAP_2020_Athaliana/results/GSE60141_WRKY30_Col-0-leaves-tnt-colamp/motifs/jaspar/pfm/GSE60141_WRKY30_Col-0-leaves-tnt-colamp ;



#/storage/scratch/JASPAR_2022/REMAP_2020_Athaliana/results/GSE60141_WRKY71_Col-0-leaves-tnt-col/motifs/jaspar/pfm/GSE60141_WRKY71_Col-0-leaves-tnt-col_peak-motifs_m3.tab         -bgfile /storage/scratch/JASPAR_2022/REMAP_2020_Athaliana/results/GSE60141_WRKY71_Col-0-leaves-tnt-col/peak-motifs/results/composition/GSE60141_WRKY71_Col-0-leaves-tnt-col_test_inclusive-1str-ovlp_2nt.txt         -t 5         -return sites > /storage/scratch/JASPAR_2022/REMAP_2020_Athaliana/results/GSE60141_WRKY71_Col-0-leaves-tnt-col/matrix_sites/GSE60141_WRKY71_Col-0-leaves-tnt-col_peak-motifs_m3.tf.sites


rule convert_RSAT_matrix_sites_to_BED:
    """
    Get the genomic coordinates (BED file) of the matrix sites.
    This rule is executed for each discovered motif.
    """
    input:
        os.path.join(config["out_dir"], "{TF}", "matrix_sites", "{TF}_peak-motifs_m{n}.tf.sites")
    output:
        os.path.join(config["out_dir"], "{TF}", "matrix_sites", "{TF}_peak-motifs_m{n}.tf.sites.bed")
    message:
        "; Obtaining genomic coordinates for sites - TF : {wildcards.TF} - Matrix number: {wildcards.n}"
    params:
        scripts_bin = config["bin"],
	genome_name = config['genome_name']
    priority:
        91
    shell:
        """
        awk -v species={params.genome_name} -f {params.scripts_bin}/sites-to-bed.awk {input} > {output}
        """


rule get_RSAT_matrix_sites_fasta:
    """
    Get the fasta sequences from genomic coordinates (BED file) of the matrix sites.
    This rule is executed for each discovered motif.
    """
    input:
        os.path.join(config["out_dir"], "{TF}", "matrix_sites", "{TF}_peak-motifs_m{n}.tf.sites.bed")
    output:
        os.path.join(config["out_dir"], "{TF}", "matrix_sites", "{TF}_peak-motifs_m{n}.tf.sites.fasta")
    message:
        "; Obtaining fasta sequences from genomic coordinates for TFBS sites - TF : {wildcards.TF} - Matrix number: {wildcards.n}"
    params:
        genome_fasta = config["genome_fasta"]
    priority:
        90
    shell:
        """
        bedtools getfasta -name -s -fi {params.genome_fasta} -bed {input} -fo {output}
        """


######################################################
## Motif centrality section:                        ##
## scan extended peaks + centrality p-value + plots ##
######################################################
rule Scan_JASPAR_PWM:
    """
    Scan the PWM in the extended peaks (501bp).
    This rule is executed for each discovered motif.
    """
    input:
        logos = os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "logos", "{TF}_peak-motifs_m{n}_logo.png"), \
        pwm = os.path.join(config["out_dir"], "{TF}", "motifs", "jaspar", "pwm", "{TF}_peak-motifs_m{n}.jaspar.pssm"), \
        peaks = os.path.join(config["out_dir"], "{TF}", "fasta", "{TF}.501bp.fa"), \
	sites = os.path.join(config["out_dir"], "{TF}", "matrix_sites", "{TF}_peak-motifs_m{n}.tf.sites.fasta")
    output:
        os.path.join(config["out_dir"], "{TF}", "scan", "501bp", "{TF}_m{n}.501bp.fa")
    message:
        "; Scanning PSSM - TF : {wildcards.TF} - Matrix number: {wildcards.n}"
    params:
        scripts_bin = config["bin"]
    priority:
        89
    shell:
        """
        {params.scripts_bin}/pwm_searchPFF {input.pwm} {input.peaks} 0.85 -b > {output}
        """


rule Calculate_centrimo_pvalue:
    """
    Calculate p-value for central enrichment.
    This rule is executed for each discovered motif.
    """
    input:
        sites = os.path.join(config["out_dir"], "{TF}", "scan", "501bp", "{TF}_m{n}.501bp.fa"),
        peaks = os.path.join(config["out_dir"], "{TF}", "extended_peaks", "{TF}.501bp.bed")
    output:
        os.path.join(config["out_dir"], "{TF}", "central_enrichment", "{TF}_m{n}.501bp.fa.sites.centrimo")
    message:
        "; Calculating central enrichment around peak summits - TF : {wildcards.TF} - Matrix number: {wildcards.n}"
    params:
        scripts_bin = config["bin"],
        centrimo_folder = os.path.join(config["out_dir"], "{TF}", "central_enrichment")
    priority:
        88
    shell:
        """
        mkdir -m 077 -p {params.centrimo_folder} ;
        nb_TFBS="$(wc -l {input.sites} | cut -d ' ' -f 1)" ;
        nb_peaks="$(wc -l {input.peaks} | cut -d " " -f 1)" ;
        {params.scripts_bin}/centrimo_pval {input.sites} ${{nb_TFBS}} ${{nb_peaks}} 250 > {output}
        """


rule generate_centrimo_plots:
    """
    Generate centrimo plots (local enrichment) around the peak summits.
    This rule is executed for each discovered motif.
    """
    input:
        fa = os.path.join(config["out_dir"], "{TF}", "scan", "501bp", "{TF}_m{n}.501bp.fa"), \
        fb = os.path.join(config["out_dir"], "{TF}", "central_enrichment", "{TF}_m{n}.501bp.fa.sites.centrimo")
    output:
        os.path.join(config["out_dir"], "{TF}", "central_enrichment", "{TF}_m{n}.501bp.fa.sites.centrimo.pdf")
    message:
        "; Scanning PSSM - TF : {wildcards.TF} - Matrix number: {wildcards.n}"
    params:
        scripts_bin = config["bin"]
    priority:
        87
    shell:
        """
        R --vanilla --slave --silent -f {params.scripts_bin}/centrimo_plot.R --args {input.fa} {output}
        """


##############################################################     
## Function that defines new wildcard "n" (motif number),   ##
## so the pipeline uses the checkpoint to run uninterupted. ##
##############################################################
def aggregate_motif_ind(wildcards):

    ## Name of the checkpoint rule
    checkpoint_output = checkpoints.RSAT_PSSM_to_JASPAR_format.get(**wildcards).output[0]

    ##
    ## TF = wildcard taken from the rule all
    ## n  = Motif number. Read the global wildcards from the checkpoint output directory
    file_names = expand(os.path.join(config["out_dir"], "{TF}", "central_enrichment", "{TF}_m{n}.501bp.fa.sites.centrimo.pdf"),
                        TF = wildcards.TF,
                        n = glob_wildcards(os.path.join(checkpoint_output, "{TF}_peak-motifs_m{n}.jaspar")).n)
    return file_names



########################################################################################
## For each dataset, select and annotate the motif with the lowest centrality p-value ##
########################################################################################

rule JASPAR_annotation_table:
    """
    Connect to JASPAR to download the annotation table of all the profiles in the current version.
    This information is used to annotate the motifs in a semi-automatic way
    """
    input:
        config["TF_Experiment_map"]
    output:
        JASPAR_ANN_TAB
    message:
        "; Generate JASPAR 2020 motif annotation table "
    params:
        scripts_bin = config["bin"],
        taxon = config["taxon"],
        out_dir = config["out_dir"]
    priority:
        86
    shell:
       """
       Rscript {params.scripts_bin}/Retrieve_matrix_information_from_JASPAR2020.R -o {params.out_dir} -t {params.taxon}
       """

       
rule choose_best_centrimo_experiment:
    """
    Select the motif with the best (lowest) centrimo p-value.
    This rule is executed for each dataset (experiment).
    """
    input:
        aggregate_motif_ind
    output:
        os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best")
    message:
        "; Selecting the most centrally enriched motif - TF : {wildcards.TF} "
    params:
        scripts_bin  = config["bin"],
        centrimo_dir = os.path.join(config["out_dir"], "{TF}", "central_enrichment"),
        nbmotifs     = config['top_motifs']
    priority:
        85
    shell:
        """
        bash {params.scripts_bin}/best_centrimo.sh -i {params.centrimo_dir} -m {params.nbmotifs} > {output}
        """


rule annotate_best_centrimo_experiment:
    """
    Assign the TF name to the selected motif.
    This rule is executed for the best experiment of each dataset (experiment). Make sure the experiment map has the experiment ID in the first field, and the TF name in the third field
    """
    input:
        tf_jaspar_map = config["TF_Experiment_map"],
        best_exp = os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best")
    output:
        os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best.TF_associated")
    message:
        "; Assigning TF name to the experiment: {wildcards.TF} "
    params:
        scripts_bin = config["bin"]
    priority:
        84
    shell:
        """
        bash {params.scripts_bin}/annotate_best_centrimo_experiment.sh {input.best_exp} {input.tf_jaspar_map} {output}
        """
# #perl {params.scripts_bin}/annotate_best_centrimo_experiment.pl --best {input.best_exp} --map {input.tf_jaspar_map} --output {output}



##########################
## Prepare curation PDF ##
##########################
rule best_centrimo_experiment_logo:
    """
    Export a PDF file with the motif logo, experiment, TF name and centrality p-value.
    This rule is executed for the best experiment of each dataset (experiment).
    """
    input:
        os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best.TF_associated")
    output:
        os.path.join(config["out_dir"], "{TF}", "central_enrichment", "selected_motif", "{TF}.501bp.fa.sites.centrimo.best.TF_associated.tex")
    message:
        "; Generating latex logo for selected motif: {wildcards.TF} "
    params:
        scripts_bin = config["bin"],
	logpval     = config['central_pvalue']
    priority:
        83
    shell:
        """
        bash {params.scripts_bin}/create_latex_logos.sh -l {params.scripts_bin}/latex_header.txt -i {input} -o {output} -t {params.logpval}
        """


rule Concat_annotated_experiments:
    """
    Concatenate the tables with the annotated experiments and the centrality p-value
    """
    input:
        MOST_ENRICHED_MOTIF_ASSOC_LOGO
    output:
        os.path.join(config["curation_dir"], "Annotated_experiments_cat.tab")
    message:
        "; Concatenating the tables with the annotated experiments and the centrality p-value "
    params: out_dir = config["out_dir"]
    priority:
        82
    shell:
        """
        ls {params.out_dir}/*/central_enrichment/selected_motif/*.501bp.fa.sites.centrimo.best.TF_associated | xargs cat > {output}
        """


rule Select_motifs_to_curate:
    """
    Select those motifs satisfying the -log10(centrality p-value) threshold.
    """
    input:
        os.path.join(config["curation_dir"], "Annotated_experiments_cat.tab")
    output:
        os.path.join(config["curation_dir"], "Selected_motifs_to_curate_log10_pval_{logpval}.tab")
    message:
        "; Selecting motifs to curate avobe threshold -log10(Centrality p-value): {{config['central_pvalue']}} "
    params:
        central_pval = config["central_pvalue"]
    priority:
        81
    shell:
        """
        cat {input} | awk -F"\t" '{{ if ($17 >= {params.central_pval}) {{ print }} }}' | uniq > {output}
        """


rule rename_jaspar_motif_header:
    """
    Write the correct name to the jaspar motifs.
    The motifs resulting from peak-motifs have a header that looks like this: 
    >peak-motifs_m3 peak-motifs_m3
    Change this header for an informative one
    >CTCF CTCF 
    """
    input:
        logos     = MOST_ENRICHED_MOTIF_ASSOC_LOGO, 
        exp_table = os.path.join(config["curation_dir"], "Selected_motifs_to_curate_log10_pval_{logpval}.tab")
    output:
        os.path.join(config["curation_dir"], "Renamed_log_pval_{logpval}.txt")
    message:
        "; Renaming jaspar motif header"
    priority:
        80
    shell:
        """
        perl -lne '@sl = split(/\\t/, $_); \
           $convert_cmd = " convert-matrix -i ".$sl[0]. " -from jaspar -to jaspar -attr id ".$sl[3]." -attr name ".$sl[3]." -return counts > ".$sl[0]."tmp" ; \
           $rename_file = "mv ".$sl[0]."tmp ".$sl[0] ; \
           system($convert_cmd); \
           system($rename_file); ' {input.exp_table} ;
        echo "Renamed motifs" > {output}
        """


rule Motifs_to_curate_PDF:
    """
    Concat the PDF of the selected motifs
    """
    input:
        os.path.join(config["curation_dir"], "Selected_motifs_to_curate_log10_pval_{logpval}.tab")
    output:
        os.path.join(config["curation_dir"], "Selected_motifs_to_curate_log10_pval_{logpval}.pdf")
    message:
        "; Concatenating PDFs with the motifs to curate "
    params: central_pval = config["central_pvalue"]
    priority:
        79
    shell:
        """
        PDF_FILES=` awk -F"\t" '{{ if($17 >= {params.central_pval}){{ print $20 }}}}' {input} | uniq | xargs `
        pdfunite $PDF_FILES {output}
        """
