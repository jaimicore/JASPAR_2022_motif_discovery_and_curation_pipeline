best=$1;
map=$2;
output=$3;

#printf "PWM\tcurrent_BASE_ID\tcurrent_VERSION\tTF_NAME\tUniprot\tTAX_ID\tclass\tfamily\tTFBSshapeID\tData\tSource\tValidation\tComment\tAddition_or_Upgrade_or_Non-validated_(A_or_U_or_N_)\tBED\tFASTA\tCentrality_pval\tCentrality_plot\tLogo_file\n" > $output

while IFS=$'\t' read -r -a line
do

  file=${line[0]}
  pval=${line[1]}
  pval=$(awk -v a=$pval 'BEGIN{print (-1 * a)}')

  file_bname=$(basename $file)
  file_bname_motif=${file_bname%%".501bp.fa.sites.centrimo"}
  file_bname_motif_parts=(${file_bname_motif//"_"/" "})

  dset=$(dirname $file)
  dset=$(dirname $dset)

  jaspar_file=$(find $dset/motifs/jaspar/pfm/ -name ${file_bname_motif_parts[0]}"_peak-motifs_"${file_bname_motif_parts[1]}".jaspar")

  current_base_id="NULL"
  current_version="NULL"

  TF=${file_bname_motif_parts[0]}

  uniprot="NULL"
  TAX_ID="NULL"
  class="NULL"
  family="NULL"
  TFBSshapeID="NULL"

  source="NULL"
  validation="NULL"
  comment="NULL"
  add_upgrade_nonval="NULL"

  sites_bed_file=$(echo ${dset}/matrix_sites/${file_bname_motif_parts[0]}"_peak-motifs_"${file_bname_motif_parts[1]}".tf.sites.bed")

  sites_fa_file=${sites_bed_file//".bed"/".fasta"}

  #centrimo_pdf_file=$(find $dset/central_enrichment -name ${file_bname_motif_parts[0]}"_"${file_bname_motif_parts[1]}"_"${file_bname_motif_parts[2]}"_"${file_bname_motif_parts[3]}".501bp.fa.sites.centrimo.pdf")

  logo_file=$(find $dset/motifs/jaspar/logos/ -name ${file_bname_motif_parts[0]}"_peak-motifs_"${file_bname_motif_parts[1]}"_logo.png")

  dset_id=$(basename $dset)

  pdf_selected=${dset}/central_enrichment/selected_motif/${file_bname_motif_parts[0]}.501bp.fa.sites.centrimo.best.TF_associated.pdf

  centrimo_png_file=$(echo ${dset}/central_enrichment/${file_bname_motif_parts[0]}"_"${file_bname_motif_parts[1]}"_501bp_fa_sites_centrimo.png")

  printf "${jaspar_file}\t${current_base_id}\t${current_version}\t${TF}\t${uniprot}\t${TAX_ID}\t${class}\t${family}\t${TFBSshapeID}\t${dset_id}\t${source}\t${validation}\t${comment}\t${add_upgrade_nonval}\t${sites_bed_file}\t${sites_fa_file}\t${pval}\t${centrimo_png_file}\t${logo_file}\t${pdf_selected}\n" >> $output

done < $best
