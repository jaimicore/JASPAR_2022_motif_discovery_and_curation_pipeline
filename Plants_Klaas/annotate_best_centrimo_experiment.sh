best=$1;
map=$2;
output=$3;

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

  jaspar_file=$(find $dset/motifs/jaspar/pfm/ -name ${file_bname_motif_parts[0]}"_"${file_bname_motif_parts[1]}"_"${file_bname_motif_parts[2]}"_peak-motifs_"${file_bname_motif_parts[3]}".jaspar")

  TF=${file_bname_motif_parts[1]}

  data_source="ChIP-seq"

  sites_bed_file=$(echo ${dset}/matrix_sites/${file_bname_motif_parts[0]}"_"${file_bname_motif_parts[1]}"_"${file_bname_motif_parts[2]}"_peak-motifs_"${file_bname_motif_parts[3]}".tf.sites.bed")

  sites_fa_file=${sites_bed_file//".bed"/".fa"}

  centrimo_pdf_file=${file}.pdf

  logo_file=$(find $dset/motifs/jaspar/logos/ -name ${file_bname_motif_parts[0]}"_"${file_bname_motif_parts[1]}"_"${file_bname_motif_parts[2]}"_peak-motifs_"${file_bname_motif_parts[3]}"_logo.png")

  dset_id=$(basename $dset)

  pdf_selected=${dset}/central_enrichment/selected_motif/${file_bname_motif_parts[0]}_${file_bname_motif_parts[1]}_${file_bname_motif_parts[2]}.501bp.fa.sites.centrimo.best.TF_associated.pdf

  png_file=$(echo ${dset}/central_enrichment/${file_bname_motif_parts[0]}"_"${file_bname_motif_parts[1]}"_"${file_bname_motif_parts[2]}"_"${file_bname_motif_parts[3]}"_501bp_fa_sites_centrimo.png")

  printf "${jaspar_file}\t${TF}\t${data_source}\t${sites_bed_file}\t${sites_fa_file}\t${pval}\t${centrimo_pdf_file}\t${logo_file}\t${dset_id}\t${pdf_selected}\t${png_file}\n" >> $output


done < $best
