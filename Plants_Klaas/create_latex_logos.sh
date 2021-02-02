#!/usr/bin/env bash

################################################################################
#
################################################################################

usage="""
\n
$0 -l <latex header file> -i <input file> -o <output file>
\n
""";

\usepackage[export]{adjustbox}
\usepackage{graphicx}

TEMP=`getopt -o ho:i:l: -- "$@"`;

if [ $? != 0 ]
then
  echo "Terminating..." >&2;
  exit 1;
fi
eval set -- "$TEMP";

input="NA";
output="NA";
latexheader="NA";
while true
do
  case "$1" in
    -i) input=$2; shift 2;;
    -o) output=$2; shift 2;;
    -l) latexheader=$2; shift 2;;
    -h) echo -e $usage; exit;;
    --) shift; break;;
    *) echo "Internal error!"; exit 1;;
  esac
done

# Checking options and parameters

if [ ! -e $input ]
then
    echo -e "\nPlease provide the input file\n";
    exit 1;
fi;

if [ ! -e $latexheader ]
then
    echo -e "\nPlease provide the latex header file\n";
    exit 1;
fi;

cat $latexheader > $output;

cat $input | \
while IFS=$'\t' read -r -a line
do

    variant=$(basename ${line[0]})
    variant=${variant%%".jaspar"}
    variant=(${variant//"_"/" "})
    tfname=${line[3]}
    #tfname=$(echo $line | cut -d ' ' -f 2 | tr '_' '-');
    exp_ID=$(echo ${line[9]}"-"${variant[4]} | tr '_' '-')
    #exp_ID=$(echo $line | cut -d ' ' -f 9 | tr '_' '-');
    centrimo_pval=${line[16]}
    #centrimo_pval=$(echo $line | cut -d ' ' -f 6);

    ## Read motif logo
    motif_logo_original=${line[18]} # change
    #motif_logo_original=$(echo $line | cut -d ' ' -f 8);
    motif_logo_original_dirname=$(dirname $motif_logo_original);
    # prevent error: ! LaTeX Error: Unknown graphics extension: .16_LE_WA_peak-motifs_m5_logo.png.
    # This removes any internal dots in the basename of the file, except for the dot in the extension (e.g. except the .png part)
    motif_logo_new=$(basename `echo $motif_logo_original` | perl -pe 's/(\.)(?=[^.]*\.)/\_/');
    motif_logo=${motif_logo_original_dirname}/${motif_logo_new}

    ## Read centrality plot
    centrimo_plot=${line[17]} # change
    #centrimo_plot_original=$(echo $line | cut -d ' ' -f 11);
    #centrimo_plot_original_dirname=$(dirname $centrimo_plot_original);
    # prevent error: ! LaTeX Error: Unknown graphics extension: .16_LE_WA_peak-motifs_m5_logo.png.
    # This removes any internal dots in the basename of the file, except for the dot in the extension (e.g. except the .png part)
    #centrimo_plot_new=$(basename `echo $centrimo_plot_original` | perl -pe 's/(\.)(?=[^.]*\.)/\_/');
    #centrimo_plot=${centrimo_plot_original_dirname}/${centrimo_plot_new}

    # this create a symbolic link so that that new motif filename will actually exist
    #ln -s $motif_logo_original $motif_logo
    #motif_pdf=$(echo $line | cut -d ' ' -f 7);


   # if (( $(echo "$centrimo_pval < 0." | bc -l) ))
    #then

    #printf "\\section*{$tfname}\n\\section*{$exp_ID}\nCentrality p-value = $centrimo_pval \\\\ \n\\includegraphics{$motif_logo}\nCentrimo plot \\\\ \n\\includegraphics[width=8cm, height=8cm]{$centrimo_plot}\n" >> $output
    echo "\\section*{$tfname}" >> $output;
    echo "\\section*{$exp_ID}" >> $output;
    echo "-log10(Centrality p-value) = $centrimo_pval \\\\" >> $output;
    echo "\\includegraphics{$motif_logo}" >> $output;
    echo "Centrimo plot \\\\" >> $output;
    echo "\\includegraphics[width=8cm, height=8cm]{$centrimo_plot}" >> $output;
    #fi;
done;
echo "\\end{document}" >> $output;
outdir=$(dirname $output);

pdflatex -output-directory $outdir $output;
