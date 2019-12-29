# create new folders
echo "Create bigwig and bedgraph folder if not existend."
mkdir -p ../bigwig
mkdir -p ../bedgraph

# filter for coverage
echo "Filter for coverage [10,150] and create bedgraph and bigwig files."
while read -r ID NAME
do 
    echo ${NAME}
    echo -ne "#chr\tstart\tend" | awk -v n=${NAME} '{print $0"\t"n}' >../bedgraph/${NAME}.bedgraph
    grep -v ratio RRBS_*_${ID}/${ID}*_mcall_rrbs.CpG.bed | perl -ane 'if($F[4]>=10 && $F[4]<=150){print $_}' | grep -v "chrM" | cut -f1,2,3,4 | /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools sort -i stdin >>../bedgraph/${NAME}.bedgraph
    /project/bioinf_meissner/src/UCSCtools/bedGraphToBigWig ../bedgraph/${NAME}.bedgraph /project/bioinf_meissner/references/${1}/${1}.chrom.sizes ../bigwig/${NAME}.bw
done <../ID2NAME.txt

