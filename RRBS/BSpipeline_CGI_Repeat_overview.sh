# create new folders
echo "Create folder if not existend."
mkdir -p overview

# annotation files
CGI="/project/bioinf_meissner/references/${1}/annotation/CGIs.bed"
Repeats="/project/bioinf_meissner/references/${1}/annotation/Repeats_Class.bed"

echo -ne "#chr\tstart\tend\titem\trate\tcount\tsample\n" >overview/overview.txt


# average methylation rate per CpG island
echo "Calculate average methylation rate at CpG islands."
for i in bedgraph/*.bedgraph; do 
    n=`basename $i | sed 's/.bedgraph//'`; 
    echo ${n}
    /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools intersect -wa -wb -a ${CGI} -b ${i} | sed 's/CpG:_.*\t/CpG_Island\t/' | /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools groupby -g 1,2,3,4 -c 5 -o mean,count | perl -ane 'if($F[5]>=3){print $_}' | awk -v n=${n} '{print $0"\t"n}' >>overview/overview.txt
done

# average methylation rate per CpG island
echo "Calculate average methylation rate at SINE, LINE, LTR."
for i in bedgraph/*.bedgraph; do 
    n=`basename $i | sed 's/.bedgraph//'`; 
    echo ${n}
    /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools intersect -wa -wb -a ${Repeats} -b ${i} | grep "LINE\|SINE\|LTR" | grep -v "?" | /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools groupby -g 1,2,3,4 -c 8 -o mean,count | perl -ane 'if($F[5]>=3){print $_}' | awk -v n=${n} '{print $0"\t"n}' >>overview/overview.txt
done


# plot metrics
echo "Plotting statistics."
R --slave --vanilla <<EOF
    require(ggplot2)
    data <- read.table('overview/overview.txt', header=T, comment.char='')

    pdf(file='overview/overview', width=10)
	ggplot(data, aes(x=sample, fill=sample, y=rate)) + geom_boxplot(width=0.7, outlier.size=0.25) + theme_classic(16) + ggtitle('Average methylation') + ylab('Methylation rate') + ylab('') + facet_wrap(~item, nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	dev.off()
EOF

echo "Done."


