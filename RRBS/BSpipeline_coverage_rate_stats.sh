# collect metrics
echo "Reading bed files."
for i in */*.bed; do
    n=`basename $i | sed 's/_03_mcall_.*.CpG.bed//'`;
    grep -v totalC $i | cut -f4,5,7 | awk -v n=$n '{print n"\t"$0}';
done | gzip -c - >stats.txt.gz

# plot metrics
echo "Plotting statistics."
R --slave --vanilla <<EOF
    require(ggplot2)
    stats <- read.table('stats.txt.gz', header=F, col.names=c('sample','rate','coverage','strand'))
    stats2 <- subset(stats, coverage >= 10 & coverage < 150)
        
    pdf(file='stats.pdf', width=10)
	ggplot(stats, aes(x=rate, color=sample)) + geom_line(stat='density', size=2) + theme_classic(16) + ggtitle('Methylation rate distribution (coverage > 0)') + xlab('Methylation rate') + ylab('Density')
	ggplot(stats, aes(x=coverage, fill=sample)) + geom_histogram(binwidth=2, color='black') + theme_classic(16) + ggtitle('Coverage distribution [0, 150] ') + xlab('Coverage') + ylab('Density') + coord_cartesian(xlim=c(0,150)) + facet_wrap(~sample, ncol=2)
	ggplot(stats, aes(x=sample, fill=strand)) + geom_bar(width=0.7, color='black') + theme_classic(16) + ggtitle('Strand coverage') + xlab('Sample') + ylab('Frequency') + scale_fill_manual(values=c('lightgrey','grey','black'))
	
	ggplot(stats2, aes(x=rate, color=sample)) + geom_line(stat='density', size=2) + theme_classic(16) + ggtitle('Methylation rate distribution (coverage [10, 150])') + xlab('Methylation rate') + ylab('Density')
	ggplot(stats2, aes(x=coverage, fill=sample)) + geom_histogram(binwidth=2, color='black') + theme_classic(16) + ggtitle('Coverage distribution (coverage [10, 150])') + xlab('Coverage') + ylab('Density') + facet_wrap(~sample, ncol=2)
	ggplot(stats2, aes(x=sample, fill=strand)) + geom_bar(width=0.7, color='black') + theme_classic(16) + ggtitle('Strand coverage (coverage [10, 150])') + xlab('Sample') + ylab('Frequency') + scale_fill_manual(values=c('lightgrey','grey','black'))
	dev.off()
EOF

# clean up
rm stats.txt.gz

echo "Done."
