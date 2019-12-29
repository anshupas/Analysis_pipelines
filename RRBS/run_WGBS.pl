#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Cwd qw(abs_path cwd chdir); # chdir: overrides internal chdir() and updates $ENV{PWD}!

## Needs implementation:
##  * RRBS 
##    - quality / adaptor trimming
##    - diversity adaptor trimming (Nugene)
##    - mapping using bismark with --rrbs *CHECK*
##    - no deduplication
## 
##  * PBAT
##    - quality / adaptor trimming
##    - random bases trimming PE/SE, see http://www.cell.com/cell/pdfExtended/S0092-8674(15)00564-4
##    - mapping using bismark --pbat 
##    - deduplication using bismark
## 
##  * WGBS
##    - single end data?

## A word about "Direction":
## 
## A directional protocol yields reads from both strands of the original sample-DNA. The first read
## of a pair (or every read for single-end sequencing) is known to be sequenced either from the
## original-top (OT) or the original-bottom (OB) strand. The second read of a pair is sequenced
## from a complementary strand, either ctOT or ctOB. 
## At the time of writing, examples of directional protocols include:
## 
##  • Illumina TruSeq DNA Methylation Kit (formerly EpiGnome)
##  • Kits from the NuGen Ovation family of products
##  • Swift Accel-NGS Methyl-seq DNA Library Kit
##  • Libraries prepared by the 'Lister' method
##
## In a non-directional protocol, the first read of a pair may come from any of the four strands: OT,
## OB, ctOT or ctOB. Examples include:
##
##  • Zymo Pico Methyl-Seq Library Kit
##  • Bioo Scientific (Perkin Elmer) NEXTflex Bisulfite-Seq Kit
##  • Libraries prepared by the 'Cokus' method
## 
## Source: 
## http://resources.qiagenbioinformatics.com/manuals/bisulfite-sequencing/current/BisulfiteSequencing_Plugin_User_Manual.pdf

BEGIN {
    $ENV{PATH} .= ":/package/sequencer/bin";
}

## global stuff, mandatory
my ($fq1,$fq2);
my $sample;
my $genome;

## presets
my $c_threads      = 36; # computing threads
my $p_threads      =  8; # processing threads
my $no_dedup       =  0; 
my $is_swift       =  0; # WGBS Swift library
my $is_pbat        =  0; # Post-bisulfite adaptor tagging (PBAT) library
my $is_nondir      =  0; # if 1, library is non-directional
my $maxins         = 1000; # max insert size for PE mapping
my $trim_off       = 10; # trim fixed number of bases, e.g. SWIFT 
my $is_gz          =  1; # compressed input?
my $adaptor        = 'illumina'; 
my $workflow       = 'wgbs';
my $dateflag       = strftime("%Y%m%d", localtime);
my $cwd            = cwd();

my %progs = (
    'base'      => '/package/sequencer',
    'bamtools'  => '/package/sequencer/bin/bam', # overlap clipping -- necessary?
    'fastqc'    => '/package/sequencer/bin/fastqc',
    'cutadapt'  => '/package/sequencer/anaconda3/envs/ngs/bin/cutadapt',
    'samtools'  => '/package/sequencer/bin/samtools',

    'bismark'   => '/package/sequencer/bin/bismark',
    'dedup'     => '/package/sequencer/bin/deduplicate_bismark',
    'methylex'  => '/package/sequencer/bin/bismark_methylation_extractor',
    'report'    => '/package/sequencer/bin/bismark2report',
);

my %genomes = (
    'hg19'  => '/project/genomes/BiSeq/bismark_hg19',
    'hg38'  => '/project/genomes/BiSeq/bismark_hg38',
    'mm9'   => '/project/genomes/BiSeq/bismark_mm9',
    'mm10'  => '/project/genomes/BiSeq/bismark_mm10',
    'jc01'  => '/project/genomes/BiSeq/bismark_musD_insert_seq_pmSEQ' # "alignment to retroelement sequence", 6 Feb 2018 12:41:34
);

##   Trimming (cutadapt) adaptor sequences:
my %adaptors = (
    'illumina'  => 'AGATCGGAAGAGC',
    'small_rna' => 'TGGAATTCTCGG',
    'nextera'   => 'CTGTCTCTTATA',
    'file'      => '/package/sequencer/screenlib/cutadapt.tab'
);

my %workflows  =(
    'w'     => 'WGBS',
    'r'     => 'RRBS',
    'wgbs'  => 'WGBS',
    'rrbs'  => 'RRBS'
);

## This should be considered, at least approximately.
## There will be some additional compression threads.
## Currently this pipeline should be run with C_PROC=36
## and for mxqsub THREADS=36 -- FIXME sklages 2018-02-14
my %bismark_thread_sets = (
    ## [4,2] = 4 bismark instances in parallel running 2 bowtie jobs (C->T/G->A) each with 2 threads
    '1'  => [4,2], # uses 4x2x2 = 16 threads
    '2'  => [4,4], # uses 4x4x2 = 32 threads
    '3'  => [8,2], # uses 8x2x2 = 32 threads - more I/O
    '4'  => [8,4]  # uses 8x4x2 = 64 threads - more I/O
);

my $USAGE = <<"USAGE";

Usage : $0 --c_threads INT|36 usw und sofort.

USAGE

GetOptions (
   'fq1=s'       => \$fq1,
   'fq2=s'       => \$fq2,
   'sample=s'    => \$sample,
   'genome=s'    => \$genome,
   'c_threads=i' => \$c_threads,
   'p_threads=i' => \$p_threads,
   'workflow=s'  => \$workflow,
   'adaptor=s'   => \$adaptor,
   'no_dedup'    => \$no_dedup,
   'swift'       => \$is_swift,
   'pbat'        => \$is_pbat,
   'nondir'      => \$is_nondir,
   'maxins'      => \$maxins

) or die "No! No. You don't understand.\nThere was something in the woods, David... and I think it's in here with us... now.\n";


## ====================================================================
## Parameter Checks
## ====================================================================

## fastq R1
if (! defined $fq1) {
    warn "[ERROR] Please supply Fastq file (R1).\n";
    die $USAGE;
} elsif (! -r $fq1 ) {
    die "[ERROR] Cannot read fastq file '$fq1': $!\n";
}

## fastq R2
if (! defined $fq2) {
    warn "[ERROR] Please supply Fastq file (R2).\n";
    die $USAGE;
} elsif (! -r $fq2 ) {
    die "[ERROR] Cannot read fastq file '$fq2': $!\n";
}

## sample prefix
if (! defined $sample) {
    warn "[ERROR] Please supply sample prefix.\n";
    die $USAGE;
} elsif ($sample =~ /\s/) {
    die "[ERROR] Don't use whitespace in sample prefix.\n";
}

## reference genome
if (! exists $genomes{$genome}) {
    warn "[ERROR] Provided reference name not supported. Use one of 'mm9','mm10','hg19' or 'hg38'.\n";
    die $USAGE;
}

## reference genome
if (! exists $adaptors{$adaptor}) {
    warn "[ERROR] Provided adpter name not supported. Use one of 'illumina','nextera' or 'small_rna'.\n";
    die $USAGE;
}

## workflow, currently WGBS and RRBS
if (! exists $workflows{$workflow}) {
    warn "[ERROR] Provided workflow not supported. Use one of 'w'(WGBS) or 'r' (RRBS).\n";
    die $USAGE;
}

## dedupping?
if ($no_dedup) {
    warn "[WARNING] Be careful not to dedup amplicon libraries.\n";
}

## really zipped?
unless ($fq1=~/\.gz$/) {
    $is_gz=0;
}

## swift/pbat are exclusive
if($is_swift && $is_pbat) {
    warn "[ERROR] You cannot use 'swift' and 'pbat' at the same time.\n";
    die $USAGE;
}

## usually 10 bases with Swift libraries
unless ($trim_off=~/\d+/ && $trim_off>0) {
    warn "[ERROR] Number of nucleotides to trim off must be positive integer, not '$trim_off'.\n";
    die $USAGE;
}

## max insert size
if($maxins<1 || $maxins=~/[^\d]/) {
    warn "[ERROR] Maximum insert size must be positive integer, not '$maxins'.\n";
    die $USAGE;
}

## ====================================================================
## Prepare Output Directory
## ====================================================================

## absolute paths to fastq input data
my $fq1_abspath = abs_path($fq1);
my $fq2_abspath = abs_path($fq2);

## WGBS_20180215_AM-WGBS-001
my $ANALYSIS_DIR = $cwd . '/' . $workflows{$workflow} . "_" . $dateflag . "_" . $sample;

mkdir($ANALYSIS_DIR,0770)
    or die "[ERROR] Could not create '$ANALYSIS_DIR': $!\n";

chdir($ANALYSIS_DIR)
    or die "[ERROR] Could not change to '$ANALYSIS_DIR': $!\n";

## be sure if we have to add gz extension to "new" filename
my $fq1_input = ($is_gz) ? "${sample}_R1.fq.gz" : "${sample}_R1.fq";
my $fq2_input = ($is_gz) ? "${sample}_R2.fq.gz" : "${sample}_R2.fq";

symlink($fq1_abspath,$fq1_input)
    or die "[ERROR] Could not create symlink '$fq1_input' in '$ANALYSIS_DIR': $!\n";
symlink($fq2_abspath,$fq2_input)
    or die "[ERROR] Could not create symlink '$fq2_input' in '$ANALYSIS_DIR': $!\n";

## Logs / Metrics should be stored separately
my $METRICS_DIR = "$ANALYSIS_DIR/fastqc";
mkdir($METRICS_DIR,0770)
    or die "[ERROR] Could not create '$ANALYSIS_DIR': $!\n";

## command logging
my $LOG_FILE = "$ANALYSIS_DIR/${sample}_workflow.log";


## Initial QC of input files
run_fastqc($fq1_input,$METRICS_DIR);
run_fastqc($fq2_input,$METRICS_DIR);


## ====================================================================
## Adaptor / Quality Trimming
## ====================================================================
## Depending on the kind of library used, we will always need some kind of trimming
## https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html#viii-notes-about-different-library-types-and-commercial-kits

## [0] = fastq 1, [1] = fastq 2
my @trimmed_fastq_set;
if($is_swift) {
    @trimmed_fastq_set = trim_swift_adaptor_PE($fq1_input,$fq2_input,$adaptor,$sample,$trim_off);
}
elsif($is_pbat) {
    @trimmed_fastq_set = (); ## DUMMY
}
else {
    @trimmed_fastq_set = trim_adaptor_PE($fq1_input,$fq2_input,$adaptor,$sample);
}

## Trimming Step Quality Control
run_fastqc($trimmed_fastq_set[0],$METRICS_DIR);
run_fastqc($trimmed_fastq_set[1],$METRICS_DIR);


## ====================================================================
## Run bismark with trimmed fastq on $genome
## ====================================================================
my $basename = "${sample}_Bismark_${genome}";
my $bam_bismark = run_bismark_alignment(
                    $trimmed_fastq_set[0], 
                    $trimmed_fastq_set[1], 
                    $genomes{$genome}, 
                    $sample,
                    $basename
);


## ====================================================================
## Deduplication (WGBS-only)
## ====================================================================
my $bam = $bam_bismark;
my $bam_nsorted;
my $bam_dedup;

## do dedup
if(!$no_dedup) {

    ## sort BAM file by names for deduplication
    $bam_nsorted = sort_bam_by_name($bam_bismark);

    ## deduplicate n-sorted BAM
    $bam_dedup = deduplicate_alignment($bam_nsorted);

    $bam = $bam_dedup;
}


## ====================================================================
## Bismark methylation extractor
## ====================================================================
## The IDs of Read 1 (A00417:7:H72GTDMXX:2:1368:9335:7075_1:N:0:ACAGTG) and 
## Read 2 (A00417:7:H72GTDMXX:2:2370:5638:14606_1:N:0:ACAGTG) are not the same. 
## This might be the result of sorting the paired-end SAM/BAM files by chromosomal 
## position which is not compatible with correct methylation extraction. 
## Please use an unsorted file instead or sort the file using 'samtools sort -n' 
## (by read name).
extract_methylation_calls($bam,$genomes{$genome});


## ====================================================================
## Finally write report
## ====================================================================
my $report = create_bismark_report($basename);


print "\n#### PIPELINE DONE ####\n\n";

print <<'EOM';

 ## Reporting 
 see https://sequencing.qcfail.com/articles/library-end-repair-reaction-introduces-methylation-biases-in-paired-end-pe-bisulfite-seq-applications/
 
 Created HTML Report $report
 
 Plus a bedcov or bigBed file for visual inspection in IGV.

 Check 'Methylation information will now be written into a genome-wide cytosine report' issue. 
 Have all file been written?
 
EOM


exit;



### ===================================================================
# Name .....: 
# Task .....: 
# Parameter : 
# Output ...:
# Returns ..:
### ===================================================================
sub log_command {
    ## Logging / program calls
    my ($desc,$cmd) = @_;
    my $date_time = strftime("%Y-%m-%d %H:%M:%S", localtime);

    open(my $log, ">>", $LOG_FILE)
        or die "[ERROR] Problem logging command: $!\n";

    print "## $desc\n[$date_time] $cmd\n\n";
    print $log "## $desc\n[$date_time] $cmd\n\n";

    close($log);
}


### ===================================================================
# Name .....: 
# Task .....: 
# Parameter : 
# Output ...:
# Returns ..:
### ===================================================================
sub create_bismark_report {
    my $basename=shift;
    ## AM-WGBS-001_Bismark_mm10.bam.alignment.report.txt
    ## AM-WGBS-001_Bismark_mm10.nsorted.deduplicated.M-bias.txt
    ## AM-WGBS-001_Bismark_mm10.nsorted.deduplicated_splitting_report.txt
    ## AM-WGBS-001_Bismark_mm10.nsorted.deduplication_report.txt
    ## ^^^^^^^^^^^ sample
    ## ^^^^^^^^^^^^^^^^^^^^^^^^ basename

    my $report_html = "${basename}_Overview.html";

    ## the lazy way :-)
    my @rep_align = glob ("${basename}*alignment.report.txt");
    my @rep_mbias = glob ("${basename}*M-bias.txt");
    my @rep_split = glob ("${basename}*splitting_report.txt");
    my @rep_dedup = glob ("${basename}*deduplication_report.txt");
    my @rep_nucl  = glob ("${basename}*nucleotide_stats.txt");

    my $cmd = "$progs{report}";
    $cmd .= " --output $report_html";
    $cmd .= " --alignment_report $rep_align[0]";
    $cmd .= " --dedup_report $rep_dedup[0]";
    $cmd .= " --splitting_report $rep_split[0]";
    $cmd .= " --mbias_report $rep_mbias[0]";

    if(@rep_nucl>0) {
        $cmd .= " --nucleotide_report $rep_nucl[0]";
    }

    log_command("Creating Bismark Report", $cmd);

    system($cmd)
        and die "[ERROR] Something went wrong while running  '$cmd'.\n";

    return $report_html;
}

### ===================================================================
# Name .....: extract_methylation_calls()
# Task .....: bismark_methylation_extractor (extract methylation calls)
# Parameter : n-sorted BAM 
#             Please use an unsorted file instead or sort the file using 
#             'samtools sort -n' (by read name).
# Output ...: Set of files
# Returns ..: -
# URL ......: -
### ===================================================================
sub extract_methylation_calls {
    my ($bam,$genome_folder)=@_;
    my $time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);

    print STDERR "[$time_stamp] Running 'bismark_methylation_extractor' on $bam ...\n";

    my $cmd = "$progs{methylex}";
    $cmd .= " --gzip";
    $cmd .= " --bedGraph";
    $cmd .= " --merge_non_CpG";
    $cmd .= " --buffer_size 20G";
    $cmd .= " --cytosine_report";
    $cmd .= " --genome_folder $genome_folder";
    $cmd .= " --paired-end";
    $cmd .= " --samtools_path $progs{base}/bin";
    $cmd .= " --parallel $p_threads";
    $cmd .= " --remove_spaces";
    $cmd .= " $bam";

    log_command("Extract Methylation Calls", $cmd);

    system($cmd)
        and die "[ERROR] Something went wrong while running  '$cmd'.\n";
}


### ===================================================================
# Name .....: sort_bam_by_name()
# Task .....: Sort BAM by name (as needed by Bismark's dedup tool)
# Parameter : BAM file
# Output ...: n-sorted BAM file
# Returns ..: n-sorted BAM file name
# URL ......: -
### ===================================================================
sub sort_bam_by_name {
    my($us_bam)=@_;

    my $time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);

    (my $ns_bam = $us_bam) =~ s/(.+)\.bam$/${1}.nsorted.bam/;

    print STDERR "[$time_stamp] Running 'samtools sort -n' on $us_bam ...\n";

    my $cmd = "$progs{samtools} sort -n --threads $p_threads -m 4G -O BAM -o $ns_bam $us_bam";

    log_command("Sort BAM file by name", $cmd);

    system($cmd)
        and die "[ERROR] Something went wrong while sorting BAM '$us_bam'.\n";

    return $ns_bam;
}


### ===================================================================
# Name .....: sort_bam_by_position()
# Task .....: Sort BAM by position
# Parameter : BAM file
# Output ...: p-sorted BAM file
# Returns ..: p-sorted BAM file name
# URL ......: -
### ===================================================================
sub sort_bam_by_position {
    my($us_bam)=@_;

    my $time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);

    (my $ns_bam = $us_bam) =~ s/(.+)\.bam$/${1}.psorted.bam/;

    print STDERR "[$time_stamp] Running 'samtools sort ' on $us_bam ...\n";

    my $cmd = "$progs{samtools} sort --threads $p_threads -m 4G -O BAM -o $ns_bam $us_bam";

    log_command("Sort BAM file by position", $cmd);

    system($cmd)
        and die "[ERROR] Something went wrong while sorting BAM '$us_bam'.\n";

    return $ns_bam;
}


### ===================================================================
# Name .....: deduplicate_alignment()
# Task .....: Deduplicate BAM file (WGBS only)
# Parameter : BAM file
# Output ...: Deduplicated BAM file
# Returns ..: Deduplicated BAM file name
# URL ......: -
### ===================================================================
sub deduplicate_alignment {
    my $bam = shift;
    my $time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);
    print STDERR "[$time_stamp] Bismark Dedupper on '$bam' ...\n";
    ## --> "${sample}_Bismark_${genome}.nsorted.deduplicated.bam
    ## AM-WGBS-001_Bismark_mm10.nsorted.deduplicated.bam
    ## AM-WGBS-001_Bismark_mm10.nsorted.deduplication_report.txt
    ##                                 ^^^^^^^^^^^^^^^^^^^^^^^^^
    my $cmd = "$progs{dedup} --paired --bam --samtools_path $progs{base}/bin $bam";

    log_command("Deduplicate alignment", $cmd);

    my $rc = system($cmd) 
        and die "[ERROR] Some problem running Bismark Deduplication on '$bam'.\n";

    ## new BAM filename
    $bam =~ s/(.+)\.bam$/${1}.deduplicated.bam/;

    unless(-e $bam) {
        print STDERR "\n[WARNING] BAM file '$bam' does not exist. Look for depuplicated version and methyl call extraction manually.\n\n";
        exit;
    }
    return $bam;
}


### ===================================================================
# Name .....: run_bismark_alignment()
# Task .....: Runs the actually alignment
# Parameter : Fastq file pair, genome folder, sample name
# Output ...: BAM file, some stats
# Returns ..: BAM file name
# URL ......: https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
### ===================================================================
sub run_bismark_alignment {
    my ($fq1,$fq2,$genome_folder,$sample,$basename)=@_;

    my $time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);
    my $seed_mismatches = 0; # default

    ## uses 4x4x2 = 32 threads
    my $b_multicore = 4;
    my $b_bowtie    = 4;

    if($c_threads<24) {
        $b_multicore = 4;
        $b_bowtie    = 2;
    }

    my $cmd = $progs{bismark};
    $cmd .= " --parallel $b_multicore";
    $cmd .= " -p $b_bowtie";
    $cmd .= " --path_to_bowtie $progs{base}/bin";
    $cmd .= " --samtools_path $progs{base}/bin";
    $cmd .= " --rg_id $sample";
    $cmd .= " --rg_sample $sample";
    $cmd .= " --output_dir $ANALYSIS_DIR";
    $cmd .= " --temp_dir $ANALYSIS_DIR/tmp.bismark";
    $cmd .= " --bowtie2";
    $cmd .= " -N $seed_mismatches";
    $cmd .= " --genome_folder $genome_folder";
    $cmd .= " --nucleotide_coverage";
    $cmd .= " --non_directional" if $is_nondir;
    $cmd .= " --maxins $maxins";
    $cmd .= " -1 $fq1";
    $cmd .= " -2 $fq2";

    print STDERR "[$time_stamp] Running Bismark $sample ...\n";

    log_command("Bismark alignment", $cmd);

    my $rc = system($cmd) 
        and die "[ERROR] Some problem running Bismark on '$sample'.\n";

    ## input
    ## AM-WGBS-001_R1.fq.gz
    ## AM-WGBS-001_R2.fq.gz
    (my $prefix = $fq1) =~ s/(?:\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//;

    ## Final output of standard bismark call is:
    ## -----------------------------------------
    ## AM-WGBS-001_R1_bismark_bt2_pe.bam
    ## AM-WGBS-001_R1_bismark_bt2_PE_report.txt
    ##               ^^^^^^^^^^^^^^^^^^^^^^^^^^
    my $bam_file = "${basename}.bam";
    my $report   = "${bam_file}.alignment.report.txt";

    ## in case Bismark has changed naming :-)
    unless(rename("${prefix}_bismark_bt2_pe.bam",$bam_file)) {
        warn "[WARNING] $!\n";
        $bam_file = glob "*.bam";
    }

    rename("${prefix}_bismark_bt2_PE_report.txt",$report)
        or warn "[WARNING] $!\n";

    return $bam_file;
}


### ===================================================================
# Name .....: run_fastqc
# Task .....: runs FastQC on a single fastq file
# Parameter : file name, output folder
# Output ...: one HTML file and one ZIP archive
# Returns ..: -
### ===================================================================
sub run_fastqc {
    my($fastq,$out_dir)=@_;
    my $time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);

    my $cmd = "$progs{fastqc}";
    $cmd .= " --outdir $out_dir";
    $cmd .= " --threads $p_threads";
    $cmd .= " --dir $cwd";
    $cmd .= " $fastq";

    print STDERR "[$time_stamp] Running FastQC on $fastq ...\n";

    log_command("FastQC", $cmd);

    system("$cmd") 
        and die "[ERROR] Some problem running FastQC on '$fastq'.\n";
}


### ===================================================================
# Name .....: 
# Task .....: 
# Parameter : 
# Output ...:
# Returns ..:
# Info .....: http://cutadapt.readthedocs.io/en/stable/recipes.html#piping-paired-end-data
### ===================================================================
sub trim_swift_adaptor_PE {
    my ($fq1,$fq2,$adaptor,$sample,$length)=@_;
    my $time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);

    my ($prefix,$postfix);
    my ($new_fq1,$new_fq2);

    if($fq1=~/^(.+)(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$/) {
        $prefix  = $1;
        $postfix = $2;
        $new_fq1 = "${prefix}_${adaptor}-swift-trim${postfix}";
    }

    if($fq2=~/^(.+)(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$/) {
        $prefix  = $1;
        $postfix = $2;
        $new_fq2 = "${prefix}_${adaptor}-swift-trim${postfix}";
    }

    ## AM-WGBS-001_adp-trim.log
    my $log = $new_fq1 . ".log";
    $log =~ s/(_R1|\.fq|\.gz)//g;

    my $cmd = $progs{cutadapt};
    $cmd .= " --quality-cutoff 20";
    $cmd .= " --overlap 5";
    $cmd .= " --minimum-length 25";
    $cmd .= " --cores " . ($c_threads-$p_threads);
    $cmd .= " --adapter $adaptors{$adaptor}";
    $cmd .= " -A $adaptors{$adaptor}";
    $cmd .= " --interleaved";
    $cmd .= " $fq1";
    $cmd .= " $fq2";
    $cmd .= " 2>.log_a";

    $cmd .= " |";

    $cmd .= " $progs{cutadapt}";
    $cmd .= " --interleaved";
    $cmd .= " --minimum-length 25";
    $cmd .= " --cores $p_threads";
    $cmd .= " --cut $length --cut -$length";
    $cmd .= " -U $length -U -$length";
    $cmd .= " --output $new_fq1";
    $cmd .= " --paired-output $new_fq2";
    $cmd .= " -";
    $cmd .= " &>.log_b";

    print STDERR "[$time_stamp] Running Cutadapt ($adaptor trimming + Swift) on $sample ...\n";

    log_command("Trimming adaptors (Swift, PE)", $cmd);

    system("$cmd")
        and die "[ERROR] Some problem running cutadapt:\n cmd = $cmd\n";

    foreach my $l (".log_a",".log_b") {
        system("cat $l >> $log");
    }

    return $new_fq1, $new_fq2;
}


### ===================================================================
# Name .....: trim_adaptor_PE
# Task .....: trim off Illumina adaptors from the 3' end
# Parameter : fq1/2, adaptor, sample name
# Output ...: trimmed files, log
# Returns ..: new file names (fq1,fq2,log file)
### ===================================================================
sub trim_adaptor_PE {
    my ($fq1,$fq2,$adaptor,$sample)=@_;
    my $time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);

    my ($prefix,$postfix);
    my ($new_fq1,$new_fq2);

    if($fq1=~/^(.+)(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$/) {
        $prefix  = $1;
        $postfix = $2;
        $new_fq1 = "${prefix}_${adaptor}-trim${postfix}";
    }

    if($fq2=~/^(.+)(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$/) {
        $prefix  = $1;
        $postfix = $2;
        $new_fq2 = "${prefix}_${adaptor}-trim${postfix}";
    }

    ## AM-WGBS-001_adp-trim.log
    my $log = $new_fq1 . ".log";
    $log =~ s/(_R1|\.fq|\.gz)//g;

    my $cmd = $progs{cutadapt};
    $cmd .= " --quality-cutoff 20";
    $cmd .= " --overlap 5";
    $cmd .= " --minimum-length 25";
    $cmd .= " --cores $c_threads";
    $cmd .= " --adapter $adaptors{$adaptor}";
    $cmd .= " -A $adaptors{$adaptor}";
    $cmd .= " --output $new_fq1";
    $cmd .= " --paired-output $new_fq2";
    $cmd .= " $fq1";
    $cmd .= " $fq2";

    print STDERR "[$time_stamp] Running Cutadapt (adaptor trimming) on $sample ...\n";

    log_command("Trimming adaptors (PE)", $cmd);

    system("$cmd &> $log")
        and die "[ERROR] Some problem running cutadapt:\n cmd = $cmd\n";

    return $new_fq1, $new_fq2;
}


### ===================================================================
# Name .....: trim_swift_PE
# Task .....: hard-clip PE data on both ends (from Swift Libraries)
# Parameter : fq1/2, length to clip, sample name
# Output ...: trimmed files, log
# Returns ..: new file names (fq1,fq2,log file)
### ===================================================================
sub trim_bases_PE {
    my ($fq1,$fq2,$length,$sample)=@_;
    my $time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);

    my ($prefix,$postfix);
    my ($new_fq1,$new_fq2);

    if($fq1=~/^(.+)(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$/) {
        $prefix  = $1;
        $postfix = $2;
        $new_fq1 = "${prefix}_${length}bp-trim${postfix}";
    }

    if($fq2=~/^(.+)(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$/) {
        $prefix  = $1;
        $postfix = $2;
        $new_fq2 = "${prefix}_${length}bp-trim${postfix}";
    }
    
    ## AM-WGBS-001_adp-trim_10bp-trim.log
    my $log = $new_fq1 . ".log";
    $log =~ s/(_R1|\.fq|\.gz)//g;

    my $cmd = $progs{cutadapt};
    $cmd .= " --minimum-length 25";
    $cmd .= " --cores $c_threads";
    $cmd .= " --cut $length --cut -$length";
    $cmd .= " -U $length -U -$length";
    $cmd .= " --output $new_fq1";
    $cmd .= " --paired-output $new_fq2";
    $cmd .= " $fq1";
    $cmd .= " $fq2";

    print STDERR "[$time_stamp] Running Cutadapt (fixed-bases trimming) on $sample ...\n";

    log_command("Trimming $length bases", $cmd);

    system("$cmd &> $log")
        and die "[ERROR] Some problem running cutadapt:\n cmd = $cmd\n";

    return $new_fq1, $new_fq2, $log;
}



__DATA__

### ===================================================================
# Name .....: 
# Task .....: 
# Parameter : 
# Output ...:
# Returns ..:
### ===================================================================
