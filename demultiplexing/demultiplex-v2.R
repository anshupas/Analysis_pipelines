# library(optparse)
# 
# option_list = list(
#   make_option(c("-r", "--rootFolder"), type="character", default=NULL, action = "store",
#               help="path of folder containing fastq files", metavar="character"),
#   make_option(c("-l", "--library"), type="character", default= NULL, 
#               help="sample id", metavar="character"),
#   make_option(c("-b", "--barcode"), type="character", default=NULL, 
#               help="barcode to be scanned", metavar="character"),
#   make_option(c("-s", "--sample"), type="character", default=NULL, 
#               help="sample name", metavar="character")
# ); 
# 
# opt_parser = OptionParser(option_list=option_list)
# opt = parse_args(opt_parser)
# print(opt$barcode)

args = commandArgs(trailingOnly = TRUE)
print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])

if (is.null(args[3])) {
  #print_help(opt_parser)
  stop("Barcode must be supplied.\n", call.=FALSE)
}

if (is.null(args[1])) {
  #print_help(opt_parser)
  stop("Specify the folder path containing fastq.gz files.\n", call.=FALSE)
}

if(is.null(args[4])){
  stop("Specify the sample name.\n", call.=FALSE)
}


library(gtools)
library(R.utils)
library(data.table)
library(ShortRead)


# generate sequences of same length as barcode with a Leveinstein distance of at max 1 (allowing only one mismatch)

bcData = cbind(c(1:8),c("GATATTAT","GCCGCTTT","GCATTCGC","GCCGCCAG","CGAGTTGG","AAATTTTC","AGCGTTTG","GGCGTCAG"))
colnames(bcData) = c("index","barcode")
bcData = data.frame(bcData)
bcLen = unique(sapply(bcData$barcode,function(x){
  nchar(as.character(x))
}))


bases=c('A','T','G','C')
kmers = permutations(n=length(bases),r=8,v=bases,repeats.allowed=T)
strList = apply(kmers,1,paste,collapse="")
#barcode = "GATATTAT"
barcode = args[3]
if(grepl("-",barcode)){
  bInd = unlist(strsplit(barcode,split = "-"))
  if(bInd[2] <= nrow(bcData)){
    bcList = bcData[bInd[1]:bInd[2],"barcode"]
  }else{
    bcList = NULL
    stop("Index higher than size of barcode list")
  }
  
}else{
  bcList = bcData[as.integer(barcode),"barcode"]
}

pVal = unique(as.character(sapply(bcList,function(bc){
  combSeq = setdiff(strList,bc)
  levDist = as.data.frame(adist(bc,combSeq))
  colnames(levDist) = combSeq
  bcVal = as.character(bc)
  probMatchv1 = colnames(levDist[which(levDist==1)])
  probMatchv2 = unlist(lapply(c(1:8), function(x){
    newBc = as.character(bc)
    substr(newBc,x,x)="N"
    newBc
  }))
  probMatch = c(bcVal,probMatchv1,probMatchv2)

  #pVal = paste(probMatch,collapse = "|")
})))

patternVal = paste(pVal,collapse = "|")

print("All probable patterns generated")

######################## Reading fastq files ########################
print("Reading fastq files .....")

readsR1 = readFastq(paste0(args[1],"/",args[2],"_R1_001.fastq.gz"))
readsR2 = readFastq(paste0(args[1],"/",args[2],"_R2_001.fastq.gz"))


######################## Paired filtering for the barcode ########################

print("Filter barcodes .....")

filterBarcode=srFilter(function(x){
  substr(sread(readsR1[x]),1,bcLen) %in% pVal | substr(sread(readsR2[x]),1,bcLen) %in% pVal
},name="Barcode Filter")

filtOb = filterBarcode(1:length(readsR1))
#filterR2 = filterBarcode(readsR2)

validIndex = which(filtOb@.Data  == TRUE)

validR1 = readsR1[filtOb]
validR2 = readsR2[filtOb]

### Trimming the barcode sequence from the sequence and also trim the quality line accordingly

print("Trimming reads and quality .....")

validR1@sread = DNAStringSet(substr(sread(validR1),bcLen+1,unique(width(sread(validR1)))))
validR1@quality = FastqQuality(substr(quality(quality(validR1)),bcLen+1,
                                      unique(width(quality(quality(validR1))))))

validR2@sread = DNAStringSet(substr(sread(validR2),bcLen+1,unique(width(sread(validR2)))))
validR2@quality = FastqQuality(substr(quality(quality(validR2)),bcLen+1,
                                      unique(width(quality(quality(validR2))))))

######################## Writing fastq file ########################

print("Writing fastq files .....")

if(grepl(pattern = "mESC",args[4])){
  orgInfo = "mouse"
}else{
  orgInfo = "human"
}

writeFastq(validR1, paste0(args[5],orgInfo,"/",args[2],"_",args[4],"_R1_001","_",".fastq.gz"), compress=TRUE)
writeFastq(validR2, paste0(args[5],orgInfo,"/",args[2],"_",args[4],"_R2_001","_",".fastq.gz"), compress=TRUE)
