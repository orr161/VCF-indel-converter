library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

#Enter input file name and output name/path in "" 
filePath = "hyphen.vcf"
outPut = "leading.vcf"

#Find base/seq function for hg19
findbase <- function(chr, pos, len=0) {
  letter <- toString(Hsapiens[[chr]][pos:(pos+len)])
  return(letter)
}

#Read the input file
read1 <- read.table(file = filePath, col.names = c("chr", "pos","name", "ref", "alt")
                    , colClasses=c("ref"="character", "alt"="character"))

#Converts the format
for(i in 1:nrow(read1)) {
    row <- read1[i,]
    alt <- read1[i,5]
    ref <- read1[i,4]
    pos <- read1[i,2]
    chr <- read1[i,1]
    if (alt == "-") {
        Len= 1
        leadingPos = (pos-1)
        read1[i,2] <- leadingPos
        leading = findbase(toString(chr), leadingPos)
        read1[i,5] = leading
        read1[i,4] = paste(leading, ref, sep="")
        }
    else if (ref == "-") {
        Len= 1
        leadingPos = (pos-1)
        read1[i,2] = leadingPos
        leading = findbase(toString(chr), leadingPos)
        read1[i,4] = leading
        read1[i,5] = paste(leading, alt, sep="")
    }
}


#Save output:
write.table(read1,file=outPut, sep = "\t", quote = FALSE, col.names = F, row.names = F)
print("Done!")