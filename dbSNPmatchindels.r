#!/usr/bin/env Rscript
##Author: Or Yaacov 
##Date: 10/22/2019

args=commandArgs(TRUE)

filePath = args[1] # path to your input vcf / vcf.gz
outPut = args[2] # path to your desired output vcf
notfound= args[3] # path to a desired output vcf for the indels that could not be validated on dbSNP and hence removed. 
dbsnp_file=args[4] # path to a dbSNP vcf that contains only the indels. You can make one by downloading the full file, All_20170710.vcf.gz, from dbSNP and then:
# zgrep "VC=DIV" All_20170710.vcf.gz | awk '{print $1,$2,$3,$4,$5,$6}' | gzip > indels.vcf.gz

library(data.table)

#Read the input file
read1 <- fread(filePath)
dbsnp<- fread(dbsnp_file)

# find the SNP based on position or Rsid:
findsnp <- function(Chr, Pos, Rs= 0) {
    #find first match
    idx<-which((dbsnp[,1] == Chr) & (dbsnp[,2] == Pos ) & (dbsnp[,3] == Rs))[1]
    dbRef<- dbsnp[idx,4]
    dbAlt<- dbsnp[idx,5]
    if (is.na(dbRef) | is.na(dbAlt)) {
    dbRef<- "nnnn"
    dbAlt<- "nnnn"   
    }
    return(c(dbRef, dbAlt))
}

#Convert the data
fix<- character()
err<- c()

print(dim(read1))

for(i in 1:nrow(read1)) {
    row <- read1[i,]
    alt <- read1[[i,5]]
    ref <- read1[[i,4]]
    rs <- read1[[i,3]]
    pos <- read1[[i,2]]
    chr <- read1[[i,1]]
    if ((nchar(alt) > 1) | (nchar(ref) > 1)) {
        snp<- findsnp(chr, pos, rs)
        if ((nchar(snp[[1]]) == nchar(alt)) & (nchar(snp[[2]]) == nchar(ref))) {
            read1[i,5] <- snp[[1]]
            read1[i,4] <- snp[[2]]
            fix<- append(fix, i)
            }
        else if ((nchar(snp[[2]]) == nchar(alt)) & (nchar(snp[[1]]) == nchar(ref))) {
            read1[i,5] <- snp[[2]]
            read1[i,4] <- snp[[1]]
            fix<- append(fix, i)
            } 
        else {
            err<- append(err, i)
            }
        }
}

leftout<- read1[err,]
read1<- read1[-err,]

print(dim(read1))
print("Fixed lines:")
print(fix)
print("Could not be found on dbSNP- removed lines: ")
print(err)


#Save output:
write.table(read1,file=outPut, sep = "\t", quote = FALSE, col.names = F, row.names = F)
write.table(leftout,file=notfound, sep = "\t", quote = FALSE, col.names = F, row.names = F)
print("Done!")
