library(Rsamtools)


samfile = "../../Ct1712017_1012_luojie_Z.target.realigned.sam"

what = c("rname","stand","pos","qwidth","seq")
?scanBamWhat()
param= scanBamWhat()
bam= scanBam(bamfile,param = param,)

install.packages("RUnit")
BiocGenerics:::testPackage('Rsamtools')

bamWhat(p4)
bamWhich(p3)
## subset of fields
bamfile = "../../Ct1712017_1012_luojie_Z.target.realigned.bam"
p3 <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth","flag"))

## tags; NM: edit distance; H1: 1-difference hits
p4 <- ScanBamParam(tag=c("NM", "XA"),what=scanBamWhat())
bam4 <- scanBam(bamfile, param=p4,asMates=TRUE)

bam4[[1]][["groupid"]] 
a =bam4[[1]][["tag"]]

a$XA[!is.na(a$XA)]
str(bam4[[1]][["tag"]])

## tagFilter
p5 <- ScanBamParam(tag=c("NM", "H1"), tagFilter=list(NM=c(2, 3, 4)))
bam5 <- scanBam(fl, param=p5)
table(bam5[[1]][["tag"]][["NM"]])

## flag utils
flag <- scanBamFlag(isUnmappedQuery=FALSE, isMinusStrand=TRUE)

p6 <- ScanBamParam(what="flag")
bam6 <- scanBam(fl, param=p6)
flag6 <- bam6[[1]][["flag"]]
head(bamFlagAsBitMatrix(flag6[1:9]))
colSums(bamFlagAsBitMatrix(flag6))
flag
bamFlagAsBitMatrix(flag)


for (i in c(1:22,"X")) {
  sprintf("java -Xmx10g -Djava.io.tmpdir=/tmp -jar %s -nct 2 -T MuTect2 -R %s -I:tumor %s -o vcf/%s.chr%s.MuTect2.vcf -hets 0.0001 -L %s","GATK","FA","sample","sample",i,i) %>% print()
}
