otu.file <- "LDM_multi_BMI.txt" # (modify)


#------------------------
# otu.table, freq
#------------------------

otu.table <- read.csv("../otu_l6_267.csv", header=TRUE,row.names=1, sep=",")

freq.table <- t( scale( t(otu.table), center=FALSE, scale=rowSums(otu.table) ) )
freq <- colMeans(freq.table) # 416


#------------------------
# tax.matrix
#------------------------

otu.names <- colnames(otu.table)
length(otu.names) # 416
#otu.names <- gsub("[.][.]", ".", x=otu.names)
#otu.names <- gsub("__[.]", "__", x=otu.names) # c__.Chloracidobacteria.
otu.names <- gsub("[.]__", "",      x=otu.names)
otu.names <- gsub("[.]p__", ";p__", x=otu.names)
otu.names <- gsub("[.]c__", ";c__", x=otu.names)
otu.names <- gsub("[.]o__", ";o__", x=otu.names)
otu.names <- gsub("[.]f__", ";f__", x=otu.names)
otu.names <- gsub("[.]g__", ";g__", x=otu.names)

tax.list = strsplit(otu.names, ";")
lenmax = max(sapply(tax.list, length)) # 8
fun = function(x) c(x, rep(NA, lenmax - length(x)))
tax.list = lapply(tax.list, fun)
tax.matrix = matrix(data=unlist(tax.list), nrow=length(tax.list), ncol=lenmax, byrow=TRUE)
dim(tax.matrix) # 416  6


#------------------------
# otu.detected (output from LDM)
#------------------------

otu.detected.full = read.table(otu.file, as.is=TRUE, header=FALSE)[,1]
otu.detected = as.numeric(otu.detected.full[seq(from=2, to=length(otu.detected.full), by=2)])
length(otu.detected) # 120


#------------------------
# tax.matrix.detected
#------------------------

tax.matrix.detected = tax.matrix[otu.detected,]
dim(tax.matrix.detected) # 120 6
(prop.overall = nrow(tax.matrix.detected)/nrow(tax.matrix)) # 0.2884615


#------------------------
# table
#------------------------
l = 2 #(modify)
filter.min.otu.bytaxa = 5

n.otu.bytaxa = as.matrix(table(tax.matrix[,l]))
taxa.names = rownames(n.otu.bytaxa)
n.taxa = length(taxa.names)

w = (matrix(rep(tax.matrix[,l], n.taxa), ncol=n.taxa)
    == matrix(rep(taxa.names, dim(tax.matrix)[1]), ncol=n.taxa, byrow=TRUE))
ww = (matrix(rep(tax.matrix.detected[,l], n.taxa), ncol=n.taxa)
      ==matrix(rep(taxa.names, dim(tax.matrix.detected)[1]), ncol=n.taxa, byrow=TRUE))
n.otu.detected = apply(ww, 2, sum, na.rm=TRUE)
freq.otu.bytaxa = apply(w, 2, function(x)sum(freq[x], na.rm=TRUE))
freq.otu.detected = apply(ww, 2, function(x)sum(freq[otu.detected[x]], na.rm=TRUE))
prop.otu.detected = n.otu.detected/n.otu.bytaxa

tab=cbind(n.otu.bytaxa, n.otu.detected, freq.otu.bytaxa, freq.otu.detected, prop.otu.detected)
tab=tab[order(tab[,5], decreasing=TRUE),]
(tab=tab[which(tab[,1]>filter.min.otu.bytaxa),])


#------------------------
# enrichment analysis
#------------------------

# test whether the proportion of detected OTUs varied significantly? 
set.seed(123)
fisher.test(cbind(tab[,2], tab[,1]-tab[,2]),simulate.p.value = TRUE) # p=0.1719

# logistic regression with Firth correction
t = c(1,4,5,6,7) # taxa of interest (modify)

Y = c(rep(1,sum(tab[-t,2])), rep(0,sum(tab[-t,1])-sum(tab[-t,2])))
ti = matrix(NA, ncol=length(t), nrow=sum(tab[,1]))
for (i in 1:length(t)) {
    Y = c(Y, rep(1,tab[t[i],2]), rep(0,tab[t[i],1]-tab[t[i],2]))
    ti[,i] = c(rep(1, tab[t[i],1]), rep(0, sum(tab[-t[i],1])))
}

library(logistf)
lf = logistf(Y ~ ti, firth=TRUE)

p.value = rep("*", nrow(tab))
p.value[t] = lf$prob[-1]


data.frame(tab[,c(1,2,5)], p.value)
