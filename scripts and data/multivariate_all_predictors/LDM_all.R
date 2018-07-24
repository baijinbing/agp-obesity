options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

if(length(args)==0){
    print("No arguments supplied")
    
    seed = 1
    
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

source('../LDM_fun.r')

otu.tab <- read.csv("../otu_l6_267.csv", header=TRUE,row.names=1, sep=",")
meta <- read.csv("../metadata_267.csv", header=TRUE,row.names=1, sep=",")
# tree <- read.tree("../tree.nwk")

dim(otu.tab)
dim(meta)

# LDM 1 - all as predictor
    
res1.ldm <- ldm(formula=otu.tab|(SEX+AGE_CAT+ANTIBIOTIC_RECODE+PROBIOTIC_RECODE) ~ DIET_RECODE+EXERCISE_RECODE+BMI_CAT_RECODE+WEIGHT_CHANGE,
                data=meta, 
                dist.type="Bray-Curtis",
                test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1,
                out.stat=TRUE, out.prefix="tmp_files/tmp", 
                n.perm.max=1000, seed=seed)
