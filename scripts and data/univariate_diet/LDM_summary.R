source('../LDM_fun.r')

res3.ldm <- test.ldm.fromfile(in.prefix="tmp_files/tmp", n.perm.available=42000,
                              test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1)  

res3.ldm$n.perm.completed
res3.ldm$global.tests.stopped
res3.ldm$otu.tests.stopped
res3.ldm$p.global.omni
otu.omni.var1 = sort(which(res3.ldm$q.otu.omni[1,] < 0.1))
(n.otu.omni.var1 = length(otu.omni.var1)) 
otu.omni.var1


