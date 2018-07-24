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
otu.omni.var2 = sort(which(res3.ldm$q.otu.omni[2,] < 0.1))
(n.otu.omni.var2 = length(otu.omni.var2)) 
otu.omni.var2
otu.omni.var3 = sort(which(res3.ldm$q.otu.omni[3,] < 0.1))
(n.otu.omni.var3 = length(otu.omni.var3))
otu.omni.var3
otu.omni.var4 = sort(which(res3.ldm$q.otu.omni[4,] < 0.1))
(n.otu.omni.var4 = length(otu.omni.var4))
otu.omni.var4
