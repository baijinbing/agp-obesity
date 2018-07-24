#!/bin/bash

for ((seed=1;seed<=42;seed++))
do
        echo "#!/bin/bash" >> script.out
        echo "#$ -N wt" >> script.out
        echo "#$ -cwd" >> script.out
        echo "#$ -j y" >> script.out
        echo "R CMD BATCH --no-save --no-restore \"--args seed=${seed} \" ./LDM_wt.R" >> script.out 
        chmod +x script.out
        qsub -q gene.q ./script.out
        rm -rf script.out
done

