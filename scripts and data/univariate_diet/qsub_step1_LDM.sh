#!/bin/bash

for ((seed=11;seed<=42;seed++))
do
        echo "#!/bin/bash" >> script.out
        echo "#$ -N diet" >> script.out
        echo "#$ -cwd" >> script.out
        echo "#$ -j y" >> script.out
        echo "R CMD BATCH --no-save --no-restore \"--args seed=${seed} \" ./LDM_diet.R" >> script.out 
        chmod +x script.out
        qsub -q gene.q ./script.out
        rm -rf script.out
done

