#!/bin/bash

        echo "#!/bin/bash" >> script.out
        echo "#$ -N sum_diet" >> script.out
        echo "#$ -cwd" >> script.out
        echo "#$ -j y" >> script.out
        echo "R CMD BATCH ./LDM_summary.R" >> script.out 
        chmod +x script.out
        qsub -q gene.q ./script.out
        rm -rf script.out

