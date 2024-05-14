#!/bin/bash

module load singularity
#assign variable to dataset called in command
in1=$1

# for slurm
cat ${in1}.trees.txt

echo "Running f4 tests"

#Run on every quadruple combination
singularity exec -B /hpcfs/ /hpcfs/groups/acad_users/containers/admixtools_7.0.2--h6a739c9_4.sif qpDstat \
-p <(echo "indivname:    ${in1}.ind  
snpname:      ${in1}.snp
genotypename: ${in1}.geno
poplistname: ${in1}.poplist.txt
popfilename: ${in1}.trees.txt
inbreed: 	YES
printsd:  YES
f4mode:   YES") \
> ${in1}.qpDstat_f4.out

#clean up output file for R
grep "result" ${in1}.qpDstat_f4.out | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > ${in1}.qpDstat_f4.R.out
