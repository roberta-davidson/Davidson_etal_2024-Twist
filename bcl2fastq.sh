#!/bin/bash

module load Singularity/3.10.5 

while getopts r:o:s:m:t: flag
do
    case "${flag}" in
        r) RUNFOLDER=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        s) SAMPLESHEET=${OPTARG};;
        m) USE_BASE_MASK=${OPTARG};;
        t) TILES=${OPTARG};;
    esac
done

if [ -z ${RUNFOLDER} ] || [ -z ${OUTPUT_DIR} ] || [ -z ${SAMPLESHEET} ]; then
    echo "Did not provide one of:"
    echo "Runfolder (-r): ${RUNFOLDER}"
    echo "Runfolder (-o): ${OUTPUT_DIR}"
    echo "Runfolder (-s): ${SAMPLESHEET}"
    exit 1
fi

if [ ! -z $USE_BASE_MASK ]; then
    base_mask_line="--use-bases-mask ${USE_BASE_MASK}"
fi

if [ ! -z $TILES ]; then
    tiles_line="--tiles ${TILES}"
fi

mkdir -p ${OUTPUT_DIR}

singularity exec -B /hpcfs/ /hpcfs/groups/acad_users/containers/bcl2fastq_v2.20.0.sif bcl2fastq \
    --runfolder-dir ${RUNFOLDER} \
    --output-dir ${OUTPUT_DIR} \
    --sample-sheet ${SAMPLESHEET} \
    ${base_mask_line} ${tiles_line}
