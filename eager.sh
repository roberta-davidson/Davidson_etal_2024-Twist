#!/bin/bash

INPUT_FILE=input_file_Twist_benchmark.tsv

FASTA=/hpcfs/groups/acad_users/Refs/Homo_sapiens/GATK/b37/human_g1k_v37_decoy.fasta
FASTA_INDEX=/hpcfs/groups/acad_users/Refs/Homo_sapiens/GATK/b37/human_g1k_v37_decoy.fasta.fai
SEQ_DICT=/hpcfs/groups/acad_users/Refs/Homo_sapiens/GATK/b37/human_g1k_v37_decoy.dict
BWA_INDEX=/hpcfs/groups/acad_users/Refs/Homo_sapiens/GATK/b37

nextflow run nf-core/eager -c phoenix.config \
	-r 2.4.6 \
	-with-singularity \
	--outdir 'Twist_filtered/' \
	--input ${INPUT_FILE} \
	--fasta ${FASTA} \
	--fasta_index ${FASTA_INDEX} \
  --bwa_index ${BWA_INDEX} \
  --seq_dict ${SEQ_DICT} \
	--complexity_filter_poly_g \
  --mergedonly \
	--mapper 'bwaaln' \
	--bwaalnn 0.01 \
	--bwaalno 2 \
  --bwaalnl 1024 \
	--dedupper 'markduplicates' \
	--clip_readlength 30 \
  --clip_min_read_quality 20 \
	--run_bam_filtering \
  --bam_mapping_quality_threshold 25 \
  --bam_unmapped_type 'discard' \
	--run_trim_bam \
	--bamutils_clip_double_stranded_half_udg_left 2 \
	--bamutils_clip_double_stranded_half_udg_right 2 \
	--run_sexdeterrmine \
	--sexdeterrmine_bedfile '/<path>/Twist_targets_march2024.pos' \
	--run_nuclear_contamination \
	--contamination_chrom_name 'X' \
	--run_mtnucratio \
  --mtnucratio_header 'MT' \
	--run_genotyping \
	--genotyping_source 'trimmed' \
	--genotyping_tool 'pileupcaller' \
	--pileupcaller_bedfile '/<path>/Twist_targets_march2024.pos' \
	--pileupcaller_snpfile '/<path>/Twist_targets_march2024.snp' \
	--pileupcaller_method 'randomHaploid' 
