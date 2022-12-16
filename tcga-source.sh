#!/usr/bin/env bash

#################### SETUP ####################

# print status
s="START"
PS4='+	${s}	${LINENO}	$( date "+%Y.%m.%d %R") | '
set -x 

# command fails -> exit
# set -o errexit
# undeclared variables -> exit
# set -o nounset
# mysqldump fails -> exit
# set -o pipefail

###############################################

wd=$(dirname $(realpath $0))

main(){
	01_kallisto_prep
	01_kallisto_run
	02_tx2gene
}

01_kallisto_prep(){
	local my_dir="${wd}/01-kallisto-dir"
	local my_dir_tmp="/home/scratch60/rbarreiro_magoh_kallisto/index"

	mkdir -p $my_dir_tmp
	# not sure if i'll use this or /home/users/tmiller//projects_current/RTC/expression/geneseq.idx
	kallisto index -i $my_dir_tmp/gencode.v29.transcripts.kallisto_index \
		${wd}/resources/gencode/*.fa.gz

	local my_dir_tmp="/home/scratch60/rbarreiro_magoh_kallisto/gdc_client"
	wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.0_Ubuntu_x64-py3.7_0.zip \
		-P $my_dir_tmp

	unzip $my_dir_tmp/*.zip

	mkdir $my_dir/output
}

01_kallisto_run(){
	local my_dir="${wd}/01-kallisto-dir"
	local my_dir_tmp="/home/scratch60/rbarreiro_magoh_kallisto2/kallisto_run"
	mkdir -p $my_dir_tmp

	local n_lines=$(wc -l < ${wd}/00-TCGA_manifest/manifest_chunks/manifest.txt)

	
	for (( i_lines = 1; i_lines <= $n_lines; i_lines++ )); do
		#get file id
		local file_id=$(awk '(NR=='$i_lines'){print $1}' ${wd}/00-TCGA_manifest/manifest_chunks/manifest_${MY_INDEX}.txt)
		mkdir -p ${my_dir_tmp}/${file_id}
		
		#create manifest	
		mkdir -p ${my_dir_tmp}/${file_id}/00_manifest
		(cat ${wd}/00-TCGA_manifest/manifest_chunks/header.txt
		awk '(NR=='$i_lines'){print}' ${wd}/00-TCGA_manifest/manifest_chunks/manifest_${MY_INDEX}.txt) \
		 > ${my_dir_tmp}/${file_id}/00_manifest/manifest.txt

		#download bam
		mkdir -p ${my_dir_tmp}/${file_id}/01_bamfile
		/home/users/rbarreiro/scripts/gdc-client download \
			-m ${my_dir_tmp}/${file_id}/00_manifest/manifest.txt \
			-t $TOKEN \
			-d ${my_dir_tmp}/${file_id}/01_bamfile

		#bam to fastq
		mkdir -p ${my_dir_tmp}/${file_id}/02_fastq
		bam_to_fastq $(find ${my_dir_tmp}/${file_id}/01_bamfile/ -name "*.bam") ${my_dir_tmp}/${file_id}/02_fastq

		#kallisto
		mkdir -p ${my_dir_tmp}/${file_id}/03_kallisto
		kallisto quant \
			-i /home/users/tmiller/projects_current/RTC/expression/geneseq.idx \
			-t 10 \
			-b 100 \
			-o ${my_dir_tmp}/${file_id}/03_kallisto ${my_dir_tmp}/${file_id}/02_fastq/*.fastq.gz
	
		#send to output
		mkdir -p ${my_dir}/output/${file_id}/
		mv ${my_dir_tmp}/${file_id}/03_kallisto/* ${my_dir}/output/${file_id}/

		#clean
		rm -rf ${my_dir_tmp}/${file_id}
	done
}

bam_to_fastq(){
	local bamfile=$1; 
	local outdir=$2; 
	local my_filename=$(basename $bamfile | sed 's/\..*//')

	/home/tools/manual/samtools-1.9/samtools fastq \
		-1 "${outdir}/${my_filename}_R1.fastq.gz" \
		-2 "${outdir}/${my_filename}_R2.fastq.gz" \
		-0 /dev/null \
		$bamfile 
}

02_tx2gene (){
  s='02-tx2gene'

  my_dir="${wd}/02-tx2gene"
  mkdir -p ${my_dir}

  # INPUT #####################################################################
  mkdir -p ${my_dir}/input/
  ln -sf ${wd}/01-kallisto/output/* ${my_dir}/input

  mkdir -p ${my_dir}/resources/gencode
  ln -sf $(realpath ${wd}/resources/gencode/gencode.v29.annotation.gtf) ${my_dir}/resources/gencode/

  mkdir -p ${my_dir}/scripts
  ln -sf ${wd}/scripts/kallisto_tx2gene.R  ${my_dir}/scripts/

  # TPM #######################################################################

  mkdir -p ${my_dir}/tmp/tx2gene
  create_dict_table ${my_dir}/resources/gencode/gencode.v29.annotation.gtf \
  	> ${my_dir}/tmp/tx2gene/gencode_tsx2gene.tsv

  # OUTPUT ####################################################################
  mkdir -p ${my_dir}/output/counts/
  
  for file in $(find ${my_dir}/input -maxdepth 1 | tail -n +2 | sort); do
    abundance_file="${file}/abundance.h5"
    sample_name=$(basename ${file} | sed 's/_R1//')
    output="${my_dir}/output/counts/${sample_name}_tx2gene_counts.tsv"
    
    Rscript ${my_dir}/scripts/kallisto_tx2gene.R \
      --tx2gene ${my_dir}/tmp/tx2gene/gencode_tsx2gene.tsv \
      --kallisto_abundance $abundance_file \
      --sample $sample_name \
      --output $output 

      set +x; jobslimit 20; set -x

  done  
  wait
}

create_dict_table(){
  printf "TXNAME\tGENEID\n"
  cat $1 \
    | grep -v "^#" \
    | sed -En 's/.*gene_id "([^"]*).*transcript_id "([^"]*).*/\2\t\1/p' \
    | sort -u
}

main
wait
