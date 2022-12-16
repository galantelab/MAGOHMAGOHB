#!/usr/bin/env bash

#################### SETUP ####################

# print status
s="START"
PS4='+	${s}	${LINENO}	$( date "+%Y.%m.%d %R") | '
set -x 

# command fails -> exit
set -o errexit
# undeclared variables -> exit
# set -o nounset
# mysqldump fails -> exit
set -o pipefail

# jobs wait
source /home/users/rbarreiro/scripts/bash/jobs_limit.sh
# HOW TO WAIT 
# put inside a loop:
# set +x; jobslimit {N}; set -x

###############################################

wd=$(dirname $(realpath $0))

main(){
	# 00_input
	# 01_fastqc
	# 02_merge_and_trim
	# 03_kallisto
	# 04_tximport
	# 05_deseq2_pre
	# 05_deseq2_run
  # 06_STAR
  07_rMATS
}

00_input(){
	s="INPUT"

	mkdir -p ${wd}/00-input/{KD,C}
	ln -sf /home/projects2/downloads/penalva/PHF5a/Penalva_L_08042016/U251_C*  ${wd}/00-input/C
	ln -sf /home/projects2/downloads/penalva/PHF5a/Penalva_L_08042016/U251_KD* ${wd}/00-input/KD

	for file in $(find ${wd}/00-input/ -name "*.gz"); do
		mv  $file $(sed 's/_S[^_]*_/_/' <<< $file);
	done
}

01_fastqc(){
	s="FASTQC"

	my_dir=${wd}/01-fastqc/
	mkdir -p ${my_dir}
 	fastqc $(find ${wd}/00-input/ -name "*.gz") -o ${my_dir} -t 20 -f fastq
}

02_merge_and_trim(){
	my_dir=${wd}/02-merge_and_trim/

	my_tmp_dir="/home/scratch60/magoh_u251_merge/"

	mkdir -p $my_tmp_dir
	for file in $(find 00-input/ -name "*.gz"); do
		new_name=$(echo $file | sed -E 's/_L00._(..)_001/_\1/' | sed 's|.*/||')
		cat $file >> ${my_tmp_dir}/$new_name
	done

  for file in $(find $my_tmp_dir -name *_R1.fastq.gz); do
    my_output_dir=${my_dir}/$(basename ${file/_R1.fastq.gz})
    mkdir -p ${my_output_dir}
    trim_galore \
    	-q 20 \
    	--fastqc \
    	--illumina \
    	--gzip \
    	--output_dir $my_output_dir \
    	--paired ${file} ${file/_R1.fastq.gz/_R2.fastq.gz} &
  done
  wait

}

02_merge_and_trim(){
	my_dir=${wd}/02-merge_and_trim/

	my_tmp_dir="/home/scratch60/magoh_u251_merge/"

	mkdir -p $my_tmp_dir
	for file in $(find 00-input/ -name "*.gz"); do
		new_name=$(echo $file | sed -E 's/_L00._(..)_001/_\1/' | sed 's|.*/||')
		cat $file >> ${my_tmp_dir}/$new_name
	done

  for file in $(find $my_tmp_dir -name *_R1.fastq.gz); do
    my_output_dir=${my_dir}/$(basename ${file/_R1.fastq.gz})
    mkdir -p ${my_output_dir}
    trim_galore \
    	-q 20 \
    	--fastqc \
    	--illumina \
    	--gzip \
    	--output_dir $my_output_dir \
    	--paired ${file} ${file/_R1.fastq.gz/_R2.fastq.gz} &
  done
  wait

}

03_kallisto(){
	local my_dir="${wd}/03-kallisto"

	# mkdir -p ${my_dir}/input/
	# for file in $(find ${wd}/02-merge_and_trim/ -name *fq.gz); do
	# 	new_name=$(basename $file | sed 's/_val_..fq.gz/.fastq.gz/')
	# 	ln -sf ${file} ${my_dir}/input/${new_name}
	# done

	# mkdir -p ${my_dir}/tmp
	# kallisto index -i ${my_dir}/tmp/gencode.v29.transcripts.kallisto_index \
	# 	/home/projects2/databases/gencode/release29/gencode.v29.transcripts.fa.gz

	mkdir -p ${my_dir}/output
	for fastq_file_r1 in $(find ${my_dir}/input -name "*_R1.fastq.gz"); do
		fastq_file_r2=$(echo $fastq_file_r1 | sed 's/_R1\./_R2./')	
		my_sample=$(basename $fastq_file_r1 | sed 's/_R..fastq.gz//')
		kallisto quant \
			-i ${my_dir}/tmp/gencode.v29.transcripts.kallisto_index \
			-t 10 \
			-b 100 \
			-o ${my_dir}/output/$my_sample/ $fastq_file_r1 $fastq_file_r2 &
	done
	wait
}

04_tximport (){
  s='04_tximport'

  my_dir="$wd/04-tximport"

  # # INPUT #####################################################################
  # mkdir -p $my_dir/input
  # ln -sf ${wd}/03-kallisto/output/* $my_dir/input

  # mkdir -p ${my_dir}/resources/gencode
  # ln -sf /home/projects2/databases/gencode/release29/gencode.v29.annotation.gtf \
  #  ${my_dir}/resources/gencode/

  # mkdir -p ${my_dir}/scripts
  # ln -sf ${wd}/../scripts/kallisto_tx2gene.R  ${my_dir}/scripts/

  # # # TPM #######################################################################

  # mkdir -p ${my_dir}/tmp/tx2gene
  # create_dict_table ${my_dir}/resources/gencode/gencode.v29.annotation.gtf > ${my_dir}/tmp/tx2gene/gencode_tsx2gene.tsv

  # # OUTPUT ####################################################################
  mkdir -p ${my_dir}/output/
  
  
  for file in $(find -L ${my_dir}/input/ -name abundance.h5); do
    abundance_file=$file
    sample_name=$(echo $file | sed -E 's|.*/([^/]*)/abundance.h5$|\1|')
    output="${my_dir}/output/${sample_name}_tx2gene_counts.tsv"
    
    Rscript ${my_dir}/scripts/kallisto_tx2gene.R \
      --tx2gene ${my_dir}/tmp/tx2gene/gencode_tsx2gene.tsv \
      --kallisto_abundance $abundance_file \
      --sample $sample_name \
      --output $output &

		# set +x; jobslimit 20; set -x
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

05_deseq2_pre(){
  s="DESEQ2"

  my_dir="${wd}/05-deseq2"
  mkdir -p ${my_dir}


  mkdir -p ${my_dir}/input/
  ln -sf ${wd}/04-tximport/output/* ${my_dir}/input/

  mkdir -p ${my_dir}/scripts/
  ln -s /home/users/rbarreiro/scripts/deseq2_batch/deseq2_wo_batch.R ${my_dir}/scripts

  mkdir -p ${my_dir}/tmp/combinations/

  #DONE MANUALLY
  #(printf "SAMPLE\tCONDITION\n"; for name in $(ls ../../../../input/); do X=${name/_tx2gene_counts.tsv};Y=$(echo $name  | cut -f3 -d '_'); printf "$X\t$Y\n"; done) > sample_info.tsv
  #for SF in EIF4A3 PHF5A SNRPB SNRPD2; do mkdir GSC_siCtl25vsGSC_si${SF}; done
  #for SF in RBMS1 SNRPG; do mkdir GSC_siCtl50vsGSC_si${SF}; done
  #for my_dir in $(find . -name "GSC*25*"); do echo ln -s /home/scratch60/rbarreiro_tfs_2020-06-28/04-deseq2/input/GSC_siCtl25_*_tx2gene_counts.tsv $my_dir/expression; done
  #for my_dir in $(find . -name "GSC*50*"); do echo ln -s /home/scratch60/rbarreiro_tfs_2020-06-28/04-deseq2/input/GSC_siCtl50_*_tx2gene_counts.tsv $my_dir/expression; done
  #for my_d in $(ls); do SF=$(basename $my_d | sed 's/.*GSC_//'); ln -s $(find /home/scratch60/rbarreiro_tfs_2020-06-28/04-deseq2/input -name "*${SF}*") $my_d/expression; done
  #for my_d in $(ls | grep -v .tsv); do head -1 sample_info.tsv > $my_d/sample_info/sample_info.tsv; for x in $(sed 's/vs/ /' <<< $my_d); do grep $x sample_info.tsv >> $my_d/sample_info/sample_info.tsv; done; done
}

05_deseq2_run(){
  my_dir="${wd}/05-deseq2"
  for combination in $(find ${my_dir}/tmp/combinations/ -mindepth 1 -maxdepth 1 | sort); do
    outpath=${my_dir}/output/$(basename $combination)
    mkdir -p $outpath
    group_a="C"
    group_b="KD"
    R4script ${my_dir}/scripts/deseq2_wo_batch.R \
      --group_a         $group_a \
      --group_b         $group_b \
      --sample_table    ${combination}/sample_info/sample_info.tsv \
      --col_categ       CONDITION \
      --input_extension _tx2gene_counts.tsv \
      --input_dir       ${combination}/expression  \
      --output          ${outpath} &

    set +x; jobslimit 36; set -x
  done
}

06_STAR(){
  s='STAR'
  my_dir="${wd}/06-STAR"
  mkdir -p ${my_dir}

  # mkdir -p ${my_dir}/input
  # ln -s  ${wd}/03-kallisto/input/* ${my_dir}/input/

  mkdir -p ${my_dir}/resources
  ln -sf /home/genomes/Homo_sapiens/hg38/hg38.25chrs.fa ${my_dir}/resources
  # cp -r $(realpath /home/genomes/Homo_sapiens/hg38/toSTAR.25chrs) ${my_dir}/resources

  mkdir -p ${my_dir}/resources/toSTAR.25chrs
  STAR \
    --runThreadN 32 \
    --runMode genomeGenerate \
    --genomeDir ${my_dir}/resources/toSTAR.25chrs \
    --genomeFastaFiles $(realpath /home/genomes/Homo_sapiens/hg38/hg38.25chrs.fa) \
    --sjdbGTFfile /home/projects2/databases/gencode/release29/gencode.v29.annotation.gtf \
    --sjdbOverhang 99 

  mkdir -p ${my_dir}/output

  my_genomeDir=$(realpath ${my_dir}/resources/toSTAR.25chrs/)
  STAR --genomeLoad LoadAndExit --genomeDir $my_genomeDir
  for file in $(find ${my_dir}/input -name "*_R1.fastq.gz"); do
    R1=$file
    R2=${file/_R1.fastq.gz/_R2.fastq.gz}
    mkdir -p ${my_dir}/output/$(basename ${file/_R1.fastq.gz/})
    cd ${my_dir}/output/$(basename ${file/_R1.fastq.gz/})
    STAR \
      --runThreadN            20 \
      --genomeDir             $my_genomeDir \
      --readFilesIn           $(realpath $R1) $(realpath $R2) \
      --readFilesCommand      zcat \
      --outSAMtype            BAM SortedByCoordinate \
      --genomeLoad            LoadAndKeep \
      --quantMode             GeneCounts \
      --outFileNamePrefix     ${my_dir}/output/$(basename ${file/_R1.fastq.gz/})/$(basename ${file/_R1.fastq.gz/})_ \
      --limitBAMsortRAM       15000000000

    cd ${my_dir}
  done
  STAR --genomeLoad Remove --genomeDir $my_genomeDir

  rm -rf ${my_dir}/resources # Index is 30GB+ better remove after done
}

07_rMATS(){
  s='RMATS'
  my_dir="${wd}/07-rMATS"

  # mkdir -p ${my_dir}/input
  # cp $(find ${wd}/06-STAR/output/ -name *.bam) ${my_dir}/input

  mkdir -p ${my_dir}/resources
  cp $(realpath /home/projects2/databases/gencode/release29/gencode.v29.annotation.gtf) ${my_dir}/resources/

  # find ${my_dir}/input -name "*_C*"  | tr "\n" "," > ${my_dir}/control.txt
  # find ${my_dir}/input -name "*_KD*" | tr "\n" "," > ${my_dir}/kd.txt

  read_len=151
  mkdir ${my_dir}/tmp/
  mkdir ${my_dir}/output/
  docker run -ti                                                    \
    -v     ${wd}:${wd}                                              \
    --user $(id -u):$(id -g)                                        \
    --rm   rmats-turbo-bioinfo                                      \
      --b1         ${my_dir}/control.txt                  \
      --b2         ${my_dir}/kd.txt                       \
      --gtf        ${my_dir}/resources/gencode.v29.annotation.gtf   \
      -t           paired                                           \
      --readLength ${read_len}                                      \
      --nthread    40                                               \
      --tstat      40                                               \
      --tmp        ${my_dir}/tmp/                                   \
      --od         ${my_dir}/output/
}

main
wait
