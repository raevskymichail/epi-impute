#!/bin/bash

## sudo apt install moreutils # for sponge function

# Download SRA accessions and annotation before run

WORKDIR="/home/mraevsky/MIPT/epi-impute/data/raw/"

cd "$WORKDIR"
accession_list="SRR_Acc_List.txt"
SraRunTable="SraRunTable.txt"

# Download hg38 blacklist
curl -o ../annotations/hg38/hg38_encode_blacklist.bed https://www.encodeproject.org/files/ENCFF419RSJ/@@download/ENCFF419RSJ.bed.gz
blackListFile="../annotations/hg38/hg38_encode_blacklist.bed"

promoters_enhancers="../annotations/hg38/fantom5_hg38_hgnc_promoters_enhancers.bed"

# Obtain prefetch files
# screen -S sra_toolkit_run
# ./sratoolkit.2.9.6-1-ubuntu64/bin/prefetch --option-file "$accession_list"

# Obtain fastq files
mkdir ./sra_prefetch/sra_fastq
fastq_path="./sra_prefetch/sra_fastq/"
fastq_dump="./sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump"

while IFS= read -r sranumber
do
	# faster_dump -O "$fastq_path" 
	"$fastq_dump" --outdir "$fastq_path" -I --split-3 --gzip --skip-technical \
		--read-filter pass --dumpbase --clip "$sranumber"

	mv "${fastq_path}${sranumber}_pass_1.fastq.gz" "${fastq_path}${sranumber}_R1.fq.gz"
	mv "${fastq_path}${sranumber}_pass_2.fastq.gz" "${fastq_path}${sranumber}_R2.fq.gz"

	echo "$sranumber downloaded!"
done < "$accession_list"


### Remove singletons reads
## Before assembly reads will often be grouped into "pairs" (ones that have mates) 
## and "singletons" (reads without a mate; 
## generally because the potential mate is of poor quality)
# while IFS= read -r sranumber
# do
# 	pigz -d "${fastq_path}${sranumber}_R1.fastq.gz"
# 	pigz -d "${fastq_path}${sranumber}_R2.fastq.gz"
# 	fastq_pair "${fastq_path}${sranumber}_R1.fastq" "${fastq_path}${sranumber}_R2.fastq" &
# 	pigz "${fastq_path}${sranumber}_R1.fq"
# 	pigz "${fastq_path}${sranumber}_R2.fq"
# 	echo "Output files generated: ${sranumber}_R1.paired.fq.gz and ${sranumber}_R2.paired.fq.gz !"
# done < "$accession_list"


# QC for a single paired-end read
example_fastq="SRR5353377_R1.fq.gz"
fastqc "./sra_prefetch/sra_fastq/$example_fastq" --outdir="$WORKDIR"

# Trim adaptors
mkdir ./sra_prefetch/sra_fastq_trim
fastq_trim_path="./sra_prefetch/sra_fastq_trim/"

trimmomatic="./Trimmomatic-0.39/trimmomatic-0.39.jar"
adapters_path="./Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:1:true"

while IFS= read -r sranumber
do
	## NGmerge trimming
	# ./NGmerge/NGmerge -n 30 -a  -1 "${fastq_path}${sranumber}_R1.fastq.gz"  -2 "${fastq_path}${sranumber}_R2.fastq.gz" -o "${fastq_trim_path}${sranumber}" -z -v
	# replace ".fastq.gz" with ".fq.gz"
	## Trimmomatic
	java -jar "$trimmomatic" PE -phred33 \
		"${fastq_path}${sranumber}_R1.fq.gz" \
		"${fastq_path}${sranumber}_R2.fq.gz" \
		"${fastq_trim_path}${sranumber}_R1.paired.fq.gz" \
		"${fastq_trim_path}${sranumber}_R1.unpaired.fq.gz" \
		"${fastq_trim_path}${sranumber}_R2.paired.fq.gz" \
		"${fastq_trim_path}${sranumber}_R2.unpaired.fq.gz" \
		ILLUMINACLIP:"$adapters_path" \
		TRAILING:3 \
		SLIDINGWINDOW:4:15 \
		MINLEN:25 
done < "$accession_list"


# Download and index hg38 genome
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mkdir ./genomes/ && cd ./genomes/
mkdir ./hg38 && cd ./hg38
curl -o hg38_genome.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip hg38_genome.fna.gz
bowtie2-build hg38_genome.fna hg38_index
cd ../../

# Align reads to hg38
n_threads=30 # number of threads
ref_genome_indexed="./genomes/hg38/hg38_index"

mkdir ./sra_prefetch/sra_sam
sra_sam_path="./sra_prefetch/sra_sam/"

mkdir ./sra_prefetch/sra_bw
mkdir ./sra_prefetch/sra_bw/by_sra
mkdir ./sra_prefetch/sra_bw/by_celltype
sra_bw_by_sra_path="./sra_prefetch/sra_bw/by_sra/"
sra_bw_by_celltype_path="./sra_prefetch/sra_bw/by_celltype/"

mkdir ./sra_prefetch/sra_bam
sra_bam_path="./sra_prefetch/sra_bam/"
mkdir ./sra_prefetch/sra_bam_merged
sra_bam_merged_path="./sra_prefetch/sra_bam_merged/"


mkdir ./sra_prefetch/bowtie2_log
bowtie2_log_path="./sra_prefetch/bowtie2_log/"

picard_jar="./picard/build/libs/picard.jar"
hg38_effective_size=2913022398
binSize=10

while IFS= read -r sranumber
do
	# better alignment results are frequently achieved with --very-sensitive
	# use -X 2000 to allow larger fragment size (default is 500)
	bowtie2 --threads "$n_threads" --very-sensitive -X 2000 --dovetail -x "$ref_genome_indexed" \
		-1 "${fastq_trim_path}${sranumber}_R1.paired.fq.gz" \
		-2 "${fastq_trim_path}${sranumber}_R2.paired.fq.gz" \
	  	-S "${sra_sam_path}${sranumber}.sam" 2> "${bowtie2_log_path}${sranumber}.bowtie2.log"
	# samtools sort -@ "$n_threads" -O bam -o "${sra_bam_path}${sranumber}.sorted.bam"
	# samtools index -@ "$n_threads" "${sra_bam_path}${sranumber}.sorted.bam"
done < "$accession_list"

while IFS= read -r sranumber
do
	# Retain only the properly paired reads and reads with quality above 30
	samtools view -SbhF 4 -f2 -q30 "${sra_sam_path}${sranumber}.sam" > "${sra_bam_path}${sranumber}.bam"
	# Filter Blacklist regions
	intersectBed -v -abam "${sra_bam_path}${sranumber}.bam" \
		 -b "${blackListFile}" | samtools sort -@ "$n_threads" > "${sra_bam_path}${sranumber}_sorted.bam"
	## Remove mitochondrial reads
	# 	samtools view -@ $n_threads -h ${sample}.bam | grep -v chrM | samtools sort -@ $n_threads -O bam -o ${sample}.bam
	## Remove PCR duplicates and mitochondrial reads
	java -Xmx2g -jar "$picard_jar" MarkDuplicates INPUT="${sra_bam_path}${sranumber}_sorted.bam" \
		OUTPUT="${sra_bam_path}${sranumber}_nodup.bam" \
		METRICS_FILE="${sra_bam_path}${sranumber}_nodup_stats" \
		REMOVE_DUPLICATES=True
	    
	samtools sort -@ "$n_threads" "${sra_bam_path}${sranumber}_nodup.bam" -o "${sra_bam_path}${sranumber}_nodup_sorted.bam" 
	samtools index -@ "$n_threads" "${sra_bam_path}${sranumber}_nodup_sorted.bam"
	samtools view -b "${sra_bam_path}${sranumber}_nodup_sorted.bam" \
		chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
		chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX > "${sra_bam_path}${sranumber}_cleaned.bam"

	rm "${sra_bam_path}${sranumber}_nodup.bam"

	samtools index -@ "$n_threads" "${sra_bam_path}${sranumber}_cleaned.bam"
	samtools idxstats "${sra_bam_path}${sranumber}_cleaned.bam" > "${sra_bam_path}${sranumber}_cleaned_chrom_stat.txt"

done < "$accession_list"

# Fix Tn5 sequence preference bias
# Don't need since ATACorrect only works for cutsites correction

while IFS= read -r sranumber
do
	## Shifting reads
	# In the first ATAC-seq paper (Buenrostro et al., 2013), 
	# all reads aligning to the + strand were offset by +4 bp, 
	# and all reads aligning to the – strand were offset −5 bp, 
	# since Tn5 transposase has been shown to bind as a dimer and insert two adaptors 
	# separated by 9 bp (Adey et al., 2010).
	# use --ATACshift
	# NOTE: if the --shift or --ATACshift option are used, only propertly paired-end reads will be used
	alignmentSieve --numberOfProcessors "max" --ATACshift \
		--bam "${sra_bam_path}${sranumber}_cleaned.bam" \
		-o "${sra_bam_path}${sranumber}_cleaned.tmp.bam"

	# the bam file needs to be sorted again
	samtools sort -@ "$n_threads" "${sra_bam_path}${sranumber}_cleaned.tmp.bam" -O bam -o "${sra_bam_path}${sranumber}_cleaned.shifted.bam"
	samtools index -@ "$n_threads" "${sra_bam_path}${sranumber}_cleaned.shifted.bam"
	rm "${sra_bam_path}${sranumber}_cleaned.tmp.bam"

done < "$accession_list"


## TEST



## count bams in promoter and enhancer regions

## SEQUENTIAL

# awk -F "," 'NR>1 {print "./sra_prefetch/sra_bam/"$1"_cleaned.bam"}' "$SraRunTable" > bams_list.txt

# # count each bam independently because of bedtools multicov limit on 1020 bams files

# cp "$promoters_enhancers" promoters_enhancers_counts.bed # create an empty file

# while IFS= read -r bam_file
# do
# 	# count reads in promoters and enhancers
# 	bedtools multicov -bams $bam_file -bed "$promoters_enhancers" > i_counts.bed
# 	awk -F "\t" '{print $5}' i_counts.bed > i_counts_col
# 	cp promoters_enhancers_counts.bed counts_in_regions.bed
# 	paste counts_in_regions.bed i_counts_col | column -s $'\t' -t > promoters_enhancers_counts.bed
# 	rm counts_in_regions.bed
# 	rm i_counts.bed
# 	rm i_counts_col
# done < bams_list.txt

# # Add header for the count matrix by adding BED colnames 
# # and converting accesion_list into tab delim variable
# # For count matrix convert SRA -> "<celltype>_SRA" (ex. "HSC_SRR5353377")
# cell_ids=$(awk -F "," 'BEGIN{OFS="\t"} {print $6"_"$1}' "$SraRunTable" | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/\t/g' )
# # # replace all blanks
# # bar=${foo// /.}
# # add header to the count matrix
# cp promoters_enhancers_counts.bed promoters_enhancers_counts.tmp.bed
# echo -e "chr\tstart\tend\tregion_name\t${cell_ids}" | cat - promoters_enhancers_counts.tmp.bed > promoters_enhancers_counts.bed
# rm promoters_enhancers_counts.tmp.bed

## PARALLEL

# test_arr=(SRR5353377 SRR5353377)
# parallel count_in_regions ::: SRR5353377 SRR5353378
# parallel count_in_regions ::: "${test_arr[@]}"

mkdir ./sra_prefetch/sra_bed
sra_bed_path="./sra_prefetch/sra_bed/"

export SraRunTable_sra=$(awk -F "," 'NR>1 {print $1}' "$SraRunTable")

promoters_enhancers_500_200_win='../annotations/hg38/fantom5_hg38_hgnc_promoters_enhancers_500_200_win.bed'

function count_in_regions {
   local sra=$1
   local regions=$2
   local bam_in=$(echo "./sra_prefetch/sra_bam/"${sra}"_cleaned.bam")
   local bed_tmp=$(echo "./sra_prefetch/sra_bed/"${sra}"_cleaned.tmp.bed")
   local bed_out=$(echo "./sra_prefetch/sra_bed/"${sra}"_cleaned.bed")
   # local bed_regions
   echo $bam_in
   bedtools multicov -bams $bam_in -bed $regions > "$bed_tmp"
   sed 's/ \+/\t/g' "$bed_tmp" > "$bed_out"
   rm "$bed_tmp"
}

export -f count_in_regions

parallel count_in_regions ::: $(echo "$SraRunTable_sra") ::: $(echo "$promoters_enhancers_500_200_win")

# parallel count_in_regions ::: $(cat SraRunTable_sra.txt)
# parallel count_in_regions :::: <(cat SraRunTable_sra.txt) 
# cat SraRunTable_sra.txt | parallel --xargs count_in_regions
# echo "$SraRunTable_sra" | parallel 'count_in_regions {}'

## ADD HEADERS AND MERGE
awk -F "," 'NR>1 {print "./sra_prefetch/sra_bed/"$1"_cleaned.bed"}' "$SraRunTable" > beds_list.txt
cp "$promoters_enhancers_500_200_win" promoters_enhancers_500_200_win_counts.bed # create an empty file

while IFS= read -r bed_file
do
	awk -F "\t" '{print $5}' "$bed_file" > counts_col
	paste -d'\t' promoters_enhancers_500_200_win_counts.bed counts_col | sponge promoters_enhancers_500_200_win_counts.bed
	rm counts_col
	echo "$bed_file"
done < beds_list.txt

# Add header for the count matrix by adding BED colnames 
# and converting accesion_list into tab delim variable
# For count matrix convert SRA -> "<celltype>_SRA" (ex. "HSC_SRR5353377")
cell_ids=$(awk -F "," 'BEGIN{OFS="\t"} {print $6"_"$1}' "$SraRunTable" | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/\t/g' )
# # replace all blanks
# bar=${foo// /.}
# add header to the count matrix
cp promoters_enhancers_500_200_win_counts.bed counts.tmp
echo -e "chr\tstart\tend\tregion_name\t${cell_ids}" | cat - counts.tmp > promoters_enhancers_500_200_win_counts.bed
rm counts.tmp



## bam to bigwig
# optionally normalize using 1x effective genome size
# effective genome size: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
# bamCoverage --numberOfProcessors "$n_threads" --binSize 10 --normalizeUsing RPGC \
#   --effectiveGenomeSize "$hg38_effective_size" --bam "${sra_bam_path}${sranumber}_cleaned.bam" -o "${sranumber}_cleaned.bw"

# for single SRA accession

while IFS= read -r sranumber
do
	bamCoverage --numberOfProcessors "max" --binSize "$binSize" \
	  --bam "${sra_bam_path}${sranumber}_cleaned.shifted.bam" -o "${sra_bw_by_sra_path}${sranumber}.bw"
done < "$accession_list"

# for each cell type

awk -F "," 'NR>1 {print $6}' "$SraRunTable" | sort | uniq > celltypes.txt # print cell type

while IFS= read -r celltype_i
do
	# find SRAs of given celltype
	export celltype_i
	# awk -F "," 'NR>1 {if ($6 == "$celltype_i") print $1}' "$SraRunTable" | head
	awk -F "," 'NR>1 {if (match($6, ENVIRON["celltype_i"])) print $1}' "$SraRunTable" > sra_celltype_i.txt
	# convert SRAs to samtools bam argument line
	arg_bams_celltype_i=$(awk -F "," '{print "./sra_prefetch/sra_bam/"$1"_cleaned.shifted.bam"}' sra_celltype_i.txt | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' )
	# merge BAMs by celltype
	samtools merge -@ "$n_threads" "${sra_bam_merged_path}${celltype_i}.merged.bam" $arg_bams_celltype_i
	samtools index -@ "$n_threads" "${sra_bam_merged_path}${celltype_i}.merged.bam"
	# convert to bigwig
	bamCoverage --numberOfProcessors "max" --binSize "$binSize" \
	  --bam "${sra_bam_merged_path}${celltype_i}.merged.bam" -o "${sra_bw_by_celltype_path}${celltype_i}.bw"

	rm sra_celltype_i.txt

done < celltypes.txt

rm celltypes.txt



# MACS2 Peaks?



# # mkdir ./sra_prefetch/sra_bam
# for file in ./sra_prefetch/sra_sam/*.sam
# do
# 	samfile=$(basename -- "$file")
# 	sranumber="${samfile%.*}"
# 	samtools view -bS "./sra_prefetch/sra_sam/$samfile" > "./sra_prefetch/sra_bam/$sranumber.bam"
# 	echo "$sranumber.sam to bam finished!"
# done



# ## Utils

# # Concat promoter and enhancer lists into one
# # Take only chr, start, end, name columns. Ignore strand and others
# awk -F "\t" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' ../annotations/hg38/fantom5_hg38_hgnc_gene_promoters.bed > ../annotations/hg38/fantom5_hg38_hgnc_promoters_enhancers.bed
# awk -F "\t" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' ../annotations/hg38/fantom5_hg38_enhancer_tss_associations.bed >> ../annotations/hg38/fantom5_hg38_hgnc_promoters_enhancers.bed

# # To convert to hg19 -> hg38 coordinated in BED
# # pip install CrossMap
# # no stand in enhancer-TSS were -> "+"
# python CrossMap.py bed GRCh37_to_GRCh38.chain.gz \
# 	./hg19/fantom5_hg19_enhancer_tss_associations.all_plus.bed \
# 	./hg38/fantom5_hg38_enhancer_tss_associations.bed

# python CrossMap.py bed GRCh37_to_GRCh38.chain.gz \
# 	./hg19/fantom5_hg19_hgnc_gene_main_promoters_500_200_win.bed \
# 	./hg38/fantom5_hg38_hgnc_gene_main_promoters_500_200_win.bed

# python CrossMap.py bed GRCh37_to_GRCh38.chain.gz \
# 	./hg19/fantom5_hg19_hgnc_gene_promoters.bed \
# 	./hg38/fantom5_hg38_hgnc_gene_promoters.bed

# # Print how many enhancers / genes were mapped or unmapped

# awk -F ";" '{print $3}' fantom5_hg38_enhancer_tss_associations.bed | sort | uniq | wc -l
# awk -F ";" '{print $3}' fantom5_hg38_enhancer_tss_associations.bed.unmap | sort | uniq | wc -l

# wc -l fantom5_hg38_hgnc_gene_main_promoters_500_200_win.bed
# wc -l fantom5_hg38_hgnc_gene_main_promoters_500_200_win.bed.unmap
