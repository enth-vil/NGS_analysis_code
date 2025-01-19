#!/bin/bash


Tumor_Raw_R1="/home/kayal/Desktop/test_pupil/raw_reads/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz"
echo "The file is located at: $Tumor_Raw_R1"

Tumor_Raw_R2="/home/kayal/Desktop/test_pupil/raw_reads/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz"
echo "The file is located at: $Tumor_Raw_R2"






cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 150 -M 152 -j 8 -o /home/kayal/Desktop/test_pupil/trimmed_reads/Tumor_R1.fastq -p /home/kayal/Desktop/test_pupil/trimmed_reads/Tumor_R2.fastq $Tumor_Raw_R1 $Tumor_Raw_R2

Tumor_fil_R1="/home/kayal/Desktop/test_pupil/trimmed_reads/Tumor_R1.fastq"

Tumor_fil_R2="/home/kayal/Desktop/test_pupil/trimmed_reads/Tumor_R2.fastq"

Normal_Raw_R1="/home/kayal/Desktop/test_pupil/raw_reads/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz"
echo "The file is located at: $Normal_Raw_R1"

Normal_Raw_R2="/home/kayal/Desktop/test_pupil/raw_reads/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz"
echo "The file is located at: $Normal_Raw_R2"


cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 150 -M 152 -j 8 -o /home/kayal/Desktop/test_pupil/trimmed_reads/Normal_R1.fastq -p /home/kayal/Desktop/test_pupil/trimmed_reads/Normal_R2.fastq $Normal_Raw_R1 $Normal_Raw_R2

Normal_fil_R1="/home/kayal/Desktop/test_pupil/trimmed_reads/Normal_R1.fastq"

Normal_fil_R2="/home/kayal/Desktop/test_pupil/trimmed_reads/Normal_R2.fastq"


ref="/home/kayal/Desktop/test_pupil/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

bowtie2-build --threads 8 $ref hg38

bowtie2 --threads 8 --rg-id tumor_BWB90210-2332_lane1 --rg SM:tumor --rg LB:lib1 --rg PL:ILLUMINA --rg PU:FS10003130_11_BWB90210-2332_1_tumor -x /home/kayal/Desktop/test_pupil/hg38 -1 $Tumor_fil_R1 -2 $Tumor_fil_R1 -S /home/kayal/Desktop/test_pupil/bam/Tum_var.sam

Tum_sam="/home/kayal/Desktop/test_pupil/bam/Tum_var.sam"

bowtie2 --threads 4 --rg-id normal_BWB90210-2332_lane1 --rg SM:normal --rg LB:lib2 --rg PL:ILLUMINA --rg PU:FS10003130_11_BWB90210-2332_1_normal -x /home/kayal/Desktop/test_pupil/hg38 -1 $Normal_fil_R1 -2 $Normal_fil_R2 -S /home/kayal/Desktop/test_pupil/bam/Norm_var.sam

Norm_sam="/home/kayal/Desktop/test_pupil/bam/Norm_var.sam"

samtools view -b -S $Tum_sam -o /home/kayal/Desktop/test_pupil/bam/Tum_var.bam

samtools view -b -S $Norm_sam -o /home/kayal/Desktop/test_pupil/bam/Norm_var.bam

Tum_bam="/home/kayal/Desktop/test_pupil/bam/Tum_var.bam"

Norm_bam="/home/kayal/Desktop/test_pupil/bam/Norm_var.bam"


samtools sort $Tum_bam -o /home/kayal/Desktop/test_pupil/bam/Tum_var.sorted.bam

samtools sort $Norm_bam -o /home/kayal/Desktop/test_pupil/bam/Norm_var.sorted.bam

Tum_bam_sort="/home/kayal/Desktop/test_pupil/bam/Tum_var.sorted.bam"

Norm_bam_sort="/home/kayal/Desktop/test_pupil/bam/Norm_var.sorted.bam"

samtools index $Norm_bam_sort

samtools index $Tum_bam_sort 

mv Tum_var.sorted.bam.bai /home/kayal/Desktop/test_pupil/bam

mv Norm_var.sorted.bam.bai /home/kayal/Desktop/test_pupil/bam

java -jar /home/kayal/mywork/picard/picard.jar MarkDuplicates I=$Tum_bam_sort O=/home/kayal/Desktop/test_pupil/dedup_bam/Tum_var_markdup.bam M=/home/kayal/Desktop/test_pupil/dedup_bam/marked_dup_metrics_Norm_var.txt REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true

java -jar /home/kayal/mywork/picard/picard.jar MarkDuplicates I=$Norm_bam_sort O=/home/kayal/Desktop/test_pupil/dedup_bam/Norm_var_markdup.bam M=/home/kayal/Desktop/test_pupil/dedup_bam/marked_dup_metrics_Norm_var.txt REMOVE_DUPLICATES=true

Tum_markdup="/home/kayal/Desktop/test_pupil/dedup_bam/Tum_var_markdup.bam"

Norm_markdup="/home/kayal/Desktop/test_pupil/dedup_bam/Norm_var_markdup.bam"

samtools index $Norm_markdup

samtools index $Tum_markdup

mv Tum_var_markdup.sorted.bam.bai /home/kayal/Desktop/test_pupil/dedup_bam

mv Norm_var_markdup.sorted.bam.bai /home/kayal/Desktop/test_pupil/dedup_bam

gatk CreateSequenceDictionary -R $ref


##To create germline variant vcf file(normal samples only)

gatk HaplotypeCaller  -R $ref -I /home/kayal/Desktop/test_pupil/dedup_bam/Norm_var_markdup.bam -O Normal.vcf.gz

gatk VariantFiltration -R $ref -V Normal.vcf.gz -O Germline.vcf.gz --filter-name "LowQual" --filter-expression "QUAL < 30.0 || DP < 10"

gatk VariantAnnotator -R $ref -V Germline.vcf.gz -O germline_resource.vcf.gz --annotation AlleleFraction

gatk IndexFeatureFile -I germline_resource.vcf.gz

##bed file was created using the given csv file (PANCAN_PDAC_100plex_ref.csv) in r and was used to filter variants

gatk Mutect2 -R $ref -I $Tum_markdup -I $Norm_markdup -tumor tumor -normal normal -L fil.bed --germline-resource germline_resource.vcf.gz -O somatic.vcf.gz

bcftools view -f PASS -R fil.bed germline_resource.vcf.gz -o filtered_background_mut.vcf

##Total number of bases was calculated to calculate background mutation rate

cat fil.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2-1 }END{print SUM}'

bcftools stats somatic.vcf.gz

##Mean background mutation level

bcftools query -f '%AF\n' filtered_background_mut.vcf > af_values.txt

awk '{sum += $1; count++} END {if (count > 0) print sum / count; else print "No values"}' af_values.txt









