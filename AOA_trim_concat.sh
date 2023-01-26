#!/bin/bash

#Place AamoA.db.nuc.fasta in 01_databases/
#makeblastdb -dbtype prot -in 01_databases/AamoA.db.nuc.fasta -out AamoA.db.nuc

#Place amoA_AOA_cluster.fasta and amoA_AOA_cluster.tab in 01_databases/
#makeblastdb -dbtype nucl -in 01_databases/amoA_AOA_cluster.fasta -out 01_databases/amoA_AOA_cluster

# --clip_R1 21 --clip_R2 22 : for AOB amoA assembly
# --clip_R1 22 --clip_R2 21 : for AOA amoA gap 
# Adjust based on the reading frame check

basedir=$(pwd)
input=${1}
clip_R1=${2:-22}
clip_R2=${3:-21}
file=${input##*/}
sample=${file%%_1\.*}

mkdir -p 03_output/$sample
cd 03_output/$sample
mkdir -p info
mkdir -p fastq
mkdir -p cutadapt_before_DADA2
mkdir -p Intermediate_reads
mkdir -p out
mkdir -p cluster

echo $sample

cp ${basedir}/${input} fastq/
cp ${basedir}/${input/_1\./_2\.} fastq/

fastqc --noextract -o info fastq/*.gz
gunzip -cd fastq/${sample}_1.fastq.gz > fastq/${sample}_1.fastq
gunzip -cd fastq/${sample}_2.fastq.gz > fastq/${sample}_2.fastq


#Uncomment this if you want to check the clipping values for trim_galore based on your adapters.
#seqtk seq -a fastq/${sample}_1.fastq > info/${sample}_1.fasta
#seqtk seq -a fastq/${sample}_2.fastq > info/${sample}_2.1.fasta
#seqtk seq -r info/${sample}_2.1.fasta > info/${sample}_2.fasta
#rm info/${sample}_2.1.fasta
#seqtk sample info/${sample}_1.fasta 20 | transeq -sequence /dev/stdin -outseq info/${sample}_1_1.fasta -frame 1 -table 11 -clean
#seqtk sample info/${sample}_1.fasta 20 | transeq -sequence /dev/stdin -outseq info/${sample}_1_2.fasta -frame 2 -table 11 -clean
#seqtk sample info/${sample}_1.fasta 20 | transeq -sequence /dev/stdin -outseq info/${sample}_1_3.fasta -frame 3 -table 11 -clean
#rm info/${sample}_1.fasta
#seqtk sample info/${sample}_2.fasta 20 | transeq -sequence /dev/stdin -outseq info/${sample}_2_1.fasta -frame 1 -table 11 -clean
#seqtk sample info/${sample}_2.fasta 20 | transeq -sequence /dev/stdin -outseq info/${sample}_2_2.fasta -frame 2 -table 11 -clean
#seqtk sample info/${sample}_2.fasta 20 | transeq -sequence /dev/stdin -outseq info/${sample}_2_3.fasta -frame 3 -table 11 -clean
#rm info/${sample}_2.fasta
#echo "Checking reading frames of forward and reverse reads"
#for g in info/*.fasta; do
#    blastp -query $g -db ${basedir}/01_databases/AamoA.db -outfmt 6 -out ${g%.*}_raw.blastout
#done
#rm info/*fasta
#fframe=$(find info/${sample}_1_*_raw.blastout -printf '%s %p\n' | sort -nr | head -1 | sed 's/.*_\([123]\)_raw.blastout/\1/')
#rframe=$(find info/${sample}_2_*_raw.blastout -printf '%s %p\n' | sort -nr | head -1 | sed 's/.*_\([123]\)_raw.blastout/\1/')
#echo -e "Forward frame : ${fframe}\nReverse frame : ${rframe}"
#echo -e "Forward frame : ${fframe}\nReverse frame : ${rframe}" > info/${sample}_raw.frames
#rm info/*blastout

trim_galore -a TCGTGGGCAGCGTCAGATGTGT -a2 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -q 15 --paired --clip_R1 $clip_R1 --clip_R2 $clip_R2 fastq/${sample}_1.fastq fastq/${sample}_2.fastq 2>&1 | tee info/${sample}.trimming
rm fastq/*fastq

mv ${sample}_1_val_1.fq cutadapt_before_DADA2
mv ${sample}_2_val_2.fq cutadapt_before_DADA2
fastqc --noextract -o info cutadapt_before_DADA2/*

Rscript ${basedir}/AOA_trim_concat.R
fastqc --noextract -o info cutadapt_before_DADA2/filtered/*

for f in cutadapt_before_DADA2/*fq; do
    gzip $f
done

gunzip -c cutadapt_before_DADA2/filtered/${sample}_1_filt.fastq.gz > ${sample}_R1.fastq
gunzip -c cutadapt_before_DADA2/filtered/${sample}_2_filt.fastq.gz > ${sample}_R2.fastq
seqtk seq -a ${sample}_R1.fastq > ${sample}_R1_1.fasta
seqtk seq -a ${sample}_R2.fastq > ${sample}_R2_1.fasta
seqtk seq -r ${sample}_R2_1.fasta > ${sample}_R2_2.fasta
fasta2tab ${sample}_R1_1.fasta > ${sample}_R1_1.tab
fasta2tab ${sample}_R2_2.fasta > ${sample}_R2_2.tab
awk -F "\t" '{print $1}' ${sample}_R1_1.tab > ${sample}_R1_2.tab
awk 'BEGIN {FS="/t" } NR { $1=NR} { print}' ${sample}_R1_2.tab > ${sample}_R1_3.tab
sed 's/^/AOA-'${sample}'-/' ${sample}_R1_3.tab > ${sample}_R1_4.tab
awk -F "\t" '{print $2}' ${sample}_R1_1.tab > ${sample}_R1_5.tab
awk -F "\t" '{print $2}' ${sample}_R2_2.tab > ${sample}_R2_3.tab
paste ${sample}_R1_5.tab ${sample}_R2_3.tab > ${sample}_1.tab
sed 's/\t//g' ${sample}_1.tab > ${sample}_2.tab
paste ${sample}_R1_4.tab ${sample}_2.tab > ${sample}_reads.tab
tab2fasta ${sample}_reads.tab > ${sample}_reads.fasta
usearch --cluster_fast ${sample}_reads.fasta -id 1 -centroids ${sample}_reads_id_1.fasta -sizeout 2>&1 | tee info/${sample}.clustering
transeq -sequence ${sample}_reads_id_1.fasta -outseq ${sample}_reads_id_1_1.fasta -frame 3 -table 11 -clean
grep -v X -B 1 ${sample}_reads_id_1_1.fasta > ${sample}_reads_id_1_2.fasta
usearch -fastx_truncate ${sample}_reads_id_1_2.fasta -trunclen 133 -fastaout ${sample}_reads_id_1_3.fasta 2>&1 | tee info/${sample}.truncated
fasta2tab ${sample}_reads_id_1_3.fasta > ${sample}_reads_id_1_4.tab
awk -F '\t' '{print $1}' ${sample}_reads_id_1_4.tab > ${sample}_reads_id_1_5.tab
sed 's/_3//' ${sample}_reads_id_1_5.tab > ${sample}_reads_id_1_6.tab
fasta2tab ${sample}_reads_id_1.fasta > ${sample}_reads_id_1.tab
awk -F '\t' '{print $2, $1}' OFS='\t' ${sample}_reads_id_1.tab > ${sample}_reads_id_1_2.tab
awk 'NR==FNR {h[$2] = $1; next} {print $1, h[$1]}' OFS='\t' ${sample}_reads_id_1_2.tab ${sample}_reads_id_1_6.tab > ${sample}_trans_reads_id_1.tab
tab2fasta ${sample}_trans_reads_id_1.tab > ${sample}_trans_reads_id_1.fasta
usearch -unoise3 ${sample}_trans_reads_id_1.fasta -zotus ${sample}_zotus.fasta -tabbedout ${sample}_unoise3.txt -minsize 2 2>&1 | tee info/${sample}.denoising
fasta2tab ${sample}_zotus.fasta > ${sample}_zotus.tab
sed 's/Zotu/amp/' ${sample}_zotus.tab > ${sample}_zotus_1.tab
awk -F '\t' '{print $1, $3}' OFS='\t' ${sample}_unoise3.txt > ${sample}_unoise3_1.tab
awk 'NR==FNR {h[$2] = $1; next} {print $1, h[$1]}' OFS='\t' ${sample}_unoise3_1.tab ${sample}_zotus_1.tab > ${sample}_list_non_singletons.tab
awk -F '\t' '{print $2}' ${sample}_list_non_singletons.tab > ${sample}_list_non_singletons_1.tab
fasta2tab ${sample}_trans_reads_id_1.fasta > ${sample}_valids_reads_id_1.tab
awk -F '\t' '{print $2, $1}' OFS='\t' ${sample}_valids_reads_id_1.tab > ${sample}_valids_reads_id_1_1.tab
awk 'NR==FNR {h[$2] = $1; next} {print $1, h[$1]}' OFS='\t' ${sample}_valids_reads_id_1_1.tab ${sample}_list_non_singletons_1.tab > ${sample}_Reads_non_singletons.tab
tab2fasta ${sample}_Reads_non_singletons.tab > ${sample}_Reads_non_singletons_gap.fasta
rm ${sample}_R1_*.*
rm ${sample}_R2_*.*
rm ${sample}_reads.fasta  
rm ${sample}_R*.fastq
rm ${sample}_reads_id_1_*.fasta
rm ${sample}_zotus.fasta
mv ${sample}_reads_id_1.fasta Intermediate_reads
mv ${sample}_trans_reads_id_1.fasta Intermediate_reads
mv ${sample}_Reads_non_singletons_gap.fasta out
rm *.tab
rm *.txt

echo "Assigning reads to taxonomic clusters"
blastn -max_target_seqs 1 -query out/${sample}_Reads_non_singletons_gap.fasta -db ${basedir}/01_databases/amoA_AOA_cluster -outfmt 6 > Blast_AOA_gap_cluster_DB_${sample}.tab
awk -F '[\t]' '{print $2}' Blast_AOA_gap_cluster_DB_${sample}.tab > Blast_AOA_gap_cluster_DB_${sample}_1.tab
awk -F '[_]' '{print $1}' Blast_AOA_gap_cluster_DB_${sample}_1.tab > Blast_AOA_gap_cluster_DB_${sample}_2.tab
awk -F '[;]' '{print $2}' Blast_AOA_gap_cluster_DB_${sample}.tab > Blast_AOA_gap_cluster_DB_${sample}_3.tab
awk -F '[=]' '{print $2}' Blast_AOA_gap_cluster_DB_${sample}_3.tab > Blast_AOA_gap_cluster_DB_${sample}_4.tab
paste Blast_AOA_gap_cluster_DB_${sample}_2.tab Blast_AOA_gap_cluster_DB_${sample}_4.tab > Blast_AOA_gap_cluster_DB_${sample}_5.tab
awk '{sums[$1] += $2} END { for (i in sums) printf("%s %s\n", i, sums[i])}' Blast_AOA_gap_cluster_DB_${sample}_5.tab > Blast_AOA_gap_cluster_DB_${sample}_6.tab
awk -F " " '{print $2, $1}' OFS='\t' Blast_AOA_gap_cluster_DB_${sample}_6.tab > Blast_AOA_gap_cluster_DB_${sample}_7.tab
awk 'NR==FNR {h[$2] = $1; next} {print $1, h[$1]}' OFS='\t' Blast_AOA_gap_cluster_DB_${sample}_7.tab ${basedir}/01_databases/amoA_AOA_cluster.tab > Blast_AOA_gap_cluster_DB_${sample}_8.tab
sed '1i Match_base_id\tAbundances-AOA-'${sample}'' Blast_AOA_gap_cluster_DB_${sample}_8.tab > Blast_AOA_gap_cluster_DB_${sample}_9.tab
sed '1i Match_base_id' ${basedir}/01_databases/amoA_AOA_cluster.tab > List_AOA_Cluster_with_header.tab
awk -F "\t" '{print $2}' Blast_AOA_gap_cluster_DB_${sample}_9.tab > Calcul_Blast_AOA_gap_cluster_DB_${sample}.tab
paste List_AOA_Cluster_with_header.tab Calcul_Blast_AOA_gap_cluster_DB_*.tab > AOA_map_to_cluster_name.tab
rm Blast_AOA_gap_cluster_DB_${sample}_*.tab
rm Calcul_Blast_AOA_gap_cluster_DB_${sample}.tab
rm List_AOA_Cluster_with_header.tab
mv Blast_AOA_gap_cluster_DB_${sample}.tab out/${sample}_amoA_AOA_cluster.blastout
mv AOA_map_to_cluster_name.tab out/${sample}_amoA_AOA_map_to_cluster.tab
