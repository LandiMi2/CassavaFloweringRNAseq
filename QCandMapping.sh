#runing fastqc
fastqc -o qc -t 10 *

#trim using trim_galore - version 0.6.10
#for T0
data="/data01/dataRepository/BioinfoData/CassavaFlowering/data/T0"
for i in "$data"/*_1.fq.gz;
do
   prefix=$(basename $i _1.fq.gz)
   trim_galore --paired --clip_R1 20 --clip_R2 20 --three_prime_clip_R1 2 --three_prime_clip_R2 2 \
   --fastqc --quality 20 -j 8 --illumina --length 100 --trim-n --polyA \
   "$data"/"${prefix}"_1.fq.gz "$data"/"${prefix}"_2.fq.gz

done

#this was done separately for T1 and T2

## mapping to reference
#  for example DSC120 samples 
data="/data01/mlandi/CassavaFlowering/dataFiltered/T0/Nolight/DSC120"
genome="/data01/mlandi/CassavaFlowering/genome"
out="/data01/mlandi/CassavaFlowering/star/T0/Nolight/DSC120"

#run star
for r1 in "$data"/*_1_val_1.fq.gz; do
    base=$(basename "$r1" _1_val_1.fq.gz)
    r2="$data/${base}_2_val_2.fq.gz"

    STAR --runThreadN 8 --genomeDir "$genome" --readFilesIn "$r1" "$r2" --readFilesCommand zcat \
    --outFilterMismatchNoverReadLmax 0.06 --outFileNamePrefix "$out/${base}_" --outSAMtype BAM SortedByCoordinate
done

## gff file downloaded from phytozome database: https://data.jgi.doe.gov/refine-download/phytozome?organism=Mesculenta&phytozome_version=14&expanded=Phytozome-671
#convert gtf to bed
convert2bed -i gtf < Mesculenta_671_v8.1.gene_exons.gtf > Mesculenta_671_v8.1.gene_exons.bed

## about strandedness, we confirm using RSeqC

#RSeqC
infer_experiment.py -r genome/Mesculenta_671_v8.1.gene_exons.bed -i star/T0/light/DSC120/T0_1_Aligned.sortedByCoord.out.bam

#below is the output
This is PairEnd Data
Fraction of reads failed to determine: 0.0547
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4739
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4714

#counts  - featureCounts v1.6.0
featureCounts -T 10 -p -B -t exon -g gene_id -a ../genome/Mesculenta_671_v8.1.gene_exons.gtf -o counts.txt bams/*.bam


