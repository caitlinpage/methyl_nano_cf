# code for processing wgbs samples into bismark
# basically all the slurm steps

# initially had tried using bwa-meth, but the output of that didn't like working with bismark (read level meth)
# so switched to bismark

# start with fastq file

# genome prep (builds index)
module load bismark

bismark_genome_preparation /researchers/caitlin.page

# align
# run from folder fastq is in
module load bismark
module load samtools
module load bowtie2

bismark --parallel 4 -o /researchers/caitlin.page/pancreatic_cancer_wgbs/bismark_res /
--genome /researchers/caitlin.page -q /
-1 SRR1646792_GSM1541788_38-Lg_rep1_BS-Seq_Homo_sapiens_Bisulfite-Seq_1.fastq.gz /
-2 SRR1646792_GSM1541788_38-Lg_rep1_BS-Seq_Homo_sapiens_Bisulfite-Seq_2.fastq.gz
# rrbs data just the one fastq file (single end)

# methylation extractor (did cytosine report but don't need that)
# did from bismark_res folder
module load bismark
module load samtools
module load bowtie

bismark_methylation_extractor --multicore 8 --cytosine_report /
--genome_folder /researchers/caitlin.page -p /
SRR1646792_GSM1541788_38-Lg_rep1_BS-Seq_Homo_sapiens_Bisulfite-Seq_1_bismark_bt2_pe.bam

# rrbs data -s (single end)

# important output file is CpG_(OB/OT)_file_name.txt
