# code for processing wgbs samples into bismark
# basically all the slurm steps

# initially had tried using bwa-meth, but the output of that didn't like working with bismark (read level meth)
# so switched to bismark
# ok rrbs and wgbs data I skipped adapter trimming (trimgalore) and deduplication (bismark) because whoopsies
# going to have thos in with emseq

# start with fastq file

# adapter trimming

module load trimgalore

trim_galore --cores 4 --paired file_r1 file_r2

## new problem: trim_galore keeps finding new errors
## and I don't have the brainspace or capacity for this
## risk is false positives and right now I don't care
## I just need the data to be processed
## so back to old way (for right now)
### then next week can say I analysed it but I don't think I trust it - found this pipeline - can you help me get it running

# genome prep (builds index)
module load bismark

bismark_genome_preparation /researchers/caitlin.page

# align
# run from folder fastq is in
module load bismark
module load samtools
module load bowtie2

bismark --parallel 4 -o /researchers/caitlin.page/pancreatic_cancer_wgbs/bismark_res
--genome /researchers/caitlin.page -q {
-1 SRR1646792_GSM1541788_38-Lg_rep1_BS-Seq_Homo_sapiens_Bisulfite-Seq_1.fastq.gz
-2 SRR1646792_GSM1541788_38-Lg_rep1_BS-Seq_Homo_sapiens_Bisulfite-Seq_2.fastq.gz}
# rrbs data just the one fastq file (single end)

bismark

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
