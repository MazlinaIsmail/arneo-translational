# arneo-translational
Methods and collection of scripts used for data analysis, including to produce the figures in manuscript.

Section 1: Processing RNAseq data on UCL Myriad
1. Run FastQC and multiqc on raw fastqs
2. Adaptor trim using Trimmomatic
3. Alignment and quantification using RSEM-STAR

Create sample list for array job submission

DIR=/path/to/fastq/folder
PROJECT="project-name"
TEMP=/path/to/R1FQ_temp_file.txt
find $DIR -name "*_R1_*fastq\.gz" | sort -n >>$TEMP; cat -n $TEMP >>/path/to/R1FQ_"$PROJECT".txt; rm $TEMP


