#####################################################################################################
################################## Download data on HKU HPC2021 #####################################
#####################################################################################################
ena-file-download-read_run-PRJEB55534-fastq_ftp-20250501-1123.sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR10114000/ERR10114000_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR10114000/ERR10114000_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR10114100/ERR10114100_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR10114100/ERR10114100_2.fastq.gz
# ... ...

#################################### FastQC v0.11.8 #################################################
fastqc ERR10114000_1.fastq.gz -o ./ERR10114000_1_fastqc -t 16
fastqc ERR10114000_2.fastq.gz -o ./ERR10114000_2_fastqc -t 16
# ... ...

#####################################################################################################
#################################  Trimmomatic v0.39  ############################################### 
#####################################################################################################
trimmomatic PE ERR10114000_1.fastq.gz ERR10114000_2.fastq.gz ERR10114000_1_Q25_P.fq.gz ERR10114000_1_Q25_UP.fq.gz ERR10114000_2_Q25_P.fq.gz ERR10114000_2_Q25_UP.fq.gz ILLUMINACLIP:/software/trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:6:true LEADING:20 TRAILING:20 HEADCROP:5 SLIDINGWINDOW:4:20 MINLEN:50 -threads 16

