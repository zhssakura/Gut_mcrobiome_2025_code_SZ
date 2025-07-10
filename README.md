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


