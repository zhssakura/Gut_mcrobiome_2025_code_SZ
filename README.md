# Workflow v1.0.0
# Download Data on HKU HPC2021
```
# Script: ena-file-download-read_run-PRJEB55534-fastq_ftp-20250501-1123.sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR10114000/ERR10114000_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR10114000/ERR10114000_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR10114100/ERR10114100_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR10114100/ERR10114100_2.fastq.gz
# ... add more as needed

# Sownload SRA file
conda activate BioSAK
# MetaT data
cd ./PRJNA1246224_MetaT/1_Raw_data
BioSAK sra -i ./PRJNA1246224_MetaT/0_Metadata/BioSAK/sra_id_AKK_Project_98_3.txt -o AKK_Projuect_98_20251109_metaT_3 -t 16 -maxsize 100G

# MetaG data
cd ./PRJNA1246224_MetaG/1_Raw_data
BioSAK sra -i ./PRJNA1246224_MetaG/0_Metadata/sra_id_AKK_Project_98_metaG.txt -o AKK_Projuect_98_20251115_metaG -t 16 -maxsize 100G
```

# Convert SRA files to FASTQ Format.
```
conda activate BioSAK
cd ./PRJNA1246224_MetaT/1_Raw_data
fasterq-dump ./PRJNA1246224_MetaT/1_Raw_data/AKK_Projuect_98_20251109_metaT_3/SRS24675716/SRR33081990/SRR33081990.sra --split-3 -O AKK_Projuect_98_20251109_metaT_3 -t AKK_Projuect_98_20251109_metaT_3/fasterq_dump_tmp
```


# Quality Control with FastQC (v0.11.8)
```
fastqc ERR10114000_1.fastq.gz -o ./ERR10114000_1_fastqc -t 16
fastqc ERR10114000_2.fastq.gz -o ./ERR10114000_2_fastqc -t 16
# ... add more as needed
```
# Trimming Reads with Trimmomatic (v0.39)
```
trimmomatic PE ERR10114000_1.fastq.gz ERR10114000_2.fastq.gz \
ERR10114000_1_P.fq.gz ERR10114000_1_UP.fq.gz \
ERR10114000_2_P.fq.gz ERR10114000_2_UP.fq.gz \
ILLUMINACLIP:/software/trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:6:true \
LEADING:20 TRAILING:20 HEADCROP:5 SLIDINGWINDOW:4:20 MINLEN:50 -threads 16
```

# Assembly with SPAdes (v4.2.0)
```
spades.py --meta --only-assembler -k 21,41,61,81,101 \
-1 ERR10114000_1_P.fastq -2 ERR10114000_2_P.fastq \
-o $out_dir -t 16
```

# Binning with MetaWRAP (v1.3.2)
```
metawrap binning -t 36 --metabat2 --maxbin2 --concoct \
-a ./3_Assembly/SPAdes/ERR10114000/scaffolds.fasta \
-o ./4_Binning/ERR10114000_metaWRAP_wd \
./2_Trimmomatic/ERR10114000_1.fastq ./2_Trimmomatic/ERR10114000_2.fastq
```

# Bin Refinement
```
metawrap bin_refinement -o ./4_Binning/ERR10114000_metaWRAP_wd/refine_wd \
-t 36 -c 50 -x 5 \
-A ./4_Binning/ERR10114000_metaWRAP_wd/metabat2_bins \
-B ./4_Binning/ERR10114000_metaWRAP_wd/maxbin2_bins \
-C ./4_Binning/ERR10114000_metaWRAP_wd/concoct_bins

# Output location:
/4_Binning/C149_metaWRAP_wd/refine_wd/metawrap_50_5_bins/
```

# GTDB-Tk (v2.4.1) — Taxonomic Classification (requires ~200GB RAM)
```
export GTDBTK_DATA_PATH=/scr/u/shanbio/Database/GTDB/release226/release226
gtdbtk classify_wf --cpus 36 --pplacer_cpus 1 \
--genome_dir refined_bin_renamed_20250612 \
--extension fna --skip_ani_screen \
--out_dir refined_bin_renamed_GTDB_r226_20250612 \
--prefix refined_bin_renamed_GTDB_r226
```

# Coverage Estimation with CoverM v0.7.0
```
coverm genome --bam-files ./4_Binning/ERR10114000_metaWRAP_wd/work_files/ERR10114000.bam \
--genome-fasta-directory ./4_Binning/ERR10114000_metaWRAP_wd/refine_wd/metawrap_50_5_bins/ \
-o metawrap_50_5_bins_coverm.tsv -t 36
```

# Quality Assessment with CheckM2
```
checkm2 predict --force --threads 36 --input refined_bin_renamed \
-x fna --output-directory refined_bin_renamed_checkm2
```

# De-replication with dRep v3.6.2
```
dRep dereplicate refined_bin_renamed_dRep_ANI97 \
-pa 0.9 -sa 0.97 -comp 50 -p 36 \
-g refined_bin_renamed_20250528/*.fna \
--genomeInfo refined_bin_renamed_checkm2_quality_for_drep.tsv \
--multiround_primary_clustering --primary_chunksize 5000 \
--greedy_secondary_clustering --run_tertiary_clustering
```

# Pathway & Metabolic Analysis with gapseq v1.4.0
```
declare -a keywords=("PWY0-1314" "GLYCOLYSIS" "PWY-5484" "fructokinase" "2.7.1.7" "PYRUVDEHYD-PWY" "PWY-8275" "PWY-5938" "PWY-5939" "PWY-8274" "PWY-6389" "P41-PWY" "PWY-5482" "PWY-5485" "PWY-5537" "PWY-5768" "PWY3O-440" "PWY-6588" "CENTFERM-PWY" "PWY-6583" "PWY-6883" "PWY-5480" "PWY-5486" "PWY-6587" "PWY-6863" "PWY-7111" "P108-PWY" "PWY-7545" "PWY-7544" "PWY-8450" "PWY-7351" "P125-PWY" "PWY-6396" "PWY-5464" "PWY4LZ-257" "GLYCOLYSIS-TCA-GLYOX-BYPASS" "PWY-5742" "PWY-6970" "PWY-5481" "PWY-5096" "PWY-5100" "P142-PWY" "PWY-5483" "PWY-5538" "PWY-5600")

for i in "${keywords[@]}"
do
  gapseq find -p "$i" -t Bacteria -b 100 -c 70 -l MetaCyc -n -y ERR10114000.fna > ERR10114000-"$i"-report.sh
done
```


# Conduct KEGG or COG annotation on assembled metagenomic scaffold reads.
```
conda activate BioSAK
cd ./11_Prokka_29/SRS24590646_Bac
BioSAK COG2020 -m P -t 36 -db_dir ~/my_DB/COG2020 -i Spades_assembly.faa
# The output file for next step is ./11_Prokka_29/SRS24590646_Bac/Spades_assembly_COG2020_wd/Spades_assembly_query_to_cog.txt
```





```
```





# Get Transcripts of Specific genes based on COG2020 Annotation Results.
# genes I am looking for are: tdcE, grcA and AdhE
```
grep -e 'Query' -e 'COG1882' -e 'COG3445' -e 'COG1012' -e 'COG1454' ./11_Prokka_29/SRS24590646_Bac/Spades_assembly_COG2020_wd/Spades_assembly_query_to_cog.txt > ./12_Annotation/SRS24590646_Spades_assembly_query_to_cog_grep_COG1882_COG3445_COG1012_COG1454.txt
```








# Merging Pathway Tables for Data Analysis
# On MacOS:
```
python3 ./PycharmProjects/Gut_microbiomes/gapseq/Add_one_column_with_file_name_argv.py -dir ./8_gapseq/ -format Pathways.tbl
cd ./8_gapseq
cat *-Pathways_MAG_added.tsv > 00_merged-Pathways_MAG_added.tsv
python3 ./PycharmProjects/Gut_microbiomes/gapseq/Reformat_matrix_using_dod_dict_of_dict_and_dict_of_list_argv.py -in ./8_gapseq/find/00_merged-Pathways_MAG_added.tsv
```

# On HPC2021:
# Add additional column to files
```
python3 /lustre1/g/pharm_jia/Shan_Zhang/softwares/PycharmProjects/Gut_microbiomes/gapseq/Add_one_column_with_file_name_argv.py -dir ./8_gapseq/ -format Pathways.tbl
python3 /lustre1/g/pharm_jia/Shan_Zhang/softwares/PycharmProjects/Gut_microbiomes/gapseq/Add_one_column_with_file_name_argv.py -dir ./8_gapseq_Arc/ -format Pathways.tbl

cd ./8_gapseq
rm -rf 00_merged-Pathways_MAG_added.tsv
cat *-Pathways_MAG_added.tsv > 00_merged-Pathways_MAG_added.tsv
scp -r shanbio@hpc2021-io1:./8_gapseq/00_merged-Pathways_MAG_added.tsv ./8_gapseq/find/0609/
```

# Additional Data Processing & Visualization
```
Follow the sequential commands provided in your script for further analyses, merging, filtering, and plotting.*
Plotting (MacOS & HPC2021)
# On MacOS:
# NMDS plot
./Scripts/Rscript_NMDS_gut_microbiome.R

#
```

