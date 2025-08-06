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
trimmomatic PE ERR10114000_1.fastq.gz ERR10114000_2.fastq.gz ERR10114000_1_P.fq.gz ERR10114000_1_UP.fq.gz ERR10114000_2_P.fq.gz ERR10114000_2_UP.fq.gz ILLUMINACLIP:/software/trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:6:true LEADING:20 TRAILING:20 HEADCROP:5 SLIDINGWINDOW:4:20 MINLEN:50 -threads 16

#####################################################################################################
###################################### SPAdes v4.2.0 ################################################
#####################################################################################################
spades.py --meta --only-assembler -k 21,41,61,81,101 -1  ERR10114000_1_P.fastq -2  ERR10114000_2_P.fastq -o $out_dir -t 16

#####################################################################################################
######################################## metaWRAP binning ###########################################
#####################################################################################################
metawrap binning -t 36 --metabat2 --maxbin2 --concoct -a ./3_Assembly/SPAdes/ERR10114000/scaffolds.fasta -o ./4_Binning/ERR10114000_metaWRAP_wd ./2_Trimmomatic/ERR10114000_1.fastq ./2_Trimmomatic/ERR10114000_2.fastq

#####################################################################################################
######################################### bin refinement ############################################
#####################################################################################################
metawrap bin_refinement -o ./4_Binning/ERR10114000_metaWRAP_wd/refine_wd -t 36 -c 50 -x 5 -A ./4_Binning/ERR10114000_metaWRAP_wd/metabat2_bins -B ./4_Binning/ERR10114000_metaWRAP_wd/maxbin2_bins -C ./4_Binning/ERR10114000_metaWRAP_wd/concoct_bins

# The outputs are here: /4_Binning/C149_metaWRAP_wd/refine_wd/metawrap_50_5_bins/

#####################################################################################################
################################ GTDB-Tk v2.4.1 (using 200GB memory) ################################
################################ GTDB release 226 ###################################################
#####################################################################################################
export GTDBTK_DATA_PATH=/scr/u/shanbio/Database/GTDB/release226/release226
gtdbtk classify_wf --cpus 36 --pplacer_cpus 1 --genome_dir refined_bin_renamed_20250612 --extension fna --skip_ani_screen --out_dir refined_bin_renamed_GTDB_r226_20250612 --prefix refined_bin_renamed_GTDB_r226

#####################################################################################################
######################## CoverM v0.7.0 estimate bin coverage across samples #########################
#####################################################################################################
coverm genome --bam-files ./4_Binning/ERR10114000_metaWRAP_wd/work_files/ERR10114000.bam --genome-fasta-directory ./4_Binning/ERR10114000_metaWRAP_wd/refine_wd/metawrap_50_5_bins/ -o metawrap_50_5_bins_coverm.tsv -t 36


#####################################################################################################
######################################### checkm2 ################################################### 
#####################################################################################################
checkm2 predict --force --threads 36 --input refined_bin_renamed -x fna --output-directory refined_bin_renamed_checkm2

#################################################################################################################################
######################## dRep v3.6.2 De-replication of microbial genomes assembled from multiple samples. #######################
#################################################################################################################################
dRep dereplicate refined_bin_renamed_dRep_ANI97 -pa 0.9 -sa 0.97 -comp 50 -p 36 -g refined_bin_renamed_20250528/*.fna --genomeInfo refined_bin_renamed_checkm2_quality_for_drep.tsv --multiround_primary_clustering --primary_chunksize 5000 --greedy_secondary_clustering --run_tertiary_clustering

########################################################################################################################################################
######################## gapseq v1.4.0 Informed prediction and analysis of bacterial metabolic pathways and genome-scale networks ######################
########################################################################################################################################################
declare -a keywords=("PWY0-1314" "GLYCOLYSIS" "PWY-5484" "fructokinase" "2.7.1.7" "PYRUVDEHYD-PWY" "PWY-8275" "PWY-5938" "PWY-5939" "PWY-8274" "PWY-6389" "P41-PWY" "PWY-5482" "PWY-5485" "PWY-5537" "PWY-5768" "PWY3O-440" "PWY-6588" "CENTFERM-PWY" "PWY-6583" "PWY-6883" "PWY-5480" "PWY-5486" "PWY-6587" "PWY-6863" "PWY-7111" "P108-PWY" "PWY-7545" "PWY-7544" "PWY-8450" "PWY-7351" "P125-PWY" "PWY-6396" "PWY-5464" "PWY4LZ-257" "GLYCOLYSIS-TCA-GLYOX-BYPASS" "PWY-5742" "PWY-6970" "PWY-5481" "PWY-5096" "PWY-5100" "P142-PWY" "PWY-5483" "PWY-5538" "PWY-5600")

for i in "${keywords[@]}"
do
gapseq find -p "$i" -t Bacteria -b 100 -c 70 -l MetaCyc -n -y ERR10114000.fna > ERR10114000-"$i"-report.sh
done



#########################################################################################################################################################
######################## Merge all *-Pathway.tbl for getting the data matrix to plot nmds, pheatmap and barplot t-test ##################################
#########################################################################################################################################################
# merge all *-Pathway.tbl

# on MacOS:
python3 ./PycharmProjects/Gut_microbiomes/gapseq/Add_one_column_with_file_name_argv.py -dir ./8_gapseq/ -format Pathways.tbl
cd ./8_gapseq
cat *-Pathways_MAG_added.tsv > 00_merged-Pathways_MAG_added.tsv
python3 ./PycharmProjects/Gut_microbiomes/gapseq/Reformat_matrix_using_dod_dict_of_dict_and_dict_of_list_argv.py -in ./8_gapseq/find/00_merged-Pathways_MAG_added.tsv


# on HPC2021 @20250610:
# add additional column to a file:
python3 /lustre1/g/pharm_jia/Shan_Zhang/softwares/PycharmProjects/Gut_microbiomes/gapseq/Add_one_column_with_file_name_argv.py -dir ./8_gapseq/ -format Pathways.tbl
python3 /lustre1/g/pharm_jia/Shan_Zhang/softwares/PycharmProjects/Gut_microbiomes/gapseq/Add_one_column_with_file_name_argv.py -dir ./8_gapseq_Arc/ -format Pathways.tbl

cd ./8_gapseq
rm -rf 00_merged-Pathways_MAG_added.tsv
cat *-Pathways_MAG_added.tsv > 00_merged-Pathways_MAG_added.tsv
scp -r shanbio@hpc2021-io1:./8_gapseq/00_merged-Pathways_MAG_added.tsv ./8_gapseq/find/0609/

# on MacOS
# get 00_merged-Pathways_MAG_added_reform_matrix_annotate.tsv
python3 ./PycharmProjects/Gut_microbiomes/gapseq/Reformat_matrix_using_dod_dict_of_dict_and_dict_of_list_argv.py -in -metadata

# get 00_merged-Pathways_MAG_added_reform_matrix_annotate.tsv transposed (00_merged-Pathways_MAG_added_reform_matrix_annotate_t.tsv), and this output file can be copy to Excel and add color. This is to test if any pwy prediction was missed.
python3 ./PycharmProjects/Gut_microbiomes/Transpose_tsv_using_pandas.py

# manually get 00_merged-Pathways_Bac_Arc_MAG_added_reform_matrix_annotate_rm_duplct_t.tsv

# get 00_merged-Pathways_MAG_added_reform_matrix_annotate_t_vs_ANI_97.tsv (only to check if any MAGs were missed).
python3 ./PycharmProjects/Gut_microbiomes/Match_ID_from_Factor_file_2_Key_file.py

# get 00_merged-Pathways_Bac_Arc_MAG_added_reform_matrix_annotate_t_pair_Cdb_ANI97.tsv and 00_merged-Pathways_Bac_Arc_MAG_added_reform_matrix_annotate_t_pair_Cdb_ANI97_pair_PWY.tsv
python3 ./PycharmProjects/Gut_microbiomes/dRep/Pair_rMAG_n_all_MAGs.py

# get 00_merged-Pathways_Bac_Arc_MAG_added_reform_matrix_annotate_t_pair_Cdb_ANI97_pair_PWY_filtered.tsv
python3 ./PycharmProjects/Gut_microbiomes/gapseq/Filter_out_lines_by_row_names.py

# CAUTION: Depending on using the 00_All_metawrap_50_5_bins_coverm_sampleID_added.tsv file or 00_All_metawrap_50_5_bins_coverm_sampleID_added_ignore_unmapped.tsv, the next output will consider / not consider unmapped reads.
# get: 00_All_metawrap_50_5_bins_coverm_sampleID_added_ignore_unmapped.tsv
python3 ./PycharmProjects/Gut_microbiomes/CoverM/Calculate_rel_abun_of_MAGs_regardless_of_unmapped_reads.py

# get 00_merged-Pathways_Bac_Arc_MAG_added_reform_matrix_annotate_t_pair_Cdb_ANI97_pair_PWY_filtered_add_rel_abun_CoverM20250522_by_coeff.tsv by adding relative abundance and multiplying the relative abundance as coefficient.
python3 ./PycharmProjects/Gut_microbiomes/Match_ID_from_Factor_file_2_Key_file_multiply_coefficient.py

# get 00_merged-Pathways_Bac_Arc_MAG_added_reform_matrix_annotate_t_pair_Cdb_ANI97_pair_PWY_filtered_add_rel_abun_CoverM20250522_by_coeff.tsv transposed (00_merged-Pathways_Bac_Arc_MAG_added_reform_matrix_annotate_t_pair_Cdb_ANI97_pair_PWY_filtered_add_rel_abun_CoverM20250522_by_coeff_t.tsv), Again, this output file can be copy to Excel and add color.
python3 ./PycharmProjects/Gut_microbiomes/Transpose_tsv_using_pandas.py

### get 00_merged-Pathways_Bac_Arc_MAG_added_reform_matrix_annotate_t_pair_Cdb_ANI97_pair_PWY_filtered_add_rel_abun_CoverM20250522_by_coeff_t_per_sample.tsv, this file will be used for nMDS
python3 ./PycharmProjects/Gut_microbiomes/gapseq/Sum_up_values_in_multi_columns.py



#############################################################################
######################## Plotting on MacOS ##################################
#############################################################################
# nmds:
./Scripts/Rscript_NMDS_gut_microbiome.R

# pheatmap
./Scripts/Pheatmap_Gut_Microbiome.R

# t-tests for multiple groups
./Scripts/Rscript_Gut_Microbiome_t_test_for_mutiple_groups_upgraded.R

#############################################################################
######################## Plotting on HPC2021 ################################
#############################################################################
# To get zOTU table by dRep based on ANI!
# When using 97% similarity, ERR10113288_bin.21.fna and ERR10113288_bin.29.fna, are both from the same sample (ERR10113288_bin) and same zOTU (cluster_172_1),
# was clustered into the same zOTU. Therefore, I manually corrected by merging their relative abundance in:
# 00_All_metawrap_50_5_bins_coverm_sampleID_added_zOTU4_RelAbun.tsv


# get zOTU table:


# Community Composition (stacked bar plot)
# my data
cd ./6_CoverM/0522/00_All_metawrap_50_5_bins_coverm_sampleID_added_zOTU4_RelAbun_Stacked_bar_plot
python3 Stacked_bar_plot.py -m metadata.txt -otu OTU_Table.txt -otu_c OTU_Taxa.txt -w 14 -hr "d" -mr "d" -o Gut_microbiome_community_composition_d.pdf -sample interested_sample.txt
python3 Stacked_bar_plot.py -m metadata.txt -otu OTU_Table.txt -otu_c OTU_Taxa.txt -w 14 -hr "d" -mr "p" -o Gut_microbiome_community_composition_p.pdf -sample interested_sample.txt
python3 Stacked_bar_plot.py -m metadata.txt -otu OTU_Table.txt -otu_c OTU_Taxa.txt -w 14 -hr "d" -mr "c" -o Gut_microbiome_community_composition_c.pdf -sample interested_sample.txt
python3 Stacked_bar_plot.py -m metadata.txt -otu OTU_Table.txt -otu_c OTU_Taxa.txt -w 14 -hr "d" -mr "o" -o Gut_microbiome_community_composition_o.pdf -sample interested_sample.txt
python3 Stacked_bar_plot.py -m metadata.txt -otu OTU_Table.txt -otu_c OTU_Taxa.txt -w 25 -hr "d" -mr "f" -o Gut_microbiome_community_composition_f.pdf -sample interested_sample.txt
python3 Stacked_bar_plot.py -m metadata.txt -otu OTU_Table.txt -otu_c OTU_Taxa.txt -w 40 -hr "d" -mr "g" -o Gut_microbiome_community_composition_g.pdf -sample interested_sample.txt


cd ./6_CoverM/0522/00_All_metawrap_50_5_bins_coverm_sampleID_added_zOTU4_RelAbun_Stacked_bar_plot/Gut_microbiome_community_composition__NMDS/
python3 NMDS.py
