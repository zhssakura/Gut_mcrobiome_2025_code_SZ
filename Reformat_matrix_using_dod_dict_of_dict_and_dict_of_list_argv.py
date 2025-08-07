# @20250605
# Apply this script to the output of gapseq -find *-Pathways.tbl files.
# Reformat the matrix below using dict of dict (Dod):

# ID	Completeness	MagID
# |CENTFERM-PWY|	0	C149_bin.3
# |P108-PWY|	28	ERR10113960_bin.44
# |P125-PWY|	60	ERR10113960_bin.44
# |P125-PWY|	15	ERR10111111_bin.11
# |PWY-3861|	33	C149_bin.3
# |PWY-6140|	0	C149_bin.3
# |PWY-6145|	17	C149_bin.3
# |PWY-6992|	40	C149_bin.3
# |PWY3O-1743|	50	C149_bin.3
# 									
# 									
# 									
# 									
# MAG_ID	C149_bin.3	ERR10111111_bin.11	ERR10113960_bin.44
# |PWY-3861|	33	NA	NA
# |PWY3O-1743|	50	NA	NA
# |PWY-6140|	0	NA	NA
# |PWY-6145|	17	NA	NA
# |PWY-6992|	40	NA	NA
# |CENTFERM-PWY|	0	NA	NA
# |P108-PWY|	NA	NA	28
# |P125-PWY|	NA	15	60

# An simplified nested dict:
# https://stackoverflow.com/questions/16333296/how-do-you-create-nested-dict-in-python
# >>> d = {}
# >>> d['dict1'] = {}
# >>> d['dict1']['innerkey'] = 'value'
# >>> d['dict1']['innerkey2'] = 'value2'
# >>> d
# {'dict1': {'innerkey': 'value', 'innerkey2': 'value2'}}

import argparse

######################  run script in command line.
parser = argparse.ArgumentParser()
parser.add_argument('-E.g',
                    required=False,
                    help='python3 Reformat_matrix_using_dod_dict_of_dict_and_dict_of_list_argv.py -in ')

parser.add_argument('-in',
                    required=True,
                    help='./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/8_gapseq/find/00_merged-Pathways_MAG_added.tsv')

parser.add_argument('-metadata',
                    required=True,
                    help='./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/0_Metadata/00_MAGs_ANI_97_vs_00_gtdb_of_all_MAGs_taxon.txt')



args = vars(parser.parse_args())
file_in  = args['in']
metadata  = args['metadata']


##################### Start of running script within pycharm #####################
# mydir = "./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/8_gapseq/find/"
# file_in             = mydir + "00_merged-Pathways_MAG_added.tsv"
# metadata           = "./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/0_Metadata/00_MAGs_ANI_97_vs_00_gtdb_of_all_MAGs_taxon.txt"
# 
# mydir = "./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/8_gapseq/find/0609/"
# file_in             = mydir + "00_merged-Pathways_MAG_added.tsv"
# metadata           = "./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/0_Metadata/00_MAGs_ANI_97_vs_00_gtdb_of_all_MAGs_taxon.txt"
# 
# mydir = "./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/8_gapseq/find/0610/"
# file_in             = mydir + "00_merged-Pathways_Bac_Arc_MAG_added.tsv"
# metadata           = "./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/0_Metadata/00_MAGs_ANI_97_vs_00_gtdb_of_all_MAGs_taxon.txt"
##################### End of running script within pycharm #####################


file_out            = file_in.strip().split('.')[0] + "_reform_matrix.tsv" # This is reformated dataframe.
file_out_annotate   = file_in.strip().split('.')[0] + "_reform_matrix_annotate.tsv"


# step 1: created dod (dict of dict)
dict_of_dict_pwy_2_comp = {}
dict_pwy_2_des = {}
set_mags = set()
for each in open(file_in):
    each_split = each.strip().split('\t')

    if not each.startswith("ID"):
        pwy_id = each_split[0]
        pwy_des = each_split[1]
        # pwy_comp = each_split[3]# by uncomment this line, it provides a final table of number (unit: %).
        pwy_comp = each_split[2].replace("true", "1").replace("false", "0")# by uncomment this line, it provides a final table of Presence/Absence represented by "1/0" (or TRUE and FALSE).
        print(pwy_comp)
        mag_id = each_split[-2]
        set_mags.add(mag_id)
        dict_pwy_2_des[pwy_id] = pwy_des

        # print(each_split)
        if pwy_id not in dict_of_dict_pwy_2_comp:
            dict_of_dict_pwy_2_comp[pwy_id] = {}

            if mag_id not in dict_of_dict_pwy_2_comp[pwy_id]:
                dict_of_dict_pwy_2_comp[pwy_id][mag_id] = pwy_comp
            else:
                dict_of_dict_pwy_2_comp[pwy_id][mag_id] += pwy_comp

        else:
            dict_of_dict_pwy_2_comp[pwy_id][mag_id] = pwy_comp
# print(pwy_comp_dict_of_dict)


# Step 2. sort set:
set_mags_sorted = sorted(set_mags)
# print(set_mags_sorted)


# Step 3: reformat the matrix using disk of list:
disk_mag_2_comp = {}

output_handel = open(file_out, "w")
# print("%s\t%s\n" % ("MAG_ID", '\t'.join(set_mags_sorted)))
output_handel.write("%s\t%s\n" % ("MAG_ID", '\t'.join(set_mags_sorted)))

for each_pwy in dict_of_dict_pwy_2_comp:
    list_mag = []
    list_comp = []
    # print(each_pwy)
    for each_mag in set_mags_sorted:

        if each_mag in dict_of_dict_pwy_2_comp[each_pwy]:
            pwy_comp = dict_of_dict_pwy_2_comp[each_pwy][each_mag]
            list_mag.append(each_mag)
            list_comp.append(pwy_comp)
        else:
            pwy_comp= "NA"
            list_mag.append(each_mag)
            list_comp.append(pwy_comp)
    # print(list_mag)
    # print(list_comp)
    disk_mag_2_comp[each_pwy] = list_comp

    # print(disk_mag_2_comp)
    # print("%s\t%s\n" % (each_pwy, '\t'.join(list_comp)))
    output_handel.write("%s\t%s\n" % (each_pwy, '\t'.join(list_comp)))
output_handel.close()
print("Done! Reformat complete and doing annotation next.")

# Step 4: Add pwy annotation and taxonomy classification to the reformated dataframe:
# Generate dict of using metadata:
dict_mag_2_name = {}
for each in open(metadata):
    # print(each)
    each_split = each.strip().split('\t')
    if not each.startswith('user_genome	'):
        mag_id = each_split[0]
        mag_name = each_split[1]
        dict_mag_2_name[mag_id] = mag_name

print("Done! Get MAG names by providing MAG ID.")

# Apply annotation and get final output file:
list_mag_names = []

output_handel = open(file_out_annotate, 'w')
for each in open(file_out):
    each_split = each.strip().split('\t')
    if each.startswith("MAG_ID"):
        for each_mag in each_split:
            if len(each_mag.strip().split("_Arc")) > 1: # This is to identify if the MAG id is added with "_Arc" by Shan manually.
                mag_id = each_mag.strip().split("_Arc")[0]
                # print(each_mag)
                if mag_id in dict_mag_2_name:
                    mag_name = dict_mag_2_name[mag_id]
                    list_mag_names.append(mag_name)
                else:
                    mag_id = "%s_%s" % (mag_id,"NA")
                    list_mag_names.append(mag_name)
            if len(each_mag.strip().split("_Arc")) == 1:
                mag_id = each_mag
                if mag_id in dict_mag_2_name:
                    mag_name = dict_mag_2_name[mag_id]
                    list_mag_names.append(mag_name)
                else:
                    mag_id = "%s_%s" % (mag_id,"NA")
                    list_mag_names.append(mag_name)
        output_handel.write('%s\t%s\n' % (each.strip(), "PWY_Name"))
        output_handel.write('%s\t%s\t%s\n' % ("Taxonomy_name",'\t'.join(list_mag_names[1:]),"PWY_Name"))
    if not each.startswith("MAG_ID"):
        # print(each)
        pwy_id = each.strip().split('\t')[0]
        if pwy_id in dict_pwy_2_des:
            pwy_des = dict_pwy_2_des[pwy_id]
        else:
            pwy_des = "NA"
        # print('%s\t%s\n' % (each.strip(), pwy_des))
        output_handel.write('%s\t%s\n' % (each.strip(), pwy_des))
output_handel.close()
print("Done! Both reformat and annotation has completed.")
print("Run /PycharmProjects/Gut_microbiomes/Transpose_tsv_using_pandas.py")

