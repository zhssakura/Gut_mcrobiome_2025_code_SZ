# 20250604
# Extract part of the file name, and add the sting to the first/last column of each the file.
# 20250611 updated.
import os
import glob
import argparse

def get_filenames_glob(path, pattern):
    filenames = []
    for filename in glob.glob(os.path.join(path, pattern)):
        if os.path.isfile(filename):
             filenames.append(os.path.basename(filename))
    return filenames



######################  run script in command line.
parser = argparse.ArgumentParser()
parser.add_argument('-E.g',
                    required=False,
                    help='python3 Add_one_column_with_file_name_argv_v2.py -dir -format')

parser.add_argument('-dir',
                    required=True,
                    help='/lustre1/g/pharm_jia/Shan_Zhang/HKU/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/6_CoverM/0522/')

parser.add_argument('-format',
                    required=True,
                    default='tsv',
                    help='Pattern of the file at the end.')


args = vars(parser.parse_args())
path  = args['dir']
format= args['format']


##################### Start of running script within pycharm #####################

# Example usage
# path = "./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/8_gapseq/"  # Current directory
# pattern = "*-Pathways.tbl" # Only get tbl files end with "Pathways"
# # pattern = "*.tbl" # Only get all tbl files
# files = get_filenames_glob(path, pattern)
# print(files)

#
# path = "./Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/8_gapseq/"
# format = "Pathways.tbl" # Only get tbl files end with "Pathways"
##################### End of running script within pycharm #####################
pattern = "%s%s" % ('*',format)
path = "%s%s" % (path, '/')
files = get_filenames_glob(path, pattern)

str_4_split_infile = format.strip().split('.')[-1]
for file_name in files:
    # print(file_name)
    file_name_prefix = file_name.strip().split(".%s" % (str_4_split_infile))[0]
    sample_ID = file_name.strip().split("_")[0]
    # print(mag_ID)
    outfile = ('%s%s%s' % (path, file_name_prefix, "_sampleID_added.tsv"))
    # print(outfile)
    # print('\n')
    output_handel = open(outfile,"w")
    for each in open('%s%s' % (path, file_name)):
        # print(each)

        each_split = each.strip().split('\t')
        if each.startswith("Genome	"):
            new_header = '\t'.join(each_split) + "\tSampleID\n"
            output_handel.write(new_header)
        else:
            # print(len(each_split))
            new_each = '%s\t%s\n' % ('\t'.join(each_split), sample_ID)
            output_handel.write(new_each)

    output_handel.close()
print("Done! you have the output files named *_sampleID_added.tsv")
