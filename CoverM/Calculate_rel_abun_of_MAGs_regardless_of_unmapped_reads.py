# @20250615 HKUST
# This script used the output of CoverM to get new relative abundance of qualified MAGs (based on CheckM 2 results) regardless of unmapped reads.
infile_coverm = "/Users/zzfanyi/Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/6_CoverM/0522/00_All_metawrap_50_5_bins_coverm_sampleID_added.tsv"
infile_checm2_qualified_MAGs = "/Users/zzfanyi/Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/7_CheckM2/quality_report_WS.tsv"
outfile_coverm = "/Users/zzfanyi/Documents/07_HKU/03_Research/01_PRJ_NAFLD_gut_microbiome/1_Databases/European_Nucleotide_Archive_ENA/PRJEB55534/6_CoverM/0522/00_All_metawrap_50_5_bins_coverm_sampleID_added_ignore_unmapped.tsv"

# st1. get a list of qualified MAG ids:
set_qualified_mags = set()
for each in open(infile_checm2_qualified_MAGs):
    each_split = each.strip().split("\t")
    if not each.startswith("Name\t"):
        MAG_id = each_split[0]
        set_qualified_mags.add(MAG_id)
print("%s%s%s" % ("There are ", len(set_qualified_mags), " qualified MAGs in CheckM2 results (comp. > 50%)."))

# st2. get dict of list:
dict_sample_2_qMAGs = dict()
dict_sample_2_qMAGs_RelAbun = dict()
dod_sample_2_qMAG_2_rel_abun = dict()
for each in open(infile_coverm):
    each_split = each.strip().split("\t")
# if not each.startswith("Genome\t"):
    MAG_id = each_split[3]
    Sample_id = each_split[2]
    MAG_rel_abun = each_split[1]

    # IF the current sample id is not in the dict, the facter should be the list of current mag when giving current sample ID:
    if Sample_id not in dict_sample_2_qMAGs:
        # only when MAGs were with completeness assessed by CheckM2 > 50%, will it be calculated for new relative abundance:
        if MAG_id in set_qualified_mags:
            qMAG_id = MAG_id
            qMAG_rel_abun = MAG_rel_abun

            dict_sample_2_qMAGs[Sample_id] = list()
            dict_sample_2_qMAGs[Sample_id].append(qMAG_id)

            dict_sample_2_qMAGs_RelAbun[Sample_id] = list()
            dict_sample_2_qMAGs_RelAbun[Sample_id].append(qMAG_rel_abun)

            dod_sample_2_qMAG_2_rel_abun[Sample_id] = {}
            dod_sample_2_qMAG_2_rel_abun[Sample_id][qMAG_id] = qMAG_rel_abun

    else:
        # only when MAGs were with completeness assessed by CheckM2 > 50%, will it be calculated for new relative abundance:
        if MAG_id in set_qualified_mags:
            qMAG_id = MAG_id
            qMAG_rel_abun = MAG_rel_abun

            dict_sample_2_qMAGs[Sample_id].append(qMAG_id)

            dict_sample_2_qMAGs_RelAbun[Sample_id].append(qMAG_rel_abun)

            dod_sample_2_qMAG_2_rel_abun[Sample_id][qMAG_id] = qMAG_rel_abun

        else:
            pass

    # Sumerize total value of rel_abun of all MAGs under current sample ID.
    # sum_qMAG_rel_abun = sum(float(i) for i in dict_sample_2_qMAGs_RelAbun[Sample_id])

    # print("%s\t%s\t%s\t%s" % (Sample_id, list_current_qMAG_rel_abun, qMAG_id, dict_sample_2_qMAGs_RelAbun[Sample_id]))
    # print('\n')
# print(Sample_id)

print(dict_sample_2_qMAGs.keys())
print(dict_sample_2_qMAGs["H2WU"])
# print(dict_sample_2_qMAGs_RelAbun)
print(dict_sample_2_qMAGs_RelAbun.keys())
print(dict_sample_2_qMAGs_RelAbun["H2WU"])
# print(dict_sample_2_qMAGs_RelAbun.values())

# print(dod_sample_2_qMAG_2_rel_abun)

# st3. Sumerise total value of rel_abun of all MAGs under each sample ID.
# make a dict of qmag to new rel_abun
dict_qmag_2_new_rel_abun = dict()
for each_sample in dict_sample_2_qMAGs_RelAbun.keys():
    # print(each_sample)
    sum_qMAG_rel_abun = sum(float(i) for i in dict_sample_2_qMAGs_RelAbun[each_sample])
    # print("%s\t%s\t%s" % (each_sample, sum_qMAG_rel_abun, dict_sample_2_qMAGs_RelAbun[each_sample]))
    for each_sample_mag in dod_sample_2_qMAG_2_rel_abun[each_sample]:
        current_mag_rel_abun = dod_sample_2_qMAG_2_rel_abun[each_sample][each_sample_mag]
        current_mag_rel_abun_new = str(100 * (float(current_mag_rel_abun) / float(sum_qMAG_rel_abun)))
        # print("%s\t%s\t%s\t%s" % (each_sample, each_sample_mag, current_mag_rel_abun, sum_qMAG_rel_abun))

        # calculate new relative abundance of each mags:
        current_to_print = ("%s\t%s\t%s" % (each_sample, each_sample_mag, current_mag_rel_abun_new))
        # print(current_to_print)
        dict_qmag_2_new_rel_abun[each_sample_mag]=current_mag_rel_abun_new
print("%s %s %s" % ("There are ", len(dict_qmag_2_new_rel_abun), " qMAGs with new rel_abun values."))


# st4. generate a new file in the format of coverm output.
output_handel = open(outfile_coverm, "w")
for each in open(infile_coverm):
    each_split = each.strip().split('\t')
    if each.startswith("Genome\t"):
        header = each
        output_handel.write(header)
    else:
        mag = each_split[0]
        sample_id = each_split[2]
        unqMAG = each_split[3]

        if unqMAG in dict_qmag_2_new_rel_abun:
            unqMAG_rel_abun = float(dict_qmag_2_new_rel_abun[unqMAG])
            unqMAG_rel_abun = round(unqMAG_rel_abun,8)

            # print("%s %s %s" % ("qMAG ", unqMAG, unqMAG_rel_abun))
            for_print = ("%s\t%s\t%s\t%s\n" % (mag, unqMAG_rel_abun, sample_id, unqMAG))
            output_handel.write(for_print)
        else:
            unqMAG_rel_abun = 0
            # print("%s %s %s" % ("unqMAG ", unqMAG, unqMAG_rel_abun))
            for_print = ("%s\t%s\t%s\t%s\n" % (mag, unqMAG_rel_abun, sample_id, unqMAG))
            output_handel.write(for_print)
output_handel.close()
print("Done!")
