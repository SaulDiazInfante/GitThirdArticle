# ==============================================================================
# This code filter the data respect to (OC) and (OB) levels.
# INPUT:
#   U1PahtsLongTime.npy,
#   U2PahtsLongTime.npy,
#   binary files with np-arrays of paths.
# OUTPUT:
#   U1PathsLongTimeFiltered.npy,
#   U2PathsLongTimeFiltered.npy,
#   binary files with the filtered data
# =============================================================================
import numpy as np
import os
import csv
#
u1_short_file_name_list = []
with open('u1_short_data_file_list.txt') as inputfile:
    for row in csv.reader(inputfile):
        u1_short_file_name_list.append(row)
u2_short_file_name_list = []
with open('u2_short_data_file_list.txt') as inputfile:
    for row in csv.reader(inputfile):
        u2_short_file_name_list.append(row)
u1_long_file_name_list = []
with open('u1_long_data_file_list.txt') as inputfile:
    for row in csv.reader(inputfile):
        u1_long_file_name_list.append(row)
u2_long_file_name_list = []
with open('u2_long_data_file_list.txt') as inputfile:
    for row in csv.reader(inputfile):
        u2_long_file_name_list.append(row)

r_OB = 9000
r_OC = 20

for j in np.arange(len(u1_long_file_name_list)):
    file_name_u1_short = u1_short_file_name_list.pop()[0]
    file_name_u2_short = u2_short_file_name_list.pop()[0]
    file_name_u1_long = u1_long_file_name_list.pop()[0]
    file_name_u2_long = u2_long_file_name_list.pop()[0]
    try:
        file_prefix_u1_short, ext = os.path.splitext(file_name_u1_short)
        file_prefix_u2_short, ext = os.path.splitext(file_name_u2_short)
        file_prefix_u1_long, ext = os.path.splitext(file_name_u1_long)
        file_prefix_u2_long, ext = os.path.splitext(file_name_u2_long)
    except:
        file_prefix_u1_short = ''
        file_prefix_u2_short = ''
        file_prefix_u1_long = ''
        file_prefix_u2_long = ''

    save_file_name_u1_short = file_prefix_u1_short + '_filtered.npy'
    save_file_name_u2_short = file_prefix_u2_short + '_filtered.npy'
    save_file_name_u1_long = file_prefix_u1_long + '_filtered.npy'
    save_file_name_u2_long = file_prefix_u1_long + '_filtered.npy'

    u_p1_short = np.load(file_name_u1_short)
    u_p1_long = np.load(file_name_u1_long)
    k = u_p1_short.shape[0]
    index = np.arange(k)
    trash_u1_short = []
    trash_u1_long = []
    for i in index:
        if u_p1_short[i, :].max() >= r_OC:
            trash_u1_short.append(i)
        if u_p1_long[i, :].max() >= r_OC:
            trash_u1_long.append(i)
    trash_u1_shortArr = np.array(trash_u1_short)
    trash_u1_longArr = np.array(trash_u1_long)
    print'\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    print '\tOC trash paths {: } of {: }'.format(trash_u1_shortArr.shape[0],
                                               u_p1_short.shape[0])
    del u_p1_short
    del u_p1_long
    #
    u_p2_short = np.load(file_name_u2_short)
    u_p2_long = np.load(file_name_u2_long)
    trash_u2_short = []
    trash_u2_long = []
    for i in index:
        if u_p2_short[i, :].max() >= r_OB:
            trash_u2_short.append(i)
        if u_p2_long[i, :].max() >= r_OB:
            trash_u2_long.append(i)
    trash_u2_shortArr = np.array(trash_u2_short)
    trash_u2_longArr = np.array(trash_u2_long)
    proportion_OB_short = trash_u2_shortArr.shape[0]
    proportion_OB_long = trash_u2_longArr.shape[0]

    trash_index_short = np.union1d(trash_u1_shortArr, trash_u2_shortArr)
    trash_index_long = np.union1d(trash_u1_longArr, trash_u2_longArr)
#
    filter_index_short = np.setdiff1d(index, trash_index_short)
    filter_index_long = np.setdiff1d(index, trash_index_long)
    filter_data_u_p2_short = u_p2_short[filter_index_short]
    filter_data_u_p2_long = u_p2_short[filter_index_long]
    np.save(save_file_name_u2_short, filter_data_u_p2_short)
    np.save(save_file_name_u2_long, filter_data_u_p2_long)
    proportion_OB_short = 1.0 - np.float(filter_data_u_p2_short.shape[0]) \
                       / np.float(u_p2_short.shape[0])
    proportion_OB_long = 1.0 - np.float(filter_data_u_p2_long.shape[0]) \
                                / np.float(u_p2_long.shape[0])
    print '\t OB short paths to trash  {: }: of {: }'. \
        format(proportion_OB_short, u_p2_short.shape[0])
    print '\t Acceptable short paths: {:d}'.\
        format(filter_data_u_p2_short.shape[0])
    print '\t Acceptable long paths: {:d}'.\
        format(filter_data_u_p2_long.shape[0])
    print '\tProportion of short rejections paths: {:1.8f}%'.\
        format(100.0 * proportion_OB_short)
    print '\tProportion of long rejections paths: {:1.8f}%'. \
        format(100.0 * proportion_OB_long)
    print'\t----------------------------------------------------------'
    del u_p2_short
    del filter_data_u_p2_short
    u_p1_short = np.load(file_name_u1_short)
    filter_data_u_p1_short = u_p1_short[filter_index_short]
    filter_data_u_p1_long = u_p1_short[filter_index_long]
    np.save(save_file_name_u1_short, filter_data_u_p1_short)
    del u_p1_short
    del filter_data_u_p1_short
