import numpy as np
import pandas
import string
from scipy import stats
import random
from fastdtw import fastdtw
import math
import matplotlib.pyplot as plt
from numpy import genfromtxt
def graph_plot(x,y,color_format,img_name):
    plt.plot(x, y, color_format)
    plt.grid()
    img_name = 'data/max_pars_trees_w_variable_bins_w_xx_segments/' + img_name
    plt.savefig(str(img_name)+'.png', bbox_inches='tight')
    # plt.show()
def normalise(array, global_max= 0,global_min= 0):
    maxi = global_max
    mini = global_min
    maxi = np.max(array)
    mini = np.min(array)
    denominator = maxi - mini
    norm = (array - mini) / denominator
    return np.asarray(norm)
def string_length(len):
    lower_case = string.ascii_lowercase
    upper_case = string.ascii_uppercase
    numbers = string.printable[0:10]
    special_chars = string.printable[62:94]
    total = list(upper_case + lower_case + numbers + special_chars)
    total.pop(75) #because biopython doesn't accept '.' as one of the characters
    if len <= 93:
        total = total[0:len]
        # print(''.join(total))
        return total
    else:
        # print(''.join(total))
        return total
def array_2_string_automated(sequence, partitions):
    protein_encoding = string_length(partitions)
    # all_ascii_chars = "".join(chr(x) for x in range(33, 127))
    # total_no_of_steps = len(all_ascii_chars)
    # protein_encoding = all_ascii_chars
    steps = np.linspace(0, 1, num=partitions)
    protein_seq = ''
    aug_seq = (sequence * 10000).astype(int)
    aug_steps = (steps*10000).astype(int)
    indexes = []
    for i in range(len(aug_seq)):
        protein_seq = protein_seq + (protein_encoding[np.where(aug_steps <= aug_seq[i])[0][-1]]) #adding the characters corresponding to the value of normalised reflectance data to 'protein seq'
    # print((protein_seq))
    # print('\n')
    return protein_seq
def map(value, leftMin, leftMax, rightMin, rightMax):
    leftSpan = leftMax - leftMin
    rightSpan = rightMax - rightMin
    valueScaled = (value - leftMin) / (leftSpan)
    return int(rightMin + (valueScaled * rightSpan))
def wavelength_wise_segmentation(original_csv, k):
    # input: a csv of reflectance values, number of partitions k
    # output: a list of k no. of bins per segment based on DTW score sum.
    my_data = original_csv
    wavelengths = my_data[0][2:]
    materials = my_data[1:,2:]
    splitted_wavelengths = np.array_split(wavelengths, k)
    # print(type(splitted_wavelengths))
    scores = []
    for i in range(k):
        # print(splitted_wavelengths[i][0], splitted_wavelengths[i][-1])
        mats = []
        for j, material in enumerate(materials):
            material_normalised = normalise(material)
            if normalise_bool:
                mats.append((np.array_split(material_normalised,k)[i]))
            else:
                mats.append(np.array_split(material, k)[i])
        A = []
        for i, i_mat in enumerate(mats):
            # colors = "bgrcmykw"
            # plt.plot((i_mat), c=colors[random.randint(0, 7)])
            for j, j_mat in enumerate(mats):
                if j < i:
                    A.append(0)
                    continue
                distance, path = fastdtw(i_mat, j_mat);
                A.append(distance);
        scores.append(sum(A))
        A = np.reshape(A, (-1, 16))
        # plt.show()
    no_of_bins = []

    percentile_ranks = [int(stats.percentileofscore(scores, item, 'rank')) for item in scores]
    # print('\n')
    # print(scores, "scores")
    # print(percentile_ranks, "percentile ranks")
    # for i, item in enumerate(percentile_ranks): # RANKS ARE DIRECTLY PROPORTIONAL TO 'DTW DISTANCE SUM' OF EACH SEGMENT
    #     no_of_bins.append(map(item,0,100,20,93))

    inverse_normed_scores = [np.round(100*(sum(scores)/float(i)),1) for i in scores]
    inverse_normed_percentile_ranks = [int(stats.percentileofscore(inverse_normed_scores, item, 'rank')) for item in inverse_normed_scores]
    # print(inverse_normed_scores, "normed scores")
    # print(inverse_normed_percentile_ranks, "normed percentile ranks")
    for i, item in enumerate(inverse_normed_percentile_ranks): # RANKS ARE INVERSELY PROPORTIONAL TO 'DTW DISTANCE SUM' OF EACH SEGMENT
        no_of_bins.append(map(item,0,100,20,93))

    print(no_of_bins,'no of bins')
    return no_of_bins
def variable_bin_sequence_generator(original_csv, k):
    my_data = original_csv
    no_of_bins_per_segment = wavelength_wise_segmentation(my_data, k)
    wavelengths = my_data[0][2:]
    materials = my_data[1:,2:]
    splitted_wavelengths = np.array_split(wavelengths, k)
    # print(type(splitted_wavelengths))
    scores = []
    all_materials_sequences = []
    mats = []
    all_mats = []
    for j, material in enumerate(materials):
        one_material_sequence = ''
        one_material_sequence_list = []
        mats = []
        names = []
        for i in range(k):
            material_normalised = normalise(material)
            if normalise_bool:
                if segment_wise_normalise_bool:
                    split_up_normalised_material = normalise(np.array_split(material_normalised, k)[i])
                    mats.append(split_up_normalised_material)
                    if swir_bool: img_name = 'local_normalised_sequence_segment_normalised_swir_variable_bins_'+str(i+1)+'of'+str(k)+'_segments'; names.append(img_name)
                    else: img_name = 'local_normalised_sequence_segment_normalised_vnir_variable_bins_'+str(i+1)+'of'+str(k)+'_segments'; names.append(img_name)
                else:
                    split_up_normalised_material = np.array_split(material_normalised, k)[i]
                    mats.append(split_up_normalised_material)
                    if swir_bool: img_name = 'local_normalised_sequence_swir_variable_bins_'+str(i+1)+'of'+str(k)+'_segments'; names.append(img_name)
                    else: img_name = 'local_normalised_sequence_vnir_variable_bins_'+str(i+1)+'of'+str(k)+'_segments'; names.append(img_name)
                one_material_sequence = one_material_sequence+(array_2_string_automated((split_up_normalised_material),no_of_bins_per_segment[i]))
            else:
                if segment_wise_normalise_bool:
                    split_up_material = normalise(np.array_split(material, k)[i])
                    mats.append(split_up_material)
                    if swir_bool: img_name = 'unnormalised_sequence_segment_normalised_swir_variable_bins_'+str(i+1)+'of'+str(k)+'_segments'; names.append(img_name)
                    else: img_name = 'unnormalised_sequence_segment_normalised_vnir_variable_bins_'+str(i+1)+'of'+str(k)+'_segments'; names.append(img_name)
                else:
                    split_up_material = np.array_split(material, k)[i]
                    mats.append(split_up_material)
                    if swir_bool: img_name = 'unnormalised_sequence_swir_variable_bins_'+str(i+1)+'of'+str(k)+'_segments'; names.append(img_name)
                    else: img_name = 'unnormalised_sequence_vnir_variable_bins_'+str(i+1)+'of'+str(k)+'_segments'; names.append(img_name)
                one_material_sequence = one_material_sequence+(array_2_string_automated((split_up_material),no_of_bins_per_segment[i]))
        all_mats.append(mats)
        all_materials_sequences.append(one_material_sequence)
        # if save_text_bool:
        #     text_file = open("data/" + file_name, "a")
        #     text_file.write(one_material_sequence)
        #     text_file.write('\n')
        #     text_file.close()
    # print(len(all_mats[0]))
    # print((np.asarray(all_mats).shape))
    if img_bool:
        for i in range(k):
            x = splitted_wavelengths[i]
            for j in range(len(materials)):
                # y = np.asarray(all_mats)[:,i][j]
                y = all_mats[j][i]
                img_name = names[i]
                graph_plot(x, y, colors[j], img_name=img_name)
            plt.clf()
    if save_text_bool:
        write_lines(all_materials_sequences)
    return all_materials_sequences
def hyperparameters(swir_bool, normalise_bool, segment_wise_normalise_bool, rgb_bool,variable_bins_segments,file_name):
    if rgb_bool:
        colors = "bgrcmykw"
    else:
        colors = ['yo','ko','ro','bo','yx','rx','bx','kx',
              'r+','y+','b+','k+','y-','k-','r-','b-']
    if normalise_bool:
        if segment_wise_normalise_bool:
            if swir_bool:
                my_data = my_data_swir
                file_name = 'txt_resp_phy_files/local_normalised_sequence_segment_normalised_swir_variable_bins_' + str(variable_bins_segments) + '_segments.txt'
            else:
                my_data = my_data_vnir
                file_name = 'txt_resp_phy_files/local_normalised_sequence_segment_normalised_vnir_variable_bins_' + str(variable_bins_segments) + '_segments.txt'
        else:
            if swir_bool:
                my_data = my_data_swir
                file_name = 'txt_resp_phy_files/local_normalised_sequence_swir_variable_bins_' + str(variable_bins_segments) + '_segments.txt'
            else:
                my_data = my_data_vnir
                file_name = 'txt_resp_phy_files/local_normalised_sequence_vnir_variable_bins_' + str(variable_bins_segments) + '_segments.txt'
    else:
        if segment_wise_normalise_bool:
            if swir_bool:
                my_data = my_data_swir
                file_name = 'txt_resp_phy_files/unnormalised_sequence_segment_normalised_swir_variable_bins_' + str(variable_bins_segments) + '_segments.txt'
            else:
                my_data = my_data_vnir
                file_name = 'txt_resp_phy_files/unnormalised_sequence_segment_normalised_vnir_variable_bins_' + str(variable_bins_segments) + '_segments.txt'
        else:
            if swir_bool:
                my_data = my_data_swir
                file_name = 'txt_resp_phy_files/unnormalised_sequence_swir_variable_bins_' + str(variable_bins_segments) + '_segments.txt'
            else:
                my_data = my_data_vnir
                file_name = 'txt_resp_phy_files/unnormalised_sequence_vnir_variable_bins_' + str(variable_bins_segments) + '_segments.txt'#
    return my_data, file_name, colors
def write_line(line):
    text_file = open(file_name, "a")
    text_file.write(line)
    text_file.write('\n')
    text_file.close()
def write_lines(line_list):
    p = open("data/txt_resp_phy_files/specimen_phy_file_DO_NOT_EDIT (copy).txt", "r")
    p_string = p.readlines()[1:]
    write_line('    16   288')
    for i, item in enumerate(p_string):
        split_objects = item.split(' ')
        first_word = split_objects[0]
        length_of_FW = len(first_word)
        gap = 15 - length_of_FW
        new_string = first_word + line_list[i].rjust(gap + len(line_list[i]))
        write_line(new_string)
    return None
#inputs#####inputs#####inputs####inputs####inputs####inputs####inputs####inputs####inputs####inputs####
my_data_swir = genfromtxt('data/SWIR.AllData_avg - Copy.csv', delimiter=',')
my_data_vnir = genfromtxt('data/VNIR.AllData_avg - Copy.csv', delimiter=';')
#hyperparameters####hyperparameters####hyperparameters####hyperparameters####hyperparameters###########
save_text_bool = 1
normalise_bool= 1
swir_bool = 1 #TRUE for SWIR, FALSE for VNIR data
img_bool = 0
rgb_bool = 0
segment_wise_normalise_bool = 1
variable_bins_segments = 3 # no.of segments to apply variable binning in

my_data, file_name, colors = hyperparameters(swir_bool, normalise_bool, segment_wise_normalise_bool, rgb_bool,variable_bins_segments)
file_name = 'data/txt_resp_phy_files/trial_delete_this_3.phy'
#################################################################################################
# for i in range(variable_bins_segments):
#     my_data, file_name, colors = hyperparameters(swir_bool, normalise_bool, segment_wise_normalise_bool, rgb_bool,
#                                                  variable_bins_segments)
#
#     variable_bins_segments = i+1 # no.of segments to apply variable binning i
#     variable_bin_sequence_generator(my_data, variable_bins_segments)
#################################################################################################
print(np.asarray(variable_bin_sequence_generator(my_data, variable_bins_segments)))
#################################################################################################

