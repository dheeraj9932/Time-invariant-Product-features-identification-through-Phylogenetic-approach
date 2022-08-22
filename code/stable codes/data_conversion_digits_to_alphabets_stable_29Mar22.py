import numpy as np
import pandas
import string
import math
import matplotlib.pyplot as plt
from numpy import genfromtxt


def normalise(array):
    norm = []
    max = np.max(array)
    min = np.min(array)
    # for i in range(len(array)):
    #     norm.append((array[i] - min) / (max- min))
    norm = (array - min) / (max - min)
    return np.asarray(norm)
def string_length(len):
    lower_case = string.ascii_lowercase
    upper_case = string.ascii_uppercase
    total = list(upper_case + lower_case)
    if len <= 52:
        total = total[0:len]
        return total
    else:
        return total
    # indexes_string = ''
    # counter = 0
    # index_list = []
    # for i in range(26):
    #     k = upper_case[i]
    #     for j in range(26):
    #         if counter == len:
    #             return index_list
    #         index_list += [k+lower_case[j]]
    #         indexes_string = indexes_string+k+lower_case[j]
    #         counter += 1
    # return total
def array_2_string_automated(sequence, save_text_bool, file_name, partitions):
    protein_encoding = string_length(partitions)
    # all_ascii_chars = "".join(chr(x) for x in range(33, 127))
    # total_no_of_steps = len(all_ascii_chars)
    # protein_encoding = all_ascii_chars
    steps = np.linspace(0, 1, num=partitions)
    protein_seq = ''
    aug_seq = (sequence * 10000).astype(int)
    aug_steps = (steps*10000).astype(int)
    print(aug_seq[0:20])
    print(aug_steps)
    print(protein_encoding)
    indexes = []
    for i in range(len(aug_seq)):
        protein_seq = protein_seq + protein_encoding[np.where(aug_steps <= (aug_seq[i]))[0][-1]] #adding the characters corresponding to the value of normalised reflectance data to 'protein seq'
    if save_text_bool:
        text_file = open("data/"+file_name, "a")
        text_file.write(protein_seq)
        text_file.write('\n')
        text_file.close()
    print((protein_seq))
    print('\n')
    return protein_seq
def data_conversion(original_csv, save_text_bool, normalise_bool,file_name,partitions):
    my_data = original_csv
    # wavelengths = my_data[0][2:]
    m1_y = my_data[1][2:];    # m1_y = np.stack((wavelengths, m1_y), axis=-1)
    m1_g = my_data[2][2:];
    m1_r = my_data[3][2:];
    m1_b = my_data[4][2:];

    m2_y = my_data[5][2:];
    m2_r = my_data[6][2:];
    m2_b = my_data[7][2:];
    m2_g = my_data[8][2:];

    m3_r = my_data[9][2:];
    m3_y = my_data[10][2:];
    m3_b = my_data[11][2:];
    m3_g = my_data[12][2:];

    m4_y = my_data[13][2:];
    m4_g = my_data[14][2:];
    m4_r = my_data[15][2:];
    m4_b = my_data[16][2:];

    mats = [m1_y,m1_g,m1_r,m1_b, m2_y,m2_r,m2_b,m2_g,m3_r,m3_y,m3_b,m3_g,m4_y,m4_g,m4_r,m4_b]
    for i in range(len(mats)):
        if normalise_bool:
            array_2_string_automated(normalise(mats[i]), save_text_bool, file_name, partitions)
        else:
            array_2_string_automated(mats[i], save_text_bool, file_name, partitions)

#inputs
my_data_swir = genfromtxt('data/SWIR.AllData_avg - Copy.csv', delimiter=',')
my_data_vnir = genfromtxt('data/VNIR.AllData_avg - Copy.csv', delimiter=';')
my_data = my_data_swir


#hyperparameters
save_text_bool = 1
normalise_bool= 1
file_name = 'txt_resp_phy_files/normalised_sequence_swir_53.txt'
partitions = 53

#execution
data_conversion(my_data, save_text_bool, normalise_bool, file_name, partitions)

#where i left off: creating a universal sequence converter. edit this python file.

