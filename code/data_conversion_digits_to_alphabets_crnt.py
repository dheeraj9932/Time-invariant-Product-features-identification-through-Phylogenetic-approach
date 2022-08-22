import numpy as np
import pandas
import string
import random
import math
import matplotlib.pyplot as plt
from numpy import genfromtxt

def graph_plot(normalise_bool):
    wavelengths = my_data_vnir[0][2:]
    m1_y = my_data[1][2:];  # m1_y = np.stack((wavelengths, m1_y), axis=-1)
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

    mats = [m1_y, m1_g, m1_r, m1_b, m2_y, m2_r, m2_b, m2_g, m3_r, m3_y, m3_b, m3_g, m4_y, m4_g, m4_r, m4_b]
    mats_string = ['m1_y', 'm1_g', 'm1_r', 'm1_b','m2_y', 'm2_r', 'm2_b', 'm2_g', 'm3_r', 'm3_y', 'm3_b', 'm3_g','m4_y', 'm4_g', 'm4_r', 'm4_b']
    global_max = 0
    global_min = 0

    colors = "bgrcmykw"
    colors = ['yo','ko','ro','bo','yx','rx','bx','kx',
              'r+','y+','b+','k+','y-','k-','r-','b-']
    for i in range(len(mats)):
        if global_max < max(mats[i]):
            global_max = max(mats[i])
        if global_min > min(mats[i]):
            global_min = min(mats[i])
    if normalise_bool:
        for i, materials in enumerate(mats):
            # plt.plot(normalise(materials, global_max, global_min), c=colors[random.randint(0, 7)], label=mats_string[i])
            plt.plot(wavelengths, normalise(materials, global_max, global_min), colors[i], label=mats_string[i])
            # plt.grid()
    else:
        for i, materials in enumerate(mats):
            # plt.plot((materials), c=colors[random.randint(0, 7)], label=mats_string[i])
            plt.plot(wavelengths, (materials), colors[i], label=mats_string[i])
            # plt.grid()
    plt.legend(loc='best')
    plt.grid()
    plt.show()
def normalise(array, global_max,global_min):
    max = global_max
    min = global_min
    max = np.max(array)
    min = np.min(array)
    denominator = max - min
    norm = (array - min) / denominator
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
        return total
    else:
        return total
def array_2_string_automated(sequence, save_text_bool, file_name, partitions):
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
    graph_plot(normalise_bool)
    global_max = 0
    global_min = 0
    for i in range(len(mats)):
        if global_max < max(mats[i]):
            global_max = max(mats[i])
        if global_min > min(mats[i]):
            global_min = min(mats[i])
    for i in range(len(mats)):
        if normalise_bool:
            array_2_string_automated(normalise(mats[i], global_max,global_min), save_text_bool, file_name, partitions)
        else:
            array_2_string_automated(mats[i], save_text_bool, file_name, partitions)
#inputs##########################################################################################
my_data_swir = genfromtxt('data/SWIR.AllData_avg - Copy.csv', delimiter=',')
my_data_vnir = genfromtxt('data/VNIR.AllData_avg - Copy.csv', delimiter=';')
my_data = my_data_vnir
#hyperparameters##########################################################################################
save_text_bool = 0
normalise_bool= 0
partitions = 10
# characters = [10, 20,26, 30, 40, 52,53,62,70,80,93]
characters = [52]
for i, chars in enumerate(characters):
    file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_'+str(chars)+'.txt'
    #execution
    data_conversion(my_data, save_text_bool, normalise_bool, file_name, chars)
# graph_plot()