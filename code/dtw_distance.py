import numpy as np
from fastdtw import fastdtw
from numpy import genfromtxt

def normalise(array):
    norm = []
    max = np.max(array)
    min = np.min(array)
    # for i in range(len(array)):
    #     norm.append((array[i] - min) / (max- min))
    norm = (array - min) / (max - min)
    return np.asarray(norm)

def data_conversion(original_csv):
    my_data_swir = original_csv
    # wavelengths = my_data_swir[0][2:]
    m1_y = my_data_swir[1][2:];    # m1_y = np.stack((wavelengths, m1_y), axis=-1)
    m1_g = my_data_swir[2][2:];
    m1_r = my_data_swir[3][2:];
    m1_b = my_data_swir[4][2:];

    m2_y = my_data_swir[5][2:];
    m2_r = my_data_swir[6][2:];
    m2_b = my_data_swir[7][2:];
    m2_g = my_data_swir[8][2:];

    m3_r = my_data_swir[9][2:];
    m3_y = my_data_swir[10][2:];
    m3_b = my_data_swir[11][2:];
    m3_g = my_data_swir[12][2:];

    m4_y = my_data_swir[13][2:];
    m4_g = my_data_swir[14][2:];
    m4_r = my_data_swir[15][2:];
    m4_b = my_data_swir[16][2:];
    A = []
    mats = [m1_y,m1_g,m1_r,m1_b, m2_y,m2_r,m2_b,m2_g,m3_r,m3_y,m3_b,m3_g,m4_y,m4_g,m4_r,m4_b]
    for i, i_mat in enumerate(mats):
           for j, j_mat in enumerate(mats):
               if j < i :
                   A.append(0)
                   continue
               distance, path = fastdtw(i_mat, j_mat);
               A.append(distance);
    A = (A - min(A)) / (max(A) - min(A))
    A = np.reshape(A, -1,16)
    return A

# my_data_swir = genfromtxt('data/VNIR.AllData_avg - Copy.csv', delimiter=';')
my_data_swir = genfromtxt('data/SWIR.AllData_avg - Copy.csv', delimiter=',')
A = data_conversion(my_data_swir)

# print(A)
# print("length of lst is " + str(len(A)))
# print("normalised lst is ###################################")
A = (A - min(A)) / (max(A) - min(A))
# print(A)
B = np.reshape(A, (-1, 16))
print(B)
np.savetxt("DTW_distance_matrix_of_global_normalised_data_seg" + str() + ".csv", B, delimiter=",")

# np.savetxt("DTW_distance_matrix_of_normalised_data.csv", B, delimiter=",")

# lost = [[0,0,0,0],
#        [lst[0],0,0,0],
#        [lst[1],lst[3],0,0],
#        [lst[2],lst[4],lst[5],0]]
# print(np.asarray(lost))