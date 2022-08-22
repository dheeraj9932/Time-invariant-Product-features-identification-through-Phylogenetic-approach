import numpy as np
import similaritymeasures
import matplotlib.pyplot as plt
from numpy import genfromtxt


my_data_swir = genfromtxt('data/SWIR.AllData_avg - Copy.csv', delimiter=',')
W = my_data_swir[0][2:-1]
m1_y = my_data_swir[1][2:-1]; m1_y = np.stack((W, m1_y), axis=-1)
m1_g = my_data_swir[2][2:-1]; m1_g= np.stack((W, m1_g), axis=-1)
m1_r = my_data_swir[3][2:-1]; m1_r=np.stack((W, m1_r), axis=-1)
m1_b = my_data_swir[4][2:-1]; m1_b=np.stack((W, m1_b), axis=-1)
# m2_y = my_data_swir[5][2:-1]; np.stack((W, m2_y), axis=-1)
# m2_r = my_data_swir[6][2:-1]; np.stack((W, m2_r), axis=-1)
# m2_b = my_data_swir[7][2:-1]; np.stack((W, m2_b), axis=-1)
# m2_g = my_data_swir[8][2:-1]; np.stack((W, m2_g), axis=-1)
#
# m3_r = my_data_swir[9][2:-1]; np.stack((W, m3_r), axis=-1)
# m3_y = my_data_swir[10][2:-1]; np.stack((W, m3_y), axis=-1)
# m3_b = my_data_swir[11][2:-1]; np.stack((W, m3_b), axis=-1)
# m3_g = my_data_swir[12][2:-1]; np.stack((W, m3_g), axis=-1)
#
# m4_y = my_data_swir[13][2:-1]; np.stack((W, m4_y), axis=-1)
# m4_g = my_data_swir[14][2:-1]; np.stack((W, m4_g), axis=-1)
# m4_r = my_data_swir[15][2:-1]; np.stack((W, m4_r), axis=-1)
# m4_b = my_data_swir[16][2:-1]; np.stack((W, m4_b), axis=-1)

A = []
B = []
mats = [m1_y,m1_g,m1_r,m1_b]
# mats = [m1_y,m1_g,m1_r,m1_b,   m2_y,m2_r,m2_b,m2_g,   m3_r,m3_y,m3_b,m3_g,   m4_y,m4_g,m4_r,m4_b]
for i in mats:
       for j in mats:
           distance, path = similaritymeasures.dtw(np.round(i,2), np.round(j,2));A.append(distance)
           try:
               distance, path = similaritymeasures.frechet_dist(np.round(i, 2), np.round(j, 2));
               B.append(distance)
           except TypeError:
               B.append(0)

# Generate random experimental data
# x = np.sort(np.random.randint(0,100,100))
# y = np.sort(np.random.randint(0,100,100))
# exp_data = np.zeros((100, 2))
# exp_data[:, 0] = x
# exp_data[:, 1] = y
#
# # Generate random numerical data
# # x = np.random
# x = np.sort(np.random.randint(0,100,100))
# y = np.sort(np.random.randint(0,100,100))
# num_data = np.zeros((100, 2))
# num_data[:, 0] = x
# num_data[:, 1] = y



# Discrete Frechet distance
# df = similaritymeasures.frechet_dist(exp_data, num_data)
#
# # Dynamic Time Warping distance
# dtw, d = similaritymeasures.dtw(exp_data, num_data)

# print the results

print(len(A))
print(len(B))
print("dtw", A)
print("df", B)
# print(pcm, df, area, cl, dtw)

# plot the data
# plt.figure()
# plt.plot(exp_data[:, 0], exp_data[:, 1])
# plt.plot(num_data[:, 0], num_data[:, 1])
# plt.show()