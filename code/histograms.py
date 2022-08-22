

import numpy as np
import pandas
import matplotlib.pyplot as plt

from numpy import genfromtxt

my_data_swir = genfromtxt('data/SWIR.AllData_avg - Copy.csv', delimiter=',')
my_data_vnir = genfromtxt('data/VNIR.AllData_avg - Copy.csv', delimiter=';')
my_data = my_data_swir

print(my_data[0])
print(my_data[1])
print(my_data[2])
print(len(my_data[0][2:]))
print(len(my_data[1][2:]))

wavelengths = my_data[0][2:]
max_wavelength = max(wavelengths)
min_wavelength = min(wavelengths)

m1_gelb     = my_data[1][3:]; m1_grau     = my_data[2][3:]; m1_rot      = my_data[3][3:]; m1_blau     = my_data[4][3:]

m2_gelb     = my_data[5][3:]; m2_rot     = my_data[6][3:]; m2_blau      = my_data[7][3:]; m2_grau     = my_data[8][3:]

m3_rot     = my_data[9][3:]; m3_gelb     = my_data[10][3:]; m3_blau      = my_data[11][3:]; m3_grau     = my_data[12][3:]

m4_gelb     = my_data[13][3:]; m4_grau     = my_data[14][3:]; m4_rot      = my_data[15][3:]; m4_blau     = my_data[16][3:]


plt.plot(m2_gelb - m3_gelb, color='y')
# plt.plot(m2_rot, color='r')
# plt.plot(m2_grau, color='gray')
# plt.plot(m2_blau, color='b')
#
# plt.plot(m3_gelb, color='y')
# plt.plot(m3_rot, color='r')
# plt.plot(m3_grau, color='gray')
# plt.plot(m3_blau, color='b')
plt.show()

# print(len(m4_gelb), len(m4_grau),
#       len(m4_rot), len(m4_blau))
# print(m4_gelb[0:3])
# print(m4_grau[0:3])
# print(m4_rot[0:3])
# print(m4_blau[0:3])

# fig, axes = plt.subplots(1, 3)
# axes[0].hist([m1_gelb], 20, color='y')
# axes[1].hist([m1_rot], 20, color='r')
# axes[2].hist(m1_gelb / m1_gelb, 20, color='g')
# plt.show()

# plt.autoscale(enable=True)
# denoms = [m1_gelb,m1_rot,m1_blau,m1_grau,m2_gelb,m2_rot,m2_blau,m2_grau,m3_gelb,m3_rot,m3_blau,m3_grau,m4_gelb,m4_rot,m4_blau,m4_grau]
# for i in range(16):
#     denom = denoms[i]
#     i=i+1
#     stringg = '/home/dheeraj/CRNT/ovgu/fraunhofer_work/01_2022/line_ratios/'+str(i)
#     plt.hist(m1_gelb/denom, 20,color='y'); plt.savefig(stringg+'_1.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m1_rot/denom, 20,color='r'); plt.savefig(stringg+'_2.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m1_blau/denom, 20,color='b');plt.savefig(stringg+'_3.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m1_grau/denom, 20,color='gray');plt.savefig(stringg+'_4.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m2_gelb/denom, 20,color='y'); plt.savefig(stringg+'_5.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m2_rot/denom, 20,color='r');plt.savefig(stringg+'_6.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m2_blau/denom, 20,color='b');plt.savefig(stringg+'_7.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m2_grau/denom, 20,color='gray');plt.savefig(stringg+'_8.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m3_gelb/denom, 20,color='y');plt.savefig(stringg+'_9.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m3_rot/denom, 20,color='r');plt.savefig(stringg+'_10.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m3_blau/denom, 20,color='b');plt.savefig(stringg+'_11.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m3_grau/denom, 20,color='gray');plt.savefig(stringg+'_12.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m4_gelb/denom, 20,color='y');plt.savefig(stringg+'_13.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m4_rot/denom, 20,color='r');plt.savefig(stringg+'_14.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m4_blau/denom, 20,color='b');plt.savefig(stringg+'_15.png');plt.show(block=False);plt.pause(0.1);plt.close()
#     plt.hist(m4_grau/denom, 20,color='gray');plt.savefig(stringg+'_16.png');plt.show(block=False);plt.pause(0.1);plt.close()
#
# # fig, axes = plt.subplot(1,2,2)
# n, bins, patches = plt.hist(m1_gelb, 20, density=True, facecolor='g', alpha=0.75)
# plt.hist(m1_gelb, 20, density=True, facecolor='green', alpha=0.75)
# plt.hist(m1_grau, 20, density=True, facecolor='gray', alpha=0.75)
# plt.hist(m1_rot, 20, density=True, facecolor='red', alpha=0.75)
# plt.hist(m1_blau, 20, density=True, facecolor='blue', alpha=0.75)
#
# plt.hist(m2_gelb, 20, density=True, facecolor='green', alpha=0.75)
# plt.hist(m2_grau, 20, density=True, facecolor='gray', alpha=0.75)
# plt.hist(m2_rot, 20, density=True, facecolor='red', alpha=0.75)
# plt.hist(m2_blau, 20, density=True, facecolor='blue', alpha=0.75)
#
# plt.hist(m3_gelb, 20, density=True, facecolor='green', alpha=0.75)
# plt.hist(m3_grau, 20, density=True, facecolor='gray', alpha=0.75)
# plt.hist(m3_rot, 20, density=True, facecolor='red', alpha=0.75)
# plt.hist(m3_blau, 20, density=True, facecolor='blue', alpha=0.75)
#
# plt.hist(m4_gelb, 20, density=True, facecolor='green', alpha=0.75)
# plt.hist(m4_grau, 20, density=True, facecolor='gray', alpha=0.75)
# plt.hist(m4_rot, 20, density=True, facecolor='red', alpha=0.75)
# plt.hist(m4_blau, 20, density=True, facecolor='blue', alpha=0.75)
# print(n)
# print(bins)
# print(patches[19])
# plt.xlabel('reflectance')
# plt.ylabel('frequency')
# plt.title('Histogram of material_4_hist')
# plt.title('Histogram of material_1_blau')
#
# plt.title('Histogram of material_2_gelb')
# plt.title('Histogram of material_2_grau')
# plt.title('Histogram of material_2_rot')
# plt.title('Histogram of material_2_blau')
#
# plt.title('Histogram of material_3_gelb')
# plt.title('Histogram of material_3_grau')
# plt.title('Histogram of material_3_rot')
# plt.title('Histogram of material_3_blau')
#
# plt.title('Histogram of material_4_gelb')
# plt.title('Histogram of material_4_grau')
# plt.title('Histogram of material_4_rot')
# plt.title('Histogram of material_4_blau')
# plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
# plt.xlim(40, 160)
# plt.ylim(0, 0.03)
# plt.grid(True)
# plt.show()