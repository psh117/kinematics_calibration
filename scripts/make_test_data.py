import csv
from copy import deepcopy
from math import pi
# q1,q2,q3,q4,q5,q6,q7,x,y,z
# read data

arm = ['panda_left','panda_right','panda_top']
set1 = [0, 0, 0, -pi/3, 0, pi/3, pi/4, 0, 0, 0, -pi/3, 0, pi/3, pi/4, 0, 0, 0, -pi/3, 0, pi/3, pi/4]
def rearrange_data_():
    f = open('./build/data_1.txt','w',newline='')
    for i in range(len(arm)):
        for j in range(len(set1)):
            f.write(str(set1[j]))
            if j < (len(set1) - 1):
                f.write('\t')
        if i < (len(arm) - 1):
            f.write('\n')
    f.close()

rearrange_data_()