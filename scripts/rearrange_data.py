import csv
from copy import deepcopy
# q1,q2,q3,q4,q5,q6,q7,x,y,z
# read data

arm = ['panda_left','panda_right','panda_top']

def rearrange_data_(arm_name):
    data_out = []
    q_info = []
    pos_info = []

    ## read pos
    f = open('./data/KETI/'+arm_name+'_no_q.csv', 'r', encoding='utf-8')
    reader = csv.reader(f)
    for line in reader:
        pos_info.append(line)
    f.close()

    ## read q
    f = open('./data/SNU/calib_pos.csv', 'r', encoding='utf-8')
    reader = csv.reader(f)
    for line in reader:
        q_info.append(line[1:])
    f.close()

    for i in range(len(pos_info)):
        new_line = [str(0.0)]
        for j in range(7):
            new_line.append(q_info[int(pos_info[i][0])-1][j])
        new_line.append(str(0.0))
        for j in range(3):
            new_line.append(str(float(pos_info[i][j+1])*0.001))
        data_out.append(new_line)

    # f = open('./src/'+arm_name+'.txt','w',encoding='utf-8',newline='')
    f = open('./build/'+arm_name+'.txt','w',newline='')
    for i in range(len(data_out)):
        for j in range(len(data_out[i])):
            f.write(data_out[i][j])
            if j < (len(data_out[i])-1):
                f.write('\t')
        if i < (len(data_out)-1):
            f.write('\n')
    f.close()

for arm in arm:
    rearrange_data_(arm)