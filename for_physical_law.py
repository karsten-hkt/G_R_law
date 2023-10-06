import numpy as np
import b_value_calculate as b_cal
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#读取数据
data_locate = '../data/data_processing.csv'
data = pd.read_csv(data_locate)
#将data专为numpy方便计算
data_numpy = np.array(data)
#计算总的完备震级
mag_line = 5
delta_M = 0.1
minnum = 50
b_all_area, b_error_all_area, Mc_all_area = b_cal.get_b_Mc(data_numpy,minnum,delta_M,mag_line)
print(r'b=%.2f+-%.3f，Mc=%f'%(b_all_area, b_error_all_area, Mc_all_area))
fig = plt.figure(dpi=300,figsize=(8,5))
b_1, b_error_1, Mc_1 = b_cal.month_b_value_calculate_and_draw(data_numpy,
                                                              minnum, delta_M, mag_line, color='red')
GR_Mc_N_counts = b_cal.count_GR_Mc_N(data_numpy, Mc_all_area, delta_M, mag_line)
a = b_cal.a_calculate(Mc_all_area, b_all_area, GR_Mc_N_counts, delta_M)
print(a)
np.savetxt('../out_file/for_law/GR_counts.txt',GR_Mc_N_counts)
plt.savefig('../out_file/for_law/GR.png')