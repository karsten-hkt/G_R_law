import numpy as np
import b_value_calculate as b_cal
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#读取数据
data_locate = '../data/zhou_data_processing.csv'
data = pd.read_csv(data_locate)
#将data专为numpy方便计算
data_numpy = np.array(data)
#计算总的完备震级
mag_line = 5
delta_M = 0.1
minnum = 50
b_all_area, b_error_all_area, Mc_all_area = b_cal.get_b_Mc(data_numpy,minnum,delta_M,mag_line)
print(r'b=%.2f+-%.3f，Mc=%f'%(b_all_area, b_error_all_area, Mc_all_area))
#按照每个月去计算完备震级并画图
fig = plt.figure(dpi=300,figsize=(8,5))
b_1, b_error_1, Mc_1 = b_cal.month_b_value_calculate_and_draw(data_numpy[(data_numpy[:,0]<6)],
                                                              minnum, delta_M, mag_line, color='red')
b_2, b_error_2, Mc_2 = b_cal.month_b_value_calculate_and_draw(data_numpy[(data_numpy[:,0]<12) & (data_numpy[:,0]>=6)],
                                                              minnum, delta_M, mag_line, color='orange')
b_3, b_error_3, Mc_3 = b_cal.month_b_value_calculate_and_draw(data_numpy[(data_numpy[:,0]>=12) & (data_numpy[:,0]<18)],
                                                              minnum, delta_M, mag_line, color='m')
b_4, b_error_4, Mc_4 = b_cal.month_b_value_calculate_and_draw(data_numpy[(data_numpy[:,0]>=18)],
                                                              minnum, delta_M, mag_line, color='g')
plt.yscale("log")
labels = ['from 02/06 to 02/11','from 02/12 to 02/17','from 02/18 to 02/23','from 02/23 to 02/28']
color = ['red', 'orange', 'm', 'g']#,'cyan','blue']
patches = [ mpatches.Patch(color=color[i], label="{:s}".format(labels[i]) ) for i in range(len(color))]
ax=plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*1.2, box.height*1.2])
ax.legend(handles=patches, bbox_to_anchor=(1,1), ncol=1,fontsize=10)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Mw',fontdict={'weight': 'normal', 'size': 15})
plt.ylabel('lg(N)',fontdict={'weight': 'normal', 'size': 12})
plt.savefig('../out_file/zhou_result/zhou_turkey_b_month.png',bbox_inches='tight')
plt.savefig('../out_file/zhou_result/zhou_turkey_b_month.pdf',bbox_inches='tight')