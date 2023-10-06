import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data_locate = '../data/data_processing.csv'
data = pd.read_csv(data_locate)
#将data专为numpy方便计算
data_numpy = np.array(data)

mag_line = 5
day_info = []
day_len = int(np.ceil(data_numpy[-1,0]))

data_af = np.zeros([len(data_numpy),2])
for i in range(len(data_numpy)):
    day = int(np.floor(data_numpy[i,0]))
    data_af[i,0] = day
    data_af[i,1] = data_numpy[i,mag_line]

for day in range(day_len):
    af_one_day = data_numpy[(data_numpy[:,0]>day) & (data_numpy[:,0]<=day+1)]
    day_count = []
    day_count.append(day)
    day_count.append(len(af_one_day))
    day_count.append(max(af_one_day[:,mag_line]))
    day_info.append(day_count)

day_info = np.array(day_info)
fig = plt.figure(figsize=(15,8),dpi=300)
ax1 = fig.add_subplot(111)
x = range(day_len)
lns1 = ax1.scatter(data_af[:,0],data_af[:,1], marker='o', c='r', alpha=0.2,s=2,label='Magnitude')
lns2, = ax1.plot(day_info[:,0],day_info[:,2],c='r',label='Magnitude Max')
ax1.set_ylabel('Magniude(M)',fontdict={'weight':'normal','size':15,'color':'red'})
ax1.set_xlabel('Postseismic accumulation time(Day)', fontdict={'weight':'normal','size':15,'color':'black'})
ax1.set_title('af_after_Mw7.6',fontsize=18)

ax2 = ax1.twinx()
lns3 = ax2.scatter(day_info[:,0],day_info[:,1],marker='o',c='b',s=6,label='Aftershock count')
ax2.set_ylabel('Aftershock count(Frequency)',fontdict={'weight':'normal','size':15,'color':'blue'})
ax2.set_xlabel('Same')

lns = [lns1,lns2,lns3]
labs = [l.get_label() for l in lns]
plt.legend(lns,labs,loc=0,fontsize=12)
plt.savefig('../out_file/magnitude_time/mt_plot.png',bbox_inches='tight')


