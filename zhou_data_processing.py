import pandas as pd
import time_cut
import area_cut

#读取数据
data_locate = '../data/ding-zhou_eqs-2023_tk-palm_v4.ctlg'
data = pd.read_csv(data_locate)
data.columns = ['Date','Latitude','Longitude','Depth','Magnitude']
#计算时间间隔
start_time_7_8 = '2023-02-06T04:17:00.00Z'
start_time_7_6 = '2023-02-06T13:24:00.00Z'
datetimeFormat = '%Y-%m-%dT%H:%M:%S.%fZ'
diff_days = time_cut.time_delta(data.iloc[:,0], start_time_7_6, datetimeFormat)
#将时间间隔插入data中
data.insert(0,'Day_interval',diff_days)

#筛选出震后的数据，深度大于0的
data = data[data['Day_interval']>=0]
data = data[data['Depth']>0]
#筛选出在特定范围内的，研究区域设置
#研究区域34/41/35/40
roi = [34,41,35,40] #lon_min,lon_max,lat_min,lat_max
data = area_cut.in_ploy(data,roi)
#对时间间隔进行升序
data = data.sort_values(by=['Day_interval'],ascending=True)
data = data.reset_index(drop=True)
#保存成新的文件
data.to_csv('../data/zhou_data_processing.csv',index=0)