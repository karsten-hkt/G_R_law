import pandas as pd
import numpy as np
import datetime
from datetime import timedelta

def time_cut(data):
    #cut the time based on the style of 2023-05-31T23:57:18
    day = data.iloc[:, 0].str.split('T').str[0]
    time = data.iloc[:, 0].str.split('T').str[1]
    data_new = pd.DataFrame(columns=['year','month','day','min','sec','lon','lat','depth','mag'])
    data_new['year'] = day.str.split('-').str[0].astype('int')
    data_new['month'] = day.str.split('-').str[1].astype('int')
    data_new['day'] = day.str.split('-').str[2].astype('int')
    data_new['min'] = time.str.split(':').str[0].astype('int')
    data_new['sec'] = time.str.split(':').str[1].astype('int')
    data_new['lon'] = data.iloc[:, 1]
    data_new['lat'] = data.iloc[:, 2]
    data_new['depth'] = data.iloc[:, 3]
    data_new['mag'] = data.iloc[:, 4]
    data_numpy = np.array(data_new)
    return data_numpy

def time_delta(time, start_time, datetimeFormat):
    #传入时间行
    #最终将输出时间间隔，以天为单位
    #默认格式为 '2019-07-06T03:19:52.260'
    diff_days = []
    for time_str in time.astype(str):
        diff = datetime.datetime.strptime(time_str, datetimeFormat) - datetime.datetime.strptime(start_time, datetimeFormat)
        diff_days.append(diff.days + diff.seconds / 24 / 60 / 60)
    ## 将获取到的时间间隔扩展到SEEQ中
    return diff_days


def split_time(aftershocks_file,tmp_time,delta_time,day_line):
    day_down = tmp_time*delta_time
    day_up = (tmp_time+1)*delta_time
    aftershocks_file = aftershocks_file[aftershocks_file[:,day_line]<day_up,:]
    aftershocks_file = aftershocks_file[aftershocks_file[:,day_line]>=day_down,:]
    print("此时距离起始时间：",day_down)
    print("此时距离截至时间：",day_up)
    return aftershocks_file


if __name__ == '__main__':
    print('This is for time trans')

