import numpy as np
import pandas as pd
import gmt_draw
import area_cut
import time_cut
import b_value_calculate as b_cal

#读取数据
data_locate = '../data/zhou_data_processing.csv'
data = pd.read_csv(data_locate)
#将data专为numpy方便计算
data_numpy = np.array(data)
## 设置参数
#roi: #lon_min,lon_max,lat_min,lat_max
#res_ratio: 分辨率
#distance：圆半径
#calculate_min_number:最少需要几个地震才去计算
para = {'day_line':0, 'lon_line':3, 'lat_line':2, 'mag_line': 5, 'delta_M':0.1, 'roi':[34,41,35,40],
        'delta_time':30, 'res_ratio':0.1, 'distance':0.2, 'calculate_min_number':50}
save_file = False#是否保存相关文件
## 输入主震震中
meca_file = '../data/meca_turkey_6.dat'
faults_file_name = '../data/fault.dat'
outgmt_name = 'zhou_Turkey_2023_%s_'%para['delta_time']
out_gmt_locate = '../out_file/zhou_result/'
X_len = round((para['roi'][1]-para['roi'][0])/para['res_ratio'])
Y_len = round((para['roi'][3]-para['roi'][2])/para['res_ratio'])
time_circle_number = np.ceil(max(data_numpy[:,para['day_line']])/para['delta_time'])
#fig = plt.figure(dpi=300,figsize=(8,5))
for tmp_time in range(round(time_circle_number)):
    print('开始第%s个时间间隔(间隔为%s)'%(tmp_time,para['delta_time']))
    data_now = time_cut.split_time(data_numpy,tmp_time,para['delta_time'],para['day_line'])
    ## 设置一些记录参数的数组
    Mc_record = np.zeros((X_len,Y_len))
    p_record = np.zeros((X_len,Y_len))
    b_record = np.zeros((X_len,Y_len))
    b_error_record = np.zeros((X_len,Y_len))
    b_record_num = np.zeros(X_len*Y_len)
    poly_record = np.zeros((X_len*Y_len,5,2))
    ## 开始计算每个矩形内部的参数
    for i in range(X_len):
        for j in range(Y_len):
            ## 当前所在矩形编号
            num = int(i*Y_len+j)
            ## 当前矩形位置roi: #lon_min,lon_max,lat_min,lat_max
            lon_min_now = para['roi'][0] + i*para['res_ratio']
            lon_max_now = para['roi'][0] + (i+1)*para['res_ratio']
            lat_min_now = para['roi'][2] + j*para['res_ratio']
            lat_max_now = para['roi'][2] + (j+1)*para['res_ratio']
            poly_record[num,:,:] = area_cut.record_polys(lon_min_now,lon_max_now,lat_min_now,lat_max_now)
            ## 周围一圈所包含的余震序列
            use_file = area_cut.circle_file(lon_min_now,lon_max_now,lat_min_now,lat_max_now,
                                            para['lon_line'],para['lat_line'],para['distance'],data_now)
            ## 需要内部余震数量大于一定值才计算
            if use_file.shape[0] >= para['calculate_min_number']:
                b_calculate_file = use_file
                ## 选择内部余震数量大于一定值才计算
                if b_calculate_file.shape[0] >= para['calculate_min_number']:
                    b_final_,b_error_record_,Mc = b_cal.get_b_Mc(b_calculate_file,para['calculate_min_number'],
                                                                 para['delta_M'],para['mag_line'])
                    b_record[i,j] = round(b_final_,5)
                    b_error_record[i,j] = round(b_error_record_,5)
                    Mc_record[i,j] = Mc
                else:
                    b_record[i,j] = np.NAN
                    b_error_record[i,j] = np.NAN
            else:
                Mc_record[i,j] = np.NAN
                p_record[i,j] = np.NAN
                b_record[i,j] = np.NAN
                b_error_record[i,j] = np.NAN
            b_record_num[num] = b_record[i,j]
    ## 是否保存需要的计算过程中的文件
    if save_file == True:
        np.savecsv('../data/b_spatial_temproal_variation/zhou_b_record_%s.txt'%tmp_time,
                   b_record, fmt='%.4f',newline='\n')
        np.savetxt('../data/b_spatial_temproal_variation/zhou_Mc_record_%s.txt'%tmp_time,
                   Mc_record, fmt='%.4f',newline='\n')
        np.savetxt('../data/b_spatial_temproal_variation/zhou_p_record_%s.txt'%tmp_time,
                   p_record, fmt='%.4f',newline='\n')
        np.savetxt('../data/b_spatial_temproal_variation/zhou_b_error_record_%s.txt'%tmp_time,
                   b_error_record, fmt='%.4f',newline='\n')
    ## 转换成gmt可读取的形式
    area_cut.ref2gmt(poly_record,b_record_num,outgmt_name,tmp_time,out_gmt_locate)
    print('结束第%s个时间间隔(间隔为%s)'%(tmp_time,para['delta_time']))

for tmp_time in range(round(time_circle_number)):
    gmt_draw.gmt_drawing(para['roi'], faults_file_name, meca_file, outgmt_name, tmp_time, out_gmt_locate)