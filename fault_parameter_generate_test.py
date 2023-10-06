import numpy as np
import pandas as pd
import location_trans as LocT
import area_cut
import b_value_calculate as b_cal
import gmt_draw

fault_name = '../data/turkey_7_8_fault.txt'
data_name = '../data/data_processing.csv'
length = 20000
fault = pd.read_csv(fault_name)
data = pd.read_csv(data_name)
data = np.array(data)
#build fault
fault_para = LocT.fault_parameter_genereate(fault.iloc[0,1],fault.iloc[0,2],fault.iloc[0,3],fault.iloc[0,5])
fault_para = np.array(fault_para)
#build buffer
buffer_build = LocT.buffer_generate(fault_para,length)
buffer_build = np.array(buffer_build)
#select points in the buffer
out_aftershocks = []
for i in range(data.shape[0]):
    point = np.zeros(2)
    point[0] = data[i, 2]
    point[1] = data[i, 3]
    if LocT.is_in_poly(point, buffer_build):
        out_aftershocks.append(data[i, :])
out_aftershocks = np.array(out_aftershocks)
#参数设置
para = {'day_line':0, 'lon_line':2, 'lat_line':3, 'dep_line':4, 'mag_line': 5, 'delta_M':0.1, 'roi':[34,41,35,40,-20,0],
        'delta_time':10, 'spacing':2, 'buf':3, 'calculate_min_number':50, 'drawing_number':5, 'method':'maxc',
        'save_file':False}
#将aftershocks投影到所在断层面内
outx, outy = LocT.simp_rotaxll(out_aftershocks[:, 2], out_aftershocks[:, 3],  fault.iloc[0,5], fault.iloc[0,1],
                               fault.iloc[0,2])
out_aftershocks[:, 2] = outx  # 距离震中距离
#生成矩形区域并计算b值
ext = [0,fault.iloc[0,3],0,20]
#创建矩形
refpts, refpolys = area_cut.ext2refpts(ext,para['spacing'])
refpts_poly =  area_cut.refpts2bufferpoly(refpts,para['spacing'],buf=para['buf'])
# 设置b值保存位置
b_value = np.zeros(refpts.shape[0])
b_error = np.zeros(refpts.shape[0])
## 计算全部的完备震级
tmp_time = 1
m_total = b_cal.get_maxc(out_aftershocks[:, para['mag_line']], para['delta_M'])
for i in range(refpts.shape[0]):
    b_value[i], b_error[i] = b_cal.eqs2bINpoly_maxc(out_aftershocks, refpts_poly[i],
                                                    para['calculate_min_number'], para['delta_M'],
                                                    para['lon_line'], para['mag_line'], para['dep_line'],
                                                    refpolys[i], tmp_time)
#将计算结果转成成实际坐标并保存成gmt格式
outgmt_name = '3d_Turkey_2023_cross_section_%s_'%para['delta_time']
out_gmt_locate = '../out_file/3d_cross_section/'
area_cut.ref2gmt_3d(refpolys, b_value, fault.iloc[0,5], fault.iloc[0,1],
                               fault.iloc[0,2], outgmt_name, tmp_time, out_gmt_locate)
gmt_draw.gmt_draw_cross_section_3d(para['roi'], outgmt_name, tmp_time, out_gmt_locate)

