import numpy as np
import pandas as pd
import location_trans as LocT
import area_cut
import b_value_calculate as b_cal
import gmt_draw
import time_cut

#读取数据
fault_name = '../data/fault_all.txt'
#这个断层的数据有所区别
data_name = '../data/data_processing_7_8.csv'
fault = pd.read_csv(fault_name)
data = pd.read_csv(data_name)
data = np.array(data)
#基础参数设置
para = {'day_line':0, 'lon_line':2, 'lat_line':3, 'dep_line':4, 'mag_line': 5, 'delta_M':0.1,
        'roi':[35.5,39,36.5,38.5,-20,0],'delta_time':10, 'spacing':2, 'buf':5, 'calculate_min_number':50,
        'drawing_number':5, 'buffer_length': 30000, 'depth':20,'projected_points_save':False}
#主震
epi_1 = [37.22,38.08,-12]
#projection
p = [160,80]
#按照时间分割
for tmp_time in range(para['drawing_number']):
    print('开始第%s个时间间隔(间隔为%s)'%(tmp_time,para['delta_time']))
    data_now = time_cut.split_time(data,tmp_time,para['delta_time'],para['day_line'])
    # 将计算结果转成成实际坐标并保存成gmt格式
    outgmt_name = '3d_Turkey_2023_cross_section_%s_' % para['delta_time']
    out_gmt_locate = '../out_file/3d_all/'
    #按顺序读取断层
    for tmp_fault in range(len(fault)):
        fault_para = LocT.fault_parameter_genereate(fault.iloc[tmp_fault, 1], fault.iloc[tmp_fault, 2],
                                                    fault.iloc[tmp_fault, 3], fault.iloc[tmp_fault, 5])
        fault_para = np.array(fault_para)
        fault_point = [fault.iloc[tmp_fault, 1], fault.iloc[tmp_fault, 2], 0]
        # build buffer
        buffer_build = LocT.buffer_generate(fault_para, para['buffer_length'])
        buffer_build = np.array(buffer_build)
        # select points in the buffer
        out_aftershocks = []
        for i in range(data_now.shape[0]):
            point = np.zeros(2)
            point[0] = data_now[i, 2]
            point[1] = data_now[i, 3]
            if LocT.is_in_poly(point, buffer_build):
                out_aftershocks.append(data_now[i, :])
        out_aftershocks = np.array(out_aftershocks)
        out_aftershocks[:,para['dep_line']] = -out_aftershocks[:,para['dep_line']] #需要转为负值使用
        # 将aftershocks投影到所在断层面内
        projected_points, strike_x, strike_y = LocT.dot_pro2plane(out_aftershocks[:,[2,3,4]], fault_point, fault.iloc[tmp_fault, 5],fault.iloc[tmp_fault, 6])
        if para['projected_points_save']:
            np.savetxt('projected_points.txt',projected_points)
        out_aftershocks[:, 2] = np.array(strike_x)/1000 # lon转换为x，转为km
        out_aftershocks[:, 3] = np.array(strike_y)/1000 # lat转换为y，转为km
        # 生成矩形区域并计算b值
        ext = [0, fault.iloc[tmp_fault, 3], 0, para['depth']]
        # 创建矩形
        refpts, refpolys = area_cut.ext2refpts(ext, para['spacing'])
        refpts_poly = area_cut.refpts2bufferpoly(refpts, para['spacing'], buf=para['buf'])
        # 设置b值保存位置
        b_value = np.zeros(refpts.shape[0])
        b_error = np.zeros(refpts.shape[0])
        ## 计算全部的完备震级
        m_total = b_cal.get_maxc(out_aftershocks[:, para['mag_line']], para['delta_M'])
        for i in range(refpts.shape[0]):
            b_value[i], b_error[i] = b_cal.eqs2bINpoly_maxc(out_aftershocks, refpts_poly[i],
                                                            para['calculate_min_number'], para['delta_M'],
                                                            2, para['mag_line'], 3,
                                                            refpolys[i], tmp_time) #这时候dep代表为在dip方向的距离
        if tmp_fault == 0:
            area_cut.ref2gmt_3d_dip(refpolys, b_value, fault.iloc[tmp_fault, 5], fault.iloc[tmp_fault, 6], fault.iloc[tmp_fault, 1],
                                fault.iloc[tmp_fault, 2], outgmt_name, tmp_time, out_gmt_locate)
        else:
            area_cut.ref2gmt_3d_dip_exist(refpolys, b_value, fault.iloc[tmp_fault, 5], fault.iloc[tmp_fault, 6], fault.iloc[tmp_fault, 1],
                                fault.iloc[tmp_fault, 2], outgmt_name, tmp_time, out_gmt_locate)
    gmt_draw.gmt_draw_cross_section_3d(para['roi'], epi_1, outgmt_name, tmp_time, out_gmt_locate, p)
