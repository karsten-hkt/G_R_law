import numpy as np
import pandas as pd
import area_cut
import time_cut
import location_trans
import b_value_calculate as b_cal
import gmt_draw

#读取数据
data_locate = '../out_file/zhou_result/zhou_aftershocks_in_buffer_0.csv'
data = pd.read_csv(data_locate)
#将data专为numpy方便计算
data_numpy = np.array(data)
#参数设置
para = {'day_line':0, 'lon_line':3, 'lat_line':2, 'dep_line':4, 'mag_line': 5, 'delta_M':0.1, 'roi':[34,41,35,40],
        'delta_time':30, 'spacing':2, 'buf':4, 'calculate_min_number':50,'ext':[-100,170,0,20], 'strike':66.6,
        'drawing_number':1, 'method':'maxc','save_file':False}
#ext 剖面图数据,#method = 'maxc'#完备震级计算方法,#drawing_number = 5#需要画多少个时间间隙,#save_file = False#是否保存相关文件
## 输入主震震中
meca_file = '../data/meca_turkey_6.dat'
epi_1 = [37.45,37.55,15]
epi_locate = '../out_file/zhou_result/epi.txt'
outgmt_name = 'zhou_Turkey_2023_cross_section_%s_'%para['delta_time']
out_gmt_locate = '../out_file/zhou_result/'
#other epi:
epi_2 = np.array([[37.45,37.55,15],
                  [36.84,37.21,26],
                  [37.64,37.73,21],
                 [37.63, 37.92,21]])
epi_2_now_x, epi_2_now_y = location_trans.simp_rotaxll(epi_2[:, 0], epi_2[:, 1], para['strike'], epi_1[0], epi_1[1])
epi_2[:, 0] = epi_2_now_x
epi_2[:, 1] = epi_2_now_y
print(epi_2)
np.savetxt(epi_locate, epi_2[:, [0, 2]], delimiter=' ', fmt='%.5f')
for tmp_time in range(para['drawing_number']):
    print('开始第%s个时间间隔(间隔为%s)'%(tmp_time,para['delta_time']))
    data_now = time_cut.split_time(data_numpy,tmp_time,para['delta_time'],para['day_line'])
    outx,outy = location_trans.simp_rotaxll(data_now[:,para['lon_line']],
                                            data_now[:,para['lat_line']],para['strike'],epi_1[0],epi_1[1])
    data_now[:,para['lon_line']] = outx #距离震中距离
    data_now[:,para['lat_line']] = outy
    ## 记录M>3的余震
    m_up_3_locate = '../out_file/cross_section/m_up_3.txt'
    M_up_3_data = data_now[data_now[:,para['mag_line']]>3,:]
    ## 剔除掉7.1级的地震
    M_up_3_data[:,para['mag_line']] = (M_up_3_data[:,para['mag_line']] /
                                       100*1.25**(M_up_3_data[:,para['mag_line']])) #为了之后gmt画图余震的大小
    np.savetxt(m_up_3_locate,M_up_3_data[:,[para['lon_line'],
                                            para['dep_line'],para['mag_line']]],delimiter=' ',fmt='%.5f')
    ## 创建矩形区域
    refpts, refpolys = area_cut.ext2refpts(para['ext'],para['spacing'])
    refpts_poly =  area_cut.refpts2bufferpoly(refpts,para['spacing'],buf=para['buf'])
    #print(refpts_poly )
    ## 创建记录b值的数组
    b_value = np.zeros(refpts.shape[0])
    b_error = np.zeros(refpts.shape[0])
    ## 计算全部的完备震级
    m_total = b_cal.get_maxc(data_now[:,para['mag_line']],para['delta_M'])
    for i in range(refpts.shape[0]):
        if para['method'] == 'fmbass':
            b_value[i],b_error[i] = b_cal.eqs2bINpoly_fmbass(data_now,refpts_poly[i],
                                                             para['calculate_min_number'],para['delta_M'],
                                                             para['lon_line'],para['mag_line'],para['dep_line'],
                                                             m_total)
        elif para['method']  == 'maxc':
            b_value[i],b_error[i] = b_cal.eqs2bINpoly_maxc(data_now,refpts_poly[i],
                                                           para['calculate_min_number'],para['delta_M'],
                                                           para['lon_line'],para['mag_line'],para['dep_line'],
                                                           refpolys[i],tmp_time)
        else:
            print('method wrong')

    ## 转成gmt格式
    area_cut.ref2gmt(refpolys,b_value, outgmt_name, tmp_time, out_gmt_locate)
    #画图
    gmt_draw.gmt_draw_cross_section(para['ext'], m_up_3_locate, epi_1, epi_locate, outgmt_name, tmp_time, out_gmt_locate)

