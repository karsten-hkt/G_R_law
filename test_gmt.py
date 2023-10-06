import gmt_draw

para = {'day_line':0, 'lon_line':2, 'lat_line':3, 'mag_line': 5, 'delta_M':0.1, 'roi':[34,41,35,40],
        'delta_time':180, 'res_ratio':0.1, 'distance':0.2, 'calculate_min_number':50}
save_file = False#是否保存相关文件
## 输入主震震中
epi_1 = [37.45, 37.55]#Mw7.8
faults_file_name = '../data/fault.dat'
epi_2 = [36.84, 37.21]#Mw7.6
meca_file = '../data/meca_turkey_6.dat'
outgmt_name = 'Turkey_2023_%s_'%para['delta_time']
out_gmt_locate = '../out_file/b_spatial_temproal_variation/'
tmp_time = 0
gmt_draw.gmt_drawing(para['roi'], faults_file_name, meca_file, outgmt_name, tmp_time, out_gmt_locate)
