import pandas as pd
import math
'''
调用函数去投影太慢了，升级一下，直接用向量投影
'''

import utm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
# 更改Matplotlib的字体设置
plt.rcParams['font.family'] = 'Hiragino Sans GB'
# 在此之后绘制图形

def strdip2plane(strike, dip):
    if dip<0 or dip>90:
        raise('dip need to [0,90]')
    if strike<0 or strike>=360:
        raise('strike need to [0,360)')
    dip = math.radians(dip)
    strike = math.radians(strike)
    n = np.array([math.cos(math.pi/2-dip)*math.sin(math.pi/2+strike),
         math.cos(math.pi/2-dip)*math.cos(math.pi/2+strike),
         math.sin(math.pi/2-dip)])
    return n

def project_points_to_plane_3d(points_latlon_alt, fault_point, fault_strike, fault_dip):
    projected_points_latlon_alt = []
    utm_points = []
    projected_points = []
    # 计算平面的法向量（使用极坐标转换）
    normal_vector = strdip2plane(fault_strike,fault_dip)

    for point_latlon_alt in points_latlon_alt:
        # 确定经纬度点所属的UTM投影区域
        utm_x, utm_y, utm_zone_number, utm_zone_letter = utm.from_latlon(point_latlon_alt[1], point_latlon_alt[0])

        # 将平面内的已知点坐标转换为UTM坐标
        fault_point_utm_x, fault_point_utm_y, _, _ = utm.from_latlon(fault_point[1], fault_point[0])

        # 计算点到平面上一点的向量V，考虑点的高度
        V = np.array([utm_x - fault_point_utm_x , utm_y - fault_point_utm_y,
                      point_latlon_alt[2]*1000 - fault_point[2]*1000])

        # 计算投影向量
        projected_vector = V - np.dot(V, normal_vector) / np.dot(normal_vector,normal_vector) * normal_vector

        # 投影点
        projected_point = [projected_vector[0] + fault_point_utm_x, projected_vector[1] + fault_point_utm_y,
                           projected_vector[2] + fault_point[2]*1000]

        # 想要得到走向上的值时候,走向上的单位向量是
        strike_line = np.array([math.sin(math.radians(fault_strike)),math.cos(math.radians(fault_strike)),0])
        x_dis = np.dot(strike_line, projected_vector)
        # 想要得到倾向上的值时候,倾向上的单位向量是
        dip_line = np.cross(strike_line, normal_vector)
        y_dis = np.dot(dip_line, projected_vector)
        # 将UTM坐标点转换为经纬度坐标点
        lat, lon = utm.to_latlon(projected_point[0], projected_point[1], utm_zone_number, utm_zone_letter)

        projected_points_latlon_alt.append([lon,lat,projected_point[2]/1000])
        utm_points.append([utm_x, utm_y, point_latlon_alt[2]*1000])
        projected_points.append(projected_point)

    return projected_points_latlon_alt, projected_points, utm_points, x_dis, y_dis

fault_name = '../data/turkey_7_8_fault.txt'
data_name = '../data/data_processing.csv'
fault = pd.read_csv(fault_name)
data = pd.read_csv(data_name)
#需要将depth的值转为负数，要不生成的坐标系并不是传统的坐标系统
data.iloc[:,4] = -data.iloc[:,4]
data = np.array(data.iloc[:,[2,3,4]])
data = data[:100,:]
# 调用函数进行投影
fault_point = [fault.iloc[0, 1], fault.iloc[0, 2], 0]
#print(fault.iloc[0, 5], fault.iloc[0,6])
projected_points_latlon_alt,project_point, utm_points, x, y =\
    project_points_to_plane_3d(data, fault_point, fault.iloc[0, 5], 70)
print(data)
print(projected_points_latlon_alt)
#print(projected_points_latlon_alt)
utm_points = np.array(utm_points)
project_point = np.array(project_point)
fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection = '3d')
ax.scatter3D(utm_points[:,0],utm_points[:,1],utm_points[:,2])
ax.scatter3D(project_point[:,0],project_point[:,1],project_point[:,2])
plt.show()
