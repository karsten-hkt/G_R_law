import utm
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from sympy import Plane, Point3D
# 更改Matplotlib的字体设置
plt.rcParams['font.family'] = 'Hiragino Sans GB'
# 在此之后绘制图形

#这个求解的太慢了，淘汰，请使用test_point2plane_improve

def strdip2plane(point, strike, dip):
    if dip<0 or dip>90:
        raise('dip need to [0,90]')
    if strike<0 or strike>=360:
        raise('dip need to [0,360)')
    dip = math.radians(dip)
    strike = math.radians(strike)
    n = (math.cos(math.pi/2-dip)*math.sin(math.pi/2+strike),
         math.cos(math.pi/2-dip)*math.cos(math.pi/2+strike),
         math.sin(math.pi/2-dip))
    P1 = Plane(Point3D(point[0],point[1],point[2]), normal_vector = n)
    return P1

def project_points_to_plane_3d_with_height(points_latlon_alt, fault_point, fault_strike, fault_dip):

    #全部输入都需要使用m计算，所以对于经纬度，要修改成utm投影，对于depth要修改成m

    projected_points_latlon_alt = []

    # 将平面内的已知点坐标转换为UTM坐标
    fault_point[0], fault_point[1], _, _ = utm.from_latlon(fault_point[1], fault_point[0])

    fault_plane = strdip2plane(fault_point, fault_strike, fault_dip)

    for point_latlon_alt in points_latlon_alt:
        # 确定经纬度点所属的UTM投影区域 lat lon
        utm_x, utm_y, utm_zone_number, utm_zone_letter = utm.from_latlon(point_latlon_alt[1], point_latlon_alt[0])

        # 点转换成Point3D格式，深度需要转换成m
        points_latlon_alt_3D = Point3D(utm_x, utm_y, -point_latlon_alt[2]*1000)

        # 计算投影点P'
        projectionPoint = fault_plane.projection(points_latlon_alt_3D)

        # 将UTM坐标点转换为经纬度坐标点
        lat,lon = utm.to_latlon(float(projectionPoint[0]),float(projectionPoint[1]), utm_zone_number, utm_zone_letter)

        # 添加投影后的高度信息
        projected_points_latlon_alt.append([lon,lat,float(projectionPoint[2])/1000])
        print([lon,lat,projectionPoint[2]/1000])
    return projected_points_latlon_alt


fault_name = '../data/turkey_7_8_fault.txt'
data_name = '../data/data_processing.csv'
fault = pd.read_csv(fault_name)
data = pd.read_csv(data_name)
data = np.array(data.iloc[:,[2,3,4]])
data = data[:10,:]
# 调用函数进行投影
fault_point = [fault.iloc[0, 1], fault.iloc[0, 2], 0]
projected_points_latlon_alt = project_points_to_plane_3d_with_height(data, fault_point,
                                                                     fault.iloc[0, 5], 80)

print(projected_points_latlon_alt)
