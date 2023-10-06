import numpy as np
import math
import matplotlib.pyplot as plt

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

def project_points_to_plane(points, fault_point, fault_strike, fault_dip):
    # 计算平面的单位法向量（使用极坐标转换）
    normal_vector = strdip2plane(fault_strike, fault_dip)
    print('法向量：',normal_vector)
    for point in points:
        # 计算点到平面上一点的向量V，考虑点的高度
        V = np.array([point[0] - fault_point[0], point[1] - fault_point[1], point[2] - fault_point[2]])

        #计算投影点向量
        projected_vector = V - np.dot(V, normal_vector) / np.dot(normal_vector,normal_vector) * normal_vector

        #投影点
        projected_point = fault_point + projected_vector
        print('投影向量',projected_vector)
        print('投影点',projected_point)
        # 想要得到走向上的值时候,走向上的单位向量是
        strike_line = np.array([math.sin(math.radians(fault_strike)), math.cos(math.radians(fault_strike)), 0])
        x_dis = np.dot(strike_line, projected_vector)
        # 想要得到倾向上的值时候,倾向上的单位向量是
        dip_line = np.cross(strike_line, normal_vector)
        print('strike_Line',strike_line)
        print('dip_line',dip_line)
        y_dis = np.dot(dip_line, projected_point)
        print('距离走向：',x_dis, '距离垂向:', y_dis)
        # 将UTM坐标点转换为经纬度坐标点

    return projected_point, x_dis, y_dis

data = np.array([[1,1,1],[2,0,0]])
fault_point = np.array([1,0,0])
strike = 90
dip = 0
projected_points_latlon_alt = project_points_to_plane(data, fault_point, strike, dip)
print(projected_points_latlon_alt)
