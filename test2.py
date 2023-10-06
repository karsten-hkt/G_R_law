import utm
import numpy as np
import matplotlib.pyplot as plt

# 更改Matplotlib的字体设置
plt.rcParams['font.family'] = 'Hiragino Sans GB'
# 在此之后绘制图形

def project_points_to_plane_3d(points_latlon_alt, plane_point, plane_direction, plane_trend):
    projected_points_latlon_alt = []

    # 计算平面的法向量（使用极坐标转换）
    plane_direction_radians = np.radians(plane_direction)
    plane_trend_radians = np.radians(plane_trend)

    normal_vector = np.array([
        np.sin(plane_trend_radians) * np.cos(plane_direction_radians),
        np.sin(plane_trend_radians) * np.sin(plane_direction_radians),
        np.cos(plane_trend_radians)
    ])

    for point_latlon_alt in points_latlon_alt:
        # 确定经纬度点所属的UTM投影区域
        utm_x, utm_y, utm_zone_number, utm_zone_letter = utm.from_latlon(point_latlon_alt[1], point_latlon_alt[0])

        # 将平面内的已知点坐标转换为UTM坐标
        plane_utm_x, plane_utm_y, _, _ = utm.from_latlon(plane_point[1], plane_point[0])

        # 计算点到平面上一点的向量V，考虑点的高度
        V = np.array([utm_x - plane_utm_x, utm_y - plane_utm_y, point_latlon_alt[2]])

        # 计算点到平面上一点的距离d（点到平面的距离）
        d = np.dot(V, normal_vector) / np.linalg.norm(normal_vector)

        # 计算投影点P'
        projected_point = np.array([utm_x, utm_y, point_latlon_alt[2]]) - d * normal_vector

        # 将UTM坐标点转换为经纬度坐标点
        lat, lon = utm.to_latlon(projected_point[0], projected_point[1], utm_zone_number, utm_zone_letter)
        projected_points_latlon_alt.append([lon,lat,point_latlon_alt[2]])
    return projected_points_latlon_alt

# 示例使用
# 定义平面内的已知点、走向角度和倾角角度
plane_point = [-74.006, 40.7128]  # 纽约市的中心点（经纬度）
plane_direction = 45  # 平面的走向角度（以度为单位）
plane_trend = 90  # 平面的倾角角度（以度为单位）

# 定义多个经纬度坐标点（包括高度信息）
points_latlon_alt = [
    [-74.007, 40.7129, -10],  # 纽约市的另一个点，高度为-10
    [-73.005, 41.7127, -5],   # 纽约市的另一个点，高度为-5
    [-76.007, 42.7129, -8],  # 纽约市的另一个点，高度为-10
    [-75.005, 43.7127, -9],
]

# 调用函数进行投影
projected_points_latlon_alt = project_points_to_plane_3d(points_latlon_alt, plane_point, plane_direction, plane_trend)



latitudes = [point[1] for point in points_latlon_alt]
longitudes = [point[0] for point in points_latlon_alt]
altitudes = [point[2] for point in points_latlon_alt]
projected_latitudes = [point[0] for point in projected_points_latlon_alt]
projected_longitudes = [point[1] for point in projected_points_latlon_alt]
projected_alt = [point[2] for point in projected_points_latlon_alt]
# 打印投影后的结果，包括高度信息
for i, projected_point in enumerate(projected_points_latlon_alt):
    print(f"点{i + 1}：经度={projected_point[0]:.6f}，纬度={projected_point[1]:.6f}，高度={projected_point[2]:.2f}")

# 绘制三维平面和投影点（与高度信息）
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 经纬度点
ax.scatter(longitudes, latitudes, altitudes, c='r', marker='^', label='经纬度点')

# 投影点
ax.scatter(projected_longitudes, projected_latitudes, projected_alt, c='b',
           marker='x', label='投影点')


# 设置坐标轴标签
ax.set_xlabel('经度')
ax.set_ylabel('纬度')
ax.set_zlabel('高度')

# 显示图形
plt.title('三维平面投影（包括高度信息）')
plt.show()

