import numpy as np
import cmath
import math
## utm转换带时候需要的参数
from shapely.geometry import LineString
import utm

K0 = 0.9996
ZONE_LETTERS = [
    (84, None), (72, 'X'), (64, 'W'), (56, 'V'), (48, 'U'), (40, 'T'),
    (32, 'S'), (24, 'R'), (16, 'Q'), (8, 'P'), (0, 'N'), (-8, 'M'), (-16, 'L'),
    (32, 'S'), (24, 'R'), (16, 'Q'), (8, 'P'), (0, 'N'), (-8, 'M'), (-16, 'L'),
    (-24, 'K'), (-32, 'J'), (-40, 'H'), (-48, 'G'), (-56, 'F'), (-64, 'E'),
    (-72, 'D'), (-80, 'C')
]
R = 6378137#地球半径
E = 0.00669438
E2 = E * E
E3 = E2 * E
E_P2 = E / (1.0 - E)
M1 = (1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256)
M2 = (3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024)
M3 = (15 * E2 / 256 + 45 * E3 / 1024)
M4 = (35 * E3 / 3072)
__all__ = ['to_latlon', 'from_latlon']

K0 = 0.9996

E = 0.00669438
E2 = E * E
E3 = E2 * E
E_P2 = E / (1.0 - E)

SQRT_E = math.sqrt(1 - E)
_E = (1 - SQRT_E) / (1 + SQRT_E)
_E2 = _E * _E
_E3 = _E2 * _E
_E4 = _E3 * _E
_E5 = _E4 * _E

M1 = (1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256)
M2 = (3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024)
M3 = (15 * E2 / 256 + 45 * E3 / 1024)
M4 = (35 * E3 / 3072)

P2 = (3. / 2 * _E - 27. / 32 * _E3 + 269. / 512 * _E5)
P3 = (21. / 16 * _E2 - 55. / 32 * _E4)
P4 = (151. / 96 * _E3 - 417. / 128 * _E5)
P5 = (1097. / 512 * _E4)

R = 6378137

ZONE_LETTERS = [
    (84, None), (72, 'X'), (64, 'W'), (56, 'V'), (48, 'U'), (40, 'T'),
    (32, 'S'), (24, 'R'), (16, 'Q'), (8, 'P'), (0, 'N'), (-8, 'M'), (-16, 'L'),
    (-24, 'K'), (-32, 'J'), (-40, 'H'), (-48, 'G'), (-56, 'F'), (-64, 'E'),
    (-72, 'D'), (-80, 'C')
]

def is_in_poly(p, poly):
    """
    :param p: [x, y]
    :param poly: [[], [], [], [], ...]
    :return:
    """
    px, py = p
    is_in = False
    for i, corner in enumerate(poly):
        next_i = i + 1 if i + 1 < len(poly) else 0
        x1, y1 = corner
        x2, y2 = poly[next_i]
        if (x1 == px and y1 == py) or (x2 == px and y2 == py):  # if point is on vertex
            is_in = True
            break
        if min(y1, y2) < py <= max(y1, y2):  # find horizontal edges of polygon
            x = x1 + (py - y1) * (x2 - x1) / (y2 - y1)
            if x == px:  # if point is on edge
                is_in = True
                break
            elif x > px:  # if point is on left-side of line
                is_in = not is_in
    return is_in

## 旋转坐标轴从而得到余震点距离震中的距离
def simp_rotaxll(x,y,strike,cx0,cy0):
    x0, y0, zn, zl = from_latlon(cy0,cx0)
    x0, y0 = x0 / 1000., y0 / 1000.
    x, y = lltoutm(y,x,force_zone_number=zn,force_zone_letter=zl)
    x, y = x / 1000., y / 1000.
    outx, outy = simp_rotax(x,y,strike,x0,y0)
    #print(d)
    return outx, outy

## umt投影转化为km
## Updated by Wanpeng Feng, to allow an input zone letter
## @Ottawa, 2016-06-20
def from_latlon(latitude, longitude, force_zone_number=None,force_zone_letter=None):
    if isinstance(latitude,(np.ndarray)):
        if sum(latitude.ravel()>84.0) > 0:
            raise OutOfRangeError('latitude out of range \
                                 (must be between 80 deg S and 84 deg N)')
        longitude[longitude>180] = longitude[longitude>180] - 3
    else:
        if not -80.0 <= latitude <= 84.0:
            raise OutOfRangeError('latitude out of range \
                                (must be between 80 deg S and 84 deg N)')
        if not -180.0 <= longitude <= 180.0:
            raise OutOfRangeError('northing out of range \
                                (must be between 180 deg W and 180 deg E)')
    # replacing math by np to allow matrix operation
    lat_rad = np.radians(latitude)
    lat_sin = np.sin(lat_rad)
    lat_cos = np.cos(lat_rad)
    lat_tan = lat_sin / lat_cos
    lat_tan2 = lat_tan * lat_tan
    lat_tan4 = lat_tan2 * lat_tan2
    if force_zone_number is None:
        zone_number = latlon_to_zone_number(latitude, longitude)
    else:
        zone_number = force_zone_number
    if force_zone_letter is None:
        zone_letter = latitude_to_zone_letter(latitude)
    else:
        zone_letter = force_zone_letter
    lon_rad = np.radians(longitude)
    central_lon = zone_number_to_central_longitude(zone_number)
    central_lon_rad = np.radians(central_lon)
    n = R / np.sqrt(1 - E * lat_sin**2)
    c = E_P2 * lat_cos**2
    a = lat_cos * (lon_rad - central_lon_rad)
    a2 = a * a
    a3 = a2 * a
    a4 = a3 * a
    a5 = a4 * a
    a6 = a5 * a
    m = R * (M1 * lat_rad -
             M2 * np.sin(2 * lat_rad) +
             M3 * np.sin(4 * lat_rad) -
             M4 * np.sin(6 * lat_rad))
    easting = K0 * n * (a +
                        a3 / 6 * (1 - lat_tan2 + c) +
                        a5 / 120 * (5 - 18 * lat_tan2 + lat_tan4 + 72 * c - 58 * E_P2)) + 500000
    northing = K0 * (m + n * lat_tan * (a2 / 2 +
                                        a4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c**2) +
                                        a6 / 720 * (61 - 58 * lat_tan2 + lat_tan4 + 600 * c - 330 * E_P2)))
    if isinstance(latitude, np.ndarray):
        northing[latitude<0] = northing[latitude<0] + 10000000
    else:
        if latitude < 0:
            northing += 10000000
    return easting, northing, zone_number, zone_letter

## 投影坐标转换
# Added by wanpeng feng, @Ottawa, 2017-01-03
# Projection convertion on an array based method
def lltoutm(latm,lonm,force_zone_number=None,force_zone_letter=None):
    #
    #
    # latm = np.array(latm)
    # lonm = np.array(lonm)
    # utmx = np.copy(lonm * 0)
    # utmy = np.copy(utmx)
    utmx = np.copy(lonm * 0)
    utmy = np.copy(utmx)
    #
    if utmx.size != 1:
        for index in range(utmx.shape[0]):
            #
            e,n,z,l = from_latlon(latm[index],lonm[index], force_zone_number=force_zone_number,\
                        force_zone_letter=force_zone_letter)
            utmx[index] = e
            utmy[index] = n
        #
        return utmx,utmy
    else:
        e, n, z, l = from_latlon(latm, lonm, force_zone_number=force_zone_number, \
                                 force_zone_letter=force_zone_letter)
        utmx = e
        utmy = n
        return utmx,utmy

## 旋转坐标系
def simp_rotax(x,y,strike,cx0,cy0):
    #
    import matplotlib.pyplot as plt
    import math
    x, y    = np.array(x), np.array(y)
    PI      = np.pi
    strkr   = (90.-strike)*PI/180.;
    istrkr  = complex(0,strkr)
    xystart = complex(cx0,cy0);
    cx      = np.copy(x)
    cy      = np.copy(y)
    for ni in range(cx.shape[0]):
        outdata = (complex(x[ni],y[ni])-xystart) * cmath.exp(-1. * istrkr)
        cx[ni] = outdata.real
        cy[ni] = outdata.imag
    return cx,cy

## utm带号
def latlon_to_zone_number(latitude, longitude):
    if 56 <= latitude <= 64 and 3 <= longitude <= 12:
        return 32
    if 72 <= latitude <= 84 and longitude >= 0:
        if longitude <= 9:
            return 31
        elif longitude <= 21:
            return 33
        elif longitude <= 33:
            return 35
        elif longitude <= 42:
            return 37
    return int((longitude + 180) / 6) + 1
def latitude_to_zone_letter(latitude):
    #
    for lat_min, zone_letter in ZONE_LETTERS:
        if latitude >= lat_min:
            return zone_letter

    return None
def to_latlon(easting, northing, zone_number, zone_letter=None, northern=None):
    #
    if not zone_letter and northern is None:
        raise ValueError('either zone_letter or northern needs to be set')

    elif zone_letter and northern is not None:
        raise ValueError('set either zone_letter or northern, but not both')
    if not isinstance(easting,(np.ndarray)):
      #
      #if not 100000 <= easting < 1000000:
      #    raise OutOfRangeError('easting out of range (must be between 100.000 m and 999.999 m)')
      #if not 0 <= northing <= 10000000:
      #    raise OutOfRangeError('northing out of range (must be between 0 m and 10.000.000 m)')
      if not 1 <= zone_number <= 60:
          raise OutOfRangeError('zone number out of range (must be between 1 and 60)')
    #
    if zone_letter:
        zone_letter = zone_letter.upper()
        if not 'C' <= zone_letter <= 'X' or zone_letter in ['I', 'O']:
            raise OutOfRangeError('zone letter out of range (must be between C and X)')
        northern = (zone_letter >= 'N')
    #
    x = easting - 500000
    y = northing

    if not northern:
        y -= 10000000

    m = y / K0
    mu = m / (R * M1)

    p_rad = (mu +\
             P2 * np.sin(2 * mu) +\
             P3 * np.sin(4 * mu) +\
             P4 * np.sin(6 * mu) +\
             P5 * np.sin(8 * mu))

    p_sin = np.sin(p_rad)
    p_sin2 = p_sin * p_sin

    p_cos = np.cos(p_rad)

    p_tan = p_sin / p_cos
    p_tan2 = p_tan * p_tan
    p_tan4 = p_tan2 * p_tan2

    ep_sin = 1 - E * p_sin2
    ep_sin_sqrt = np.sqrt(1 - E * p_sin2)

    n = R / ep_sin_sqrt
    r = (1 - E) / ep_sin

    c = _E * p_cos**2
    c2 = c * c

    d = x / (n * K0)
    d2 = d * d
    d3 = d2 * d
    d4 = d3 * d
    d5 = d4 * d
    d6 = d5 * d

    latitude = (p_rad - (p_tan / r) *
                (d2 / 2 -
                 d4 / 24 * (5 + 3 * p_tan2 + 10 * c - 4 * c2 - 9 * E_P2)) +
                 d6 / 720 * (61 + 90 * p_tan2 + 298 * c + 45 * p_tan4 - 252 * E_P2 - 3 * c2))

    longitude = (d -
                 d3 / 6 * (1 + 2 * p_tan2 + c) +
                 d5 / 120 * (5 - 2 * c + 28 * p_tan2 - 3 * c2 + 8 * E_P2 + 24 * p_tan4)) / p_cos

    return (np.degrees(latitude),
            np.degrees(longitude) + zone_number_to_central_longitude(zone_number))

def utmtoll(ex,ny,force_zone_number,force_zone_letter=None):
    #
    olat = np.copy(ex)
    olon = np.copy(ny)
    olat = olat * 0.
    olon = olon * 0.
    if olat.size != 1:
        for index in range(olat.shape[0]):
            olat[index],olon[index] = to_latlon(ex[index], ny[index],
                force_zone_number, zone_letter=force_zone_letter)
    else:
        olat, olon = to_latlon(ex, ny,
                force_zone_number, zone_letter=force_zone_letter)
        #
    return olat,olon

def fault_parameter_genereate(fault_lon,fault_lat,fault_length,fault_strike):
    zone = None
    nl = None
    fault_para = []
    # 生成x，y
    ux, uy, a, b = from_latlon(fault_lat, fault_lon, force_zone_number=zone,
                                                   force_zone_letter=nl)
    fault_para.append([fault_lon, fault_lat])
    # ux for easting, uy for northing
    # 加上断层的长度以及走向
    ux = ux + math.sin(fault_strike * math.pi / 180) * 1000 * fault_length
    uy = uy + math.cos(fault_strike * math.pi / 180) * 1000 * fault_length
    # start to build fault one by one
    lat, lon = to_latlon(ux, uy, a, b, northern=None)
    fault_para.append([lon, lat])
    return fault_para

def buffer_generate(fault_para, length):
    nl = latitude_to_zone_letter(fault_para[0,1])
    zone = latlon_to_zone_number(fault_para[0,1], fault_para[0,0])
    ux,uy = lltoutm(fault_para[:, 1], fault_para[:, 0], force_zone_number=zone, force_zone_letter=nl)
    fault_line = [[ux[0],uy[0]],[ux[1],uy[1]]]
    line = LineString(fault_line)
    buffer = line.buffer(length)
    x2, y2 = buffer.boundary.xy
    lat,lon = utmtoll(x2,y2,force_zone_number=zone,force_zone_letter=nl)
    buffer_build = []
    for i in range(len(lat)):
        buffer_build.append([lon[i],lat[i]])
    return buffer_build

def zone_number_to_central_longitude(zone_number):
    return (zone_number - 1) * 6 - 180 + 3

def length2ll(length, cx0, cy0, strike):
    nl = latitude_to_zone_letter(cy0)
    zone = latlon_to_zone_number(cy0, cx0)
    # 由于lltoutm内部判断了ndarray的存在，因此在此处需要吧cy0，cx0转换成list的形式
    ux0, uy0 = lltoutm(cy0, cx0, force_zone_number=zone, force_zone_letter=nl)
    ux = np.zeros(len(length))
    uy = np.zeros(len(length))
    for i in range(len(length)):
        ux[i] = ux0 + math.sin(strike * math.pi / 180) * 1000 * length[i]
        uy[i] = uy0 + math.cos(strike * math.pi / 180) * 1000 * length[i]
    lat, lon = utmtoll(ux, uy, force_zone_number=zone, force_zone_letter=nl)
    out_utm = [lon,lat]
    return out_utm

def center24points2ll(now_poly, c0, strike_vector, dip_vector):
    '''

    :param now_poly: 输入的矩阵
    :param c0: 中心点
    :param strike_vector: 走向方向的向量
    :param dip_vector: 倾向方向的向量
    :return:
    四个方形点所在的ll坐标
    '''
    nl = latitude_to_zone_letter(c0[1])
    zone = latlon_to_zone_number(c0[1], c0[0])
    # 由于lltoutm内部判断了ndarray的存在，因此在此处需要吧cy0，cx0转换成list的形式
    ux0, uy0 = lltoutm(c0[1], c0[0], force_zone_number=zone, force_zone_letter=nl)
    u0 = np.array([ux0,uy0,c0[2]*1000])#m为单位的utm坐标
    points = []
    for i in range(len(now_poly)):
        point = u0 + now_poly[i,0] * 1000 * strike_vector + now_poly[i,1]* 1000 * dip_vector
        points.append(point)
    points = np.array(points)
    lat, lon = utmtoll(points[:,0], points[:,1], force_zone_number=zone, force_zone_letter=nl)
    points[:,[0,1]] = np.array([lon,lat]).T
    points[:,2] = points[:,2] / 1000
    return points

def simp_rotaxll2utm(x,y,strike,cx0,cy0):
    outx, outy = simp_rotaxll(x,y,strike,cx0,cy0)
    out_utm = length2ll(outx,cx0,cy0,strike)
    return out_utm

def haversine(lat1, lon1, lat2, lon2):
    '''
    两个经纬度之间的距离，km
    :param lat1:
    :param lon1:
    :param lat2:
    :param lon2:
    :return:
    '''
    # 将经纬度转换为弧度
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    # 地球半径（单位：千米）
    radius = 6371.0

    # Haversine公式
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # 计算距离
    distance = radius * c

    return distance

def strdip2plane(strike, dip):
    if dip<0 or dip>180:
        raise('dip need to [0,180]')
    if strike<0 or strike>=360:
        raise('dip need to [0,360)')
    dip = math.radians(dip)
    strike = math.radians(strike)
    n = np.array([math.cos(math.pi/2-dip)*math.sin(math.pi/2+strike),
         math.cos(math.pi/2-dip)*math.cos(math.pi/2+strike),
         math.sin(math.pi/2-dip)])
    return n

def dot_pro2plane(points_latlon_alt, fault_point, fault_strike, fault_dip):
    '''

    :param points_latlon_alt: 需要投影点点
    :param fault_point: 断层上点点（深度需要为负数）
    :param fault_strike: 断层走向
    :param fault_dip: 断层倾角
    :return:
    projected_points_latlon_alt: 投影后在平面上所在点
    x_dis: 距离走向的距离
    y_dis: 距离垂向的距离
    '''
    projected_points_latlon_alt = []
    strike_xs = []
    strike_ys = []

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
        strike_xs.append(x_dis)
        # 想要得到倾向上的值时候,倾向上的单位向量是
        dip_line = np.cross(strike_line, normal_vector)
        y_dis = np.dot(dip_line, projected_vector)
        strike_ys.append(y_dis)
        # 将UTM坐标点转换为经纬度坐标点
        lat, lon = utm.to_latlon(projected_point[0], projected_point[1], utm_zone_number, utm_zone_letter)
        projected_points_latlon_alt.append([lon,lat,projected_point[2]/1000])

    return projected_points_latlon_alt, strike_xs, strike_ys

if __name__ == '__main__':
    print('This is for location trans')
