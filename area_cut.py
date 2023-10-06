import numpy as np
import inpoly
import sys
sys.path.append('../src/PYTHONPATH')
import pSAR
import math
from shapely.geometry import LineString
import location_trans as LocT

def ext2polygon(ext):
    outpoly = np.zeros([5,2])
    outpoly[0,:] = [ext[0],ext[2]]
    outpoly[1,:] = [ext[0],ext[3]]
    outpoly[2,:] = [ext[1],ext[3]]
    outpoly[3,:] = [ext[1],ext[2]]
    outpoly[4,:] = [ext[0],ext[2]]
    return outpoly

def in_ploy(data, lat_line, lon_line, roi):
    polys = ext2polygon(roi)
    ## 使用inploy函数
    IN, ON = inpoly.inpoly2(data.iloc[:, [lon_line, lat_line]], polys)
    ## 选出在内部的点
    data_in_ploy = data[IN==True]
    return data_in_ploy

## 选取在一定半径内部的余震
## 需要传入此时的经度，纬度，经度所在行，纬度所在行，圆半径，余震序列
def circle_file(lon_min,lon_max,lat_min,lat_max,lon_line,lat_line,distance,aftershocks_file):
    x_mid = (lon_min+lon_max)/2
    y_mid = (lat_min+lat_max)/2
    use_file = ((x_mid-aftershocks_file[:,lon_line])**2+(y_mid-aftershocks_file[:,lat_line])**2)**0.5
    af_file = aftershocks_file[use_file[:]<=distance,:]
    return af_file

## 将矩形转换成gmt可以读取的格式
## 需要传入矩形，b值，输出名字，此时时间间隔
def ref2gmt(points_polys,bvalues,outgmt,tmp_time,gmt_locate):
    file_name = gmt_locate + outgmt + ('%s'%tmp_time) + '.gmt'
    fid = open(file_name,'w')
    for i in range(bvalues.shape[0]):
        if math.isnan(bvalues[i]) == False:
            fid.write('>-Z%4.3f\n' % bvalues[i])
            for j in range(5):
                fid.write('%f %f\n' % (points_polys[i][j,0],points_polys[i][j,1]))
    fid.close()
    return

def ref2gmt_3d(points_polys, bvalues, strike, cx0, cy0, outgmt,tmp_time,gmt_locate):
    #add new in an unexisting file
    '''
    数据格式
    # 数据头段
    >  -Z1.800000e-01
     # 经度 纬度 深度
     99.3785  34.5324  -0.0010
     99.3574  34.5279  -0.0010
     99.3569  34.5297  -1.9900
     99.3779  34.5342  -1.9900
    >  -Z5.000000e-02
     99.3574  34.5279  -0.0010
     99.3363  34.5233  -0.0010
     99.3358  34.5252  -1.9900
     99.3569  34.5297  -1.9900
    ...
    '''
    file_name = gmt_locate + outgmt + ('%s'%tmp_time) + '.gmt'
    fid = open(file_name,'w')
    for i in range(bvalues.shape[0]):
        if math.isnan(bvalues[i]) == False:
            fid.write('>-Z%4.3f\n' % bvalues[i])
            now_poly = np.array(points_polys[i])
            out_af_utm = LocT.length2ll(now_poly[:,0], cx0, cy0, strike)
            for j in range(5):
                fid.write('%f %f %f\n' % (out_af_utm[0][j],out_af_utm[1][j],-points_polys[i][j,1]))
    fid.close()
    return

def ref2gmt_3d_dip(points_polys, bvalues, strike, dip, cx0, cy0, outgmt,tmp_time,gmt_locate):
    #add new in an unexisting file
    '''
    数据格式
    # 数据头段
    >  -Z1.800000e-01
     # 经度 纬度 深度
     99.3785  34.5324  -0.0010
     99.3574  34.5279  -0.0010
     99.3569  34.5297  -1.9900
     99.3779  34.5342  -1.9900
    >  -Z5.000000e-02
     99.3574  34.5279  -0.0010
     99.3363  34.5233  -0.0010
     99.3358  34.5252  -1.9900
     99.3569  34.5297  -1.9900
    ...
    '''
    file_name = gmt_locate + outgmt + ('%s'%tmp_time) + '.gmt'
    fid = open(file_name,'w')
    strike_vector = np.array([math.sin(math.radians(strike)), math.cos(math.radians(strike)), 0])
    normal_vector = LocT.strdip2plane(strike,dip)
    dip_vector = np.cross(strike_vector, normal_vector)
    c0 = [cx0,cy0, 0]
    for i in range(bvalues.shape[0]):
        if math.isnan(bvalues[i]) == False:
            fid.write('>-Z%4.3f\n' % bvalues[i])
            now_poly = np.array(points_polys[i])
            out_af_utm = LocT.center24points2ll(now_poly, c0, strike_vector, dip_vector)
            for j in range(5):
                fid.write('%f %f %f\n' % (out_af_utm[j][0],out_af_utm[j][1],out_af_utm[j][2]))
    fid.close()
    return


def ref2gmt_3d_exist(points_polys, bvalues, strike, cx0, cy0, outgmt,tmp_time,gmt_locate):
    #add new in the existing file
    '''
    数据格式
    # 数据头段
    >  -Z1.800000e-01
     # 经度 纬度 深度
     99.3785  34.5324  -0.0010
     99.3574  34.5279  -0.0010
     99.3569  34.5297  -1.9900
     99.3779  34.5342  -1.9900
    >  -Z5.000000e-02
     99.3574  34.5279  -0.0010
     99.3363  34.5233  -0.0010
     99.3358  34.5252  -1.9900
     99.3569  34.5297  -1.9900
    ...
    '''
    file_name = gmt_locate + outgmt + ('%s'%tmp_time) + '.gmt'
    fid = open(file_name,'a')
    for i in range(bvalues.shape[0]):
        if math.isnan(bvalues[i]) == False:
            fid.write('>-Z%4.3f\n' % bvalues[i])
            now_poly = np.array(points_polys[i])
            out_af_utm = LocT.length2ll(now_poly[:,0], cx0, cy0, strike)
            for j in range(5):
                fid.write('%f %f %f\n' % (out_af_utm[0][j],out_af_utm[1][j],-points_polys[i][j,1]))
    fid.close()
    return

def ref2gmt_3d_dip_exist(points_polys, bvalues, strike, dip, cx0, cy0, outgmt,tmp_time,gmt_locate):
    #add new in an unexisting file
    '''
    数据格式
    # 数据头段
    >  -Z1.800000e-01
     # 经度 纬度 深度
     99.3785  34.5324  -0.0010
     99.3574  34.5279  -0.0010
     99.3569  34.5297  -1.9900
     99.3779  34.5342  -1.9900
    >  -Z5.000000e-02
     99.3574  34.5279  -0.0010
     99.3363  34.5233  -0.0010
     99.3358  34.5252  -1.9900
     99.3569  34.5297  -1.9900
    ...
    '''
    file_name = gmt_locate + outgmt + ('%s'%tmp_time) + '.gmt'
    fid = open(file_name,'a')
    strike_vector = np.array([math.sin(math.radians(strike)), math.cos(math.radians(strike)), 0])
    normal_vector = LocT.strdip2plane(strike,dip)
    dip_vector = np.cross(strike_vector, normal_vector)
    c0 = [cx0,cy0,0]
    for i in range(bvalues.shape[0]):
        if math.isnan(bvalues[i]) == False:
            fid.write('>-Z%4.3f\n' % bvalues[i])
            now_poly = np.array(points_polys[i])
            out_af_utm = LocT.center24points2ll(now_poly, c0, strike_vector, dip_vector)
            for j in range(5):
                fid.write('%f %f %f\n' % (out_af_utm[j][0],out_af_utm[j][1],out_af_utm[j][2]))
    fid.close()
    return


def record_polys(lon_min_now,lon_max_now,lat_min_now,lat_max_now):
    out_data = np.zeros((5,2))
    out_data[0,:] = [lon_min_now,lat_min_now]
    out_data[1,:] = [lon_min_now,lat_max_now]
    out_data[2,:] = [lon_max_now,lat_max_now]
    out_data[3,:] = [lon_max_now,lat_min_now]
    out_data[4,:] = [lon_min_now,lat_min_now]
    return out_data

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
#

## 绘图的区块
def ext2refpts(ext,spacing):
    nx = int((ext[1]-ext[0])/spacing)+2
    ny = int((ext[3]-ext[2])/spacing)+2
    outdata = []
    polys = []
    for i in range(nx):
        for j in range(ny):
            outdata.append([ext[0]+i*spacing,ext[2]+j*spacing])
            cext = [ext[0]+i*spacing - spacing*0.5,ext[0]+i*spacing + spacing*0.5,ext[2]+j*spacing-0.5*spacing,ext[2]+j*spacing+0.5*spacing]
            polys.append(ext2polygon(cext))
    return np.array(outdata),polys


## 一定半径范围内的区块
def refpts2bufferpoly(refpts,spacing,buf=1):
    outdata = []
    for i in range(refpts.shape[0]):
        #
        cdata = np.zeros([5,2])
        tmp_ext = [refpts[i,0]-spacing*buf,refpts[i,0]+spacing*buf,refpts[i,1]-spacing*buf,refpts[i,1]+spacing*buf]
        polygons = ext2polygon(tmp_ext)
        outdata.append(polygons)
    return outdata

def fault_parameter_genereate(fault_lon,fault_lat,fault_length,fault_strike):
    zone = None
    nl = None
    fault_para = []
    # 生成x，y
    ux, uy, a, b = pSAR.utm_conversion.from_latlon(fault_lat, fault_lon, force_zone_number=zone,
                                                   force_zone_letter=nl)
    fault_para.append([fault_lon, fault_lat])
    # ux for easting, uy for northing
    # 加上断层的长度以及走向
    ux = ux + math.sin(fault_strike * math.pi / 180) * 1000 * fault_length
    uy = uy + math.cos(fault_strike * math.pi / 180) * 1000 * fault_length
    # start to build fault one by one
    lat, lon = pSAR.utm_conversion.to_latlon(ux, uy, a, b, northern=None)
    fault_para.append([lon, lat])
    return fault_para

def buffer_generate(fault_para, length):
    nl = pSAR.utm_conversion.latitude_to_zone_letter(fault_para[0,1])
    zone = pSAR.utm_conversion.latlon_to_zone_number(fault_para[0,1], fault_para[0,0])
    ux,uy = pSAR.utm_conversion.lltoutm(fault_para[:, 1], fault_para[:, 0], force_zone_number=zone, force_zone_letter=nl)
    fault_line = [[ux[0],uy[0]],[ux[1],uy[1]]]
    line = LineString(fault_line)
    buffer = line.buffer(length)
    x2, y2 = buffer.boundary.xy
    lat,lon = pSAR.utm_conversion.utmtoll(x2,y2,force_zone_number=zone,force_zone_letter=nl)
    buffer_build = []
    for i in range(len(lat)):
        buffer_build.append([lon[i],lat[i]])
    return buffer_build

if __name__ == '__main__':
    print('This is for time trans')