from sympy import Plane, Point, Point3D
import numpy as np
import pandas as pd
import math


def strdip2plane(point, strike, dip):
    if dip < 0 or dip > 90:
        raise ('dip need to [0,90]')
    if strike < 0 or strike >= 360:
        raise ('dip need to [0,360)')
    dip = math.radians(dip)
    strike = math.radians(strike)
    n = (math.cos(math.pi / 2 - dip) * math.sin(math.pi / 2 + strike),
         math.cos(math.pi / 2 - dip) * math.cos(math.pi / 2 + strike),
         math.sin(math.pi / 2 - dip))
    P1 = Plane(Point3D(point[0], point[1], point[2]), normal_vector=n)
    return P1


fault_name = '../data/turkey_7_8_fault.txt'
data_name = '../data/data_processing.csv'
fault = pd.read_csv(fault_name)
data = pd.read_csv(data_name)
data = np.array(data)
print(data)
print(fault)
fault_point = [fault.iloc[0, 1], fault.iloc[0, 2], 0]
p = Point3D(data[0, 2], data[0, 3], data[0, 4])
p1 = strdip2plane(fault_point, fault.iloc[0, 5], fault.iloc[0, 6])

projectionPoint = p1.projection(p)
print(p1.equation())
