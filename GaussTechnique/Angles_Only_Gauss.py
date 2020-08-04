# -*- coding: utf-8 -*-
"""
Created on Wed May 27 15:31:10 2020

This will be the Gauss Angles only routine used to produce an initial orbit determination
starting from angles only observations


@author: quamo
"""



#initialize orekit and JVM
import orekit
orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime
setup_orekit_curdir()

# Math
from math import radians, pi, degrees
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Geometry and Frames
from org.hipparchus.geometry.euclidean.threed import Vector3D, SphericalCoordinates
from org.orekit.data import DataProvidersManager, ZipJarCrawler
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory

# Time and Coordinates
from org.orekit.time import TimeScalesFactory, AbsoluteDate, DateComponents, TimeComponents
from org.orekit.utils import IERSConventions, Constants, PVCoordinates, PVCoordinatesProvider, AbsolutePVCoordinates

# Orbits
from org.orekit.orbits import KeplerianOrbit, CartesianOrbit, PositionAngle, OrbitType

# Propagators
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from java.io import File

# Measurements
from org.orekit.estimation.measurements import AngularRaDec

# Ground stations
from org.orekit.estimation.measurements import GroundStation
from org.orekit.estimation.measurements import ObservableSatellite

# CSV
import csv



##############################################################################


# Create Earth and Earth Frames
earthFrame = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                       Constants.WGS84_EARTH_FLATTENING,
                                       earthFrame)


# CREATE THE GROUND STATION
# First create the Topocentric Frame
station_coord = GeodeticPoint(radians(40), radians(-110), 2000.0)
station_frame = TopocentricFrame(earth, station_coord, 'MyGndStation')

# Now create the Ground station itself and the satellite object
GndStation = GroundStation(station_frame)
Sat = ObservableSatellite(1)

# SET THE DAY OF THE OBSERVATIONS
yr_mo_day = (2012, 8, 20)


# OPEN THE CSV FILE CONTAINING THE OBSERVATIONS
data = open('Observations.csv')

csv_data = csv.reader(data)

data_lines = list(csv_data)

# LOOP OVER THE FILE AND CREATE THE OBSERVATION OBJECTS
utc = TimeScalesFactory.getUTC()
dates = []
obs = []
i = 0
for line in data_lines[1::]:
    print(line)
    dates.append(AbsoluteDate(yr_mo_day[0], yr_mo_day[1], yr_mo_day[2],
                              int(line[1]), int(line[2]), float(line[3]), utc))
    
    obs.append(AngularRaDec(GndStation, earthFrame, dates[i], [radians(float(line[4])), 
                                                               radians(float(line[5]))], [0.0, 0.0], [0.0, 0.0], Sat))
    i += 1




# GAUSS ORBIT DETERMINATION ROUTINE
# Here's the actual orbit determination algorithm


tau1 = obs[0].getDate().durationFrom(obs[1].getDate())
tau3 = obs[2].getDate().durationFrom(obs[1].getDate())


a1 = tau3 / (tau3 - tau1)

a1_u = tau3*((tau3 - tau1)**2 - tau3**2) / (6*(tau3 - tau1))

a3 = - tau1 / (tau3 - tau1)

a3_u = - tau1 * ((tau3 - tau1)**2 - tau1**2) / (6*(tau3 - tau1))




































# def AnglesOnlyGauss(obs_table, obs_station):
    
#     obs = []
#     for item in obs_table:
#         obs.append(AnglesOnlyObs(item[0:6], item[6], item[7], item[8], item[9], obs_station))
     
#     JD1 = obs[0].date_JD
#     JD2 = obs[1].date_JD
#     JD3 = obs[2].date_JD
    
#     tau1 = JD1 - JD2
#     tau3 = JD3 - JD2
    

    
#     L1 = obs[0].Get_LineOfSight_UnitVector()
#     L2 = obs[1].Get_LineOfSight_UnitVector()
#     L3 = obs[2].Get_LineOfSight_UnitVector()
    
#     L = np.hstack((L1, L2, L3)) # horizontal concatenation to form L matrix
    
#     L_inv = np.linalg.inv(L)
    
#     print(L_inv)
    