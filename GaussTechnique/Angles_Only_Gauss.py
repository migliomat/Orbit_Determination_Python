# -*- coding: utf-8 -*-
"""
Created on Wed May 27 15:31:10 2020

This will be the Gauss Angles only routine used to produce an initial orbit determination
starting from angles only observations

Python 3.5 or above necessary


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

# Constants, Time and Coordinates
from org.orekit.time import TimeScalesFactory, AbsoluteDate, DateComponents, TimeComponents
from org.orekit.utils import IERSConventions, Constants, PVCoordinates, PVCoordinatesProvider, AbsolutePVCoordinates

mu = Constants.IERS2010_EARTH_MU

# Orbits
from org.orekit.orbits import KeplerianOrbit, CartesianOrbit, PositionAngle, OrbitType

# Propagators
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from java.io import File

# Measurements
from org.orekit.estimation.measurements import AngularRaDec
from Orbit_Determination_Python.InitialOrbitDetermination.ObservationTypes import Get_LineOfSight_UnitVector

# Ground stations
from org.orekit.estimation.measurements import GroundStation
from org.orekit.estimation.measurements import ObservableSatellite

# CSV
import csv



##############################################################################
##############################################################################
# Create the accessory objects and handle the observations
##############################################################################
##############################################################################

# Create the GCRF (ECI) and ITRF (ECEF) Frames
GCRF_Frame = FramesFactory.getGCRF()
J2000_Frame = FramesFactory.getEME2000()
ITRF_Frame = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                       Constants.WGS84_EARTH_FLATTENING,
                                       ITRF_Frame)


# CREATE THE GROUND STATION
# First create the Topocentric Frame
station_coord = GeodeticPoint(radians(40), radians(-110), 2000.0)
station_frame = TopocentricFrame(earth, station_coord, 'MyGndStation')

# Now create the Ground station itself and the satellite object
GndStation = GroundStation(station_frame)
Sat = ObservableSatellite(1) # create the observable satellite object, name it 1 as default

# SET THE DAY OF THE OBSERVATIONS
yr_mo_day = (2012, 8, 20)
utc = TimeScalesFactory.getUTC()

# OPEN THE CSV FILE CONTAINING THE OBSERVATIONS
data = open('Observations.csv')

csv_data = csv.reader(data)

data_lines = list(csv_data)

# LOOP OVER THE FILE AND CREATE THE OBSERVATION OBJECTS
# AngularRaDec objects hold Ra and Dec in this order

dates = []
obs = []

# I use the list method here to append rows and then transpose it, the next for
# loop displays a different method using np.arrays
LOSversors = []


i = 0
for line in data_lines[1::]:
    
    dates.append(AbsoluteDate(yr_mo_day[0], yr_mo_day[1], yr_mo_day[2],
                              int(line[1]), int(line[2]), float(line[3]), utc))
    
    obs.append(AngularRaDec(GndStation, ITRF_Frame, dates[i], [radians(float(line[4])), 
                                                               radians(float(line[5]))], [0.0, 0.0], [0.0, 0.0], Sat))
    
    LOSver = Get_LineOfSight_UnitVector(obs[i].getObservedValue()[0], obs[i].getObservedValue()[1])
    LOSver_column = LOSver.toArray()
    

    LOSversors.append(LOSver_column)
    
    i += 1

# Could not find a better way to get the vectors in as column vectors, will have to investigate further
LOSarray = np.array(LOSversors).T


##############################################################################
##############################################################################
# GAUSS ORBIT DETERMINATION ROUTINE
##############################################################################
##############################################################################


tau1 = obs[0].getDate().durationFrom(obs[1].getDate())
tau3 = obs[2].getDate().durationFrom(obs[1].getDate())


a1 = tau3 / (tau3 - tau1)

a1_u = tau3*((tau3 - tau1)**2 - tau3**2) / (6*(tau3 - tau1))

a3 = - tau1 / (tau3 - tau1)

a3_u = - tau1 * ((tau3 - tau1)**2 - tau1**2) / (6*(tau3 - tau1))

# Get the ground station position vector in GCRF (ECI) at each of the observation dates
# and insert it in a 3x3 matrix
r_site_GCRF_mat = np.empty(shape=[3,0])

for date in dates:
    
    GndStationPos_GCRF = station_frame.getPVCoordinates(date, GCRF_Frame).getPosition()
    
    Pos = np.array([[GndStationPos_GCRF.getX()], [GndStationPos_GCRF.getY()], 
                                     [GndStationPos_GCRF.getZ()]])
    
    r_site_GCRF_mat = np.append(r_site_GCRF_mat, Pos, axis=1 ) # append along the columns (axis=1)


# We now need to solve the linear system:
#
# LOSarray [[c1*rho1],[c2*rho2],[c3*rho3]] = r_site_GCRF_mat*[[-c1], [-c2], [-c3]]
#
# Our unknowns are the rho values, that is the range (or slant range) information to each observation
# We are specifically interested in the rho2 value

LOSarray_inv = np.linalg.inv(LOSarray)

M = LOSarray_inv @ r_site_GCRF_mat # @ is the matrix multiplication operator 

d1 = M[1][0]*a1 - M[1][1] + M[1][2]*a3 #given python indexing all indices must be reduced by 1 wrt the written eq
d2 = M[1][0]*a1_u + M[1][2]*a3_u

C = np.dot(LOSarray[:,1], r_site_GCRF_mat[:,1])

#### 8th degree polynomial in r2

P = list(range(9))

P[0] = 1 #8th
P[1] = 0 #7th
P[2] = -(d1**2 + 2*C*d1 + np.linalg.norm(r_site_GCRF_mat[:,1])**2) #6th
P[3] = 0 #5th
P[4] = 0 #4th
P[5] = -2*mu*(C*d2 + d1*d2) #3rds
P[6] = 0 #2nd
P[7] = 0 #1st
P[8] = -mu**2 * d2**2 #0th

roots = np.roots(P)

PosRoots = []
for item in roots:
    if np.isreal(item) and item > 0:
        PosRoots.append(item)

if len(PosRoots) == 1:
    r2 = PosRoots[0].real
else:
    pass

u = mu / r2**3

c1 = -(-a1 -a1_u*u)
c2 = -1
c3 = -(-a3 -a3_u*u)

B = M @ np.array([[-c1], [-c2], [-c3]])

A = np.identity(3)
A[0,0] = c1
A[1,1] = c2
A[2,2] = c3

slant_ranges = np.linalg.solve(A,B)

print(slant_ranges)

###############################################################################
###############################################################################
# END OF GAUSS ORBIT DETERMINATION ROUTINE
###############################################################################
###############################################################################
 


































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
    