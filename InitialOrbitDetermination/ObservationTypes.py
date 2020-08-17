# -*- coding: utf-8 -*-
"""
Created on Sun May 31 16:27:01 2020

@author: quamo
"""


import numpy as np
from org.hipparchus.geometry.euclidean.threed import Vector3D

# class AnglesOnlyObs ():
    
    
#     def __init__(self, date, RA, DE, RAunit, DEunit, site):
#         self.date_Greg = date
#         self.date_JD = self.Convert_Greg_to_JD ()
#         self.RAunit = RAunit
#         self.DEunit = DEunit
#         self.site = (np.radians(site[0]), np.radians(site[1]), site[2])
        
#         self.Get_r_site_ECEF()
        
        
#         if self.RAunit == "HH:MM:SS":
#             self.RA = RA
#             self.HHMMSS_to_deg()
#         elif self.RAunit == "deg":
#             self.RA = float(RA)
        
#         if self.DEunit == "DD:MM:SS":
#             self.DE = DE
#             self.DDMMSS_to_deg()
#         elif self.DEunit == "deg":
#             self.DE = float(DE)
        
#     def __str__(self):
#         return ("Observation at {} JD, RA {}, DE {}, r_site_ECEF {}".format(
#             self.date_JD, self.RA, self.DE, self.r_site_ECEF/1000))
        
#     def HHMMSS_to_deg (self):
        
#         self.RA =  360/ 86400 * (float(self.RA[0:2])*3600 + float(self.RA[2:4])*60 + float(self.RA[4::])) 
#         print(self.RA)
        
        
#     def DDMMSS_to_deg (self):
        
#         self.DE = float(self.DE[0:3]) + float(self.DE[0] + self.DE[3:5])/60 + float(self.DE[0] + self.DE[5::])/3600
#         print(self.DE)
        
    
#     def Convert_Greg_to_JD (self):
#          return TimeConversions.Gregorian_to_JD (self.date_Greg)
    

        
#     def Get_r_site_ECEF (self):
        
#         self.r_site_ECEF = CC.LatLongAlt_to_rSiteECEF(self.site)
        
        
            
        
def Get_LineOfSight_UnitVector(RA, DE):
  
    """
    Calculates the line of sight versor given the Declination and Right Ascension
    values for that observation  
    Returns a numpy array for each of these
    """
  
    # L = np.mat(np.array([[np.cos(DE) * np.cos(RA)],
    #                      [np.cos(DE) * np.sin(RA)],
    #                      [np.sin(DE)]]))  
    
    # due to the Java Python interface the Vector3D class will create problems if 
    # the float types are not explicitly cast
    L = Vector3D(float(np.cos(DE) * np.cos(RA)),
                 float(np.cos(DE) * np.sin(RA)),
                 float(np.sin(DE))) 
    
    
    return L    
        
        
