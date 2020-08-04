# -*- coding: utf-8 -*-
"""
Created on Wed May 27 15:50:23 2020


Site Track application as proposed by Vallado in Initial Orbit Determination

@author: quamo
"""


import math
import numpy as np


R_Earth = 6378136.3 # meters
e_Earth = 0.081819221456


def SITE_TRACK (Site_gdLat, Site_Long, h_ell, rho, beta, el, rhoDot, betaDot, elDot):
    
    # Reassign the angular measurements as radians and range as meters
    Site_gdLat = np.radians(Site_gdLat) #rad
    Site_Long  = np.radians(Site_Long) #rad
    h_ell      = h_ell * 1000 # m
    rho        = rho * 1000 # m
    beta       = np.radians(beta) #radians
    el         = np.radians(el) #radians
    rhoDot     = rhoDot * 1000 #m/sec
    betaDot    = np.radians(betaDot) #rad/sec
    elDot      = np.radians(elDot) #rad/sec
    
    
    C_Earth = R_Earth / (math.sqrt(1-math.pow(e_Earth, 2)*(math.sin(Site_gdLat))**2))
    S_Earth = C_Earth*(1-e_Earth**2)
    
    r_delta = (C_Earth + h_ell)*math.cos(Site_gdLat)
    r_K     = (S_Earth + h_ell)*math.sin(Site_gdLat)
    
    r_site_ECEF = np.mat(np.array([[r_delta*np.cos(Site_Long)], [r_delta*np.sin(Site_Long)], [r_K ]]))
    
    rho_SEZ = np.mat(np.array([[-rho*np.cos(el)*np.cos(beta)],
                       [ rho*np.cos(el)*np.sin(beta)],
                       [ rho*np.sin(el)]]))
    
    
    rhoDot_SEZ = np.mat(np.array([[-rhoDot*np.cos(el)*np.cos(beta) + rho*np.sin(el)*np.cos(beta)*elDot + rho*np.cos(el)*np.sin(beta)*betaDot],
                          [ rhoDot*np.cos(el)*np.sin(beta) - rho*np.sin(el)*np.sin(beta)*elDot + rho*np.cos(el)*np.cos(beta)*betaDot],
                          [ rhoDot*np.sin(el) + rho*np.cos(el)*elDot]]))
    
    SEZ_2_ECEF = np.mat(np.array([[np.sin(Site_gdLat)*np.cos(Site_Long), -np.sin(Site_Long), np.cos(Site_gdLat)*np.cos(Site_Long)],
                          [np.sin(Site_gdLat)*np.sin(Site_Long),  np.cos(Site_Long), np.cos(Site_gdLat)*np.sin(Site_Long)],
                          [-np.cos(Site_gdLat),                   0,                 np.sin(Site_gdLat)                  ]]))
    
    
    rho_ECEF = SEZ_2_ECEF * rho_SEZ
    rhoDot_ECEF = SEZ_2_ECEF * rhoDot_SEZ
    
    r_ECEF = (r_site_ECEF + rho_ECEF) / 1000 # now transorm to km
    v_ECEF = rhoDot_ECEF / 1000 # now transform to km/s
    
    return r_ECEF, v_ECEF
    
    
    
    
    