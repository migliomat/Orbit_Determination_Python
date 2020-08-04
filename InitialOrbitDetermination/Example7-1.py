# -*- coding: utf-8 -*-
"""
Created on Thu May 28 13:18:37 2020

@author: quamo
"""


# Testing the SITE TRACK ALGORITHM SO FAR

from SITE_TRACK import SITE_TRACK


Site_gdLat = 39.007 #deg
Site_Long  = -104.883 #deg
h_ell = 2.187 #km
rho = 604.68 #km
beta = 205.6 #deg
el = 30.7 #deg
rhoDot = 2.08 #km/s
betaDot = 0.15 #deg/s
elDot = 0.17 #deg/s


results = SITE_TRACK (Site_gdLat, Site_Long, h_ell, rho, beta, el, rhoDot, betaDot, elDot)

print(f"r = {results[0]} \n v = {results[1]}")