# Orbit_Determination_Python

A repository containing scripts developing initial orbit determination routines as taken from Vallado. 

SITE-TRACK is an orbit determination routine which starting from range, azimuth, elevation and their respective range measurements provides an initial orbit determination (position and velocity of a satellite).

Angles_Only_Gauss uses Gauss angles only method, which provides an initial orbit determination starting from 3 lines of sight vectors and no range or rate measurements. Orekit is used as low level library for handling time and coordinate frames. 
