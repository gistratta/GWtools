import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

"""
   
   Author: Giulia Stratta
   
   Purpose: read the GW skymap, compute the sky coordinates
   of maximum probability point in the map, define a grid 
   of FOV to be visualize on the skymap
   
   Usage: 
   ipython
   run findprobatradec
   
"""
print ' '
print ' Please, type here skymap file with full path: '
print ' '
skymapfile = raw_input(' skymap file: ')
hpx = hp.read_map(skymapfile)
# skymap number of pixels
npix = len(hpx)
# nside=lateral resolution of HEALPix map (nside=2^x con x = max order)
nside = hp.npix2nside(npix)
deg2perpix = hp.nside2pixarea(nside, degrees=True)
deg2perpix

"""
# plot skymap
#moll = hp.mollview(hpx)
#plt.plot(moll)
#plt.show()

# find pixel that contain RA, Dec
ra = 2.865000
dec = -64.27300
theta = 0.5 * np.pi - np.deg2rad(dec)
phi = np.deg2rad(ra)
ipix = hp.ang2pix(nside, theta, phi)
ipix
"""

# find highest probability pixel (ipix_max). What is the probability inside it?
ipix_max = np.argmax(hpx)
ipix_max 
hpx[ipix_max]
print ' '
print ' highest prob. value = ',hpx[ipix_max] 

# Where is the highest probability pixel on the sky? Use pix2ang.
theta, phi = hp.pix2ang(nside, ipix_max)
ra = np.rad2deg(phi)
dec = np.rad2deg(0.5 * np.pi - theta)
print ' highest prob coord. ra,dec = ',ra, dec
print '  '



# Delta RA(deg) and Delta Dec(deg) of FoV
raw=3.0
decw=3.0

# area to cover centered on max prob. pix (sides are multiple of raw and decw)
ampRA=6.
ranRA=int(ampRA/raw)
ampDec=10.
ranDec=int(ampDec/decw)

k=0
for i in range(-ranRA,ranRA):
  for j in range(-ranDec,ranDec):

	dec1=(dec+j*decw)+decw/2.
	dec2=(dec+j*decw)-decw/2.
	decw0=dec2-dec1
	dec0=(dec2+dec1)/2.
	raw01=raw/np.cos(dec1*np.pi/180.)
	raw02=raw/np.cos(dec2*np.pi/180.)
	raw0=raw/np.cos(dec0*np.pi/180.)
	ra0=ra+i*raw0
	ra11=(ra+i*raw01)-raw01/2.
	ra12=(ra+i*raw01)+raw01/2.
	ra21=(ra+i*raw02)-raw02/2.
	ra22=(ra+i*raw02)+raw02/2.
#	print 'ICRS;POLYGON ',ra11,dec1,ra12,dec1,ra22,dec2,ra21,dec2,' # text = {prob.',prob,'}'
	k=k+1
	theta1 = 0.5 * np.pi - np.deg2rad(dec1)
	theta2 = 0.5 * np.pi - np.deg2rad(dec2)
	phi11 = np.deg2rad(ra11)
	phi12 = np.deg2rad(ra12)
	phi21 = np.deg2rad(ra21)
	phi22 = np.deg2rad(ra22)
	xyz11 = hp.ang2vec(theta1, phi11)
	xyz12 = hp.ang2vec(theta1, phi12)
	xyz22 = hp.ang2vec(theta2, phi22)
	xyz21 = hp.ang2vec(theta2, phi21)
	xyz=[(xyz11),(xyz12),(xyz22),(xyz21)]
	ipix_poly = hp.query_polygon(nside, xyz)
	pro=hpx[ipix_poly].sum()
#	print pro
	print 'ICRS;POLYGON ',ra11,dec1,ra12,dec1,ra22,dec2,ra21,dec2,' # text = {prob.,ra,dec',pro,ra0,dec0,'}'

#	xwidth1=xyz12-xyz11
#	xwidth2=xyz22-xyz21
#	ywidth1=xyz22-xyz11
#	ywidth2=xyz21-xyz12
#	print 'ICRS;BOX ',ra0,dec0,xwidth1,ywidth1,xwidth2,ywidth2,45,' # text = {prob.,ra,dec',pro,ra0,dec0,'}'



"""
    PART TO BE TEST...
    
    
table=np.array(16,10)
table=np.concatenate(['ICRS;POLYGON ',raw11,dec1,raw12,dec1,raw22,dec2,raw21,dec2,' # text = {source',i,j,'}'],axis=1)
table=table[~np.isnan(table).any(axis=1)]
np.savetxt('skycoord.reg', table)


print '  '
print 'Lista dei centroidi e delle probabilita '

for i in range(-2,2):
  for j in range(-2,2):
	dec1=(dec+j*decw)+decw/2.
	dec2=(dec+j*decw)-decw/2.
	dec0=(dec2+dec1)/2.
	raw0=(raw/2)*np.cos(dec0*np.pi/180.)
	ra0=ra+i*raw0
	print ra0,dec0,pro
	

#pixels that are contained within the circle ra, dec ,radius
# RA, Dec, and radius of circle in degrees (use: ipix_disc)
ra = 2.865000
dec = -6.427300
radiusdeg = 3

# Spherical polar coordinates and radius of circle in radians
theta = 0.5 * np.pi - np.deg2rad(dec)
phi = np.deg2rad(ra)
radius = np.deg2rad(radiusdeg)

# Cartesian coordinates of center of circle
xyz = hp.ang2vec(theta, phi)

# Array of indices of pixels inside circle
ipix_disc = hp.query_disc(nside, xyz, radius)

# Probability that source is within circle
probcirc=hpx[ipix_disc].sum()
probcirc
print ' '
print 'Prob that source is within circle of', radiusdeg, 'deg, centered at RA=', ra, 'and Dec.=',dec, 'is', probcirc

#Similarly, we can use the query_polygon function 
# qui definisco l'x,y,z  cartesiani dei vertici di un poligono
#xyz = [[-0.67, -0.41, -0.59],[-0.68, -0.40,-0.60],[-0.69,-0.39,-0.60],[-0.70,-0.40,-0.58]]
# e trovo i pixel contenuti nel poligono
#ipix_poly=hp.query_polygon(nside, xyz)
# Probability that source is within polygon
#hpx[ipix_poly].sum()

"""


