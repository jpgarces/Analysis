#!/usr/bin/python
 
 
######
# PYTHON 2.7 VERSION!
# This is intended to be as simple as possible, in order to share with others
# With this script you can make a moment 0 and PV plot
# Specifically, it is created for ALMA Data (but probably works for others...)
# You can do many more cool things, see Plunkett et al. 2015 Nature paper, Figure 3
########!!!!!########
#TO DO:
# 1) Try ploting separately the pv diagram and the moment 0 map. If that can be done, then te error asociates with fig3.show() comes from the creation of a figure variable. 
# 2) Try solving this error by finding a new way to plot both plots together. Look into matplotlib documentation.
# 3) Enhance the script in any possible way. Make it more friendly.
# 4) Modify the script to plot moment map or/and pv diagram.
########!!!!!########

import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
#plt.ion()
#plt.ioff()
from astropy.io import fits #for inspecting Fits datacubes
import numpy as np  #for mathy calculations; you need an up-to-date numpy
import aplpy ##JP## I get a WARNING here
from spectral_cube import SpectralCube #Get this here: https://spectral-cube.readthedocs.io/en/latest/installing.html
from matplotlib import lines
 
from astropy import wcs
from astropy import coordinates
from astropy import units as u
from astropy.coordinates import SkyCoord
 
import os
import sys
 
from pvextractor.gui import PVSlicer #Get this here: https://pypi.org/project/pvextractor/
from pvextractor import extract_pv_slice
from pvextractor.geometry import Path
import matplotlib.ticker as ticker

##Needed documentation##
#https://aplpy.readthedocs.io/en/stable/fitsfigure/quick_reference.html
#https://matplotlib.org/3.1.1/tutorials/index.html#advanced
#https://docs.astropy.org/en/stable/
#https://spectral-cube.readthedocs.io/en/latest/
#https://docs.astropy.org/en/stable/wcs/
#https://pvextractor.readthedocs.io/en/latest/programmatic.html#defining-a-path
writechan = False ## CHOOSE WHETHER TO WRITE IMAGE TO A FILE (True), OR DISPLAY DIRECTLY (False).
#Input Data
source=43
scoordra='18:30:04.4985'
scoorddec='-02:02:48.3805'
angle=132.7
majorax=3.2
 
 
### OPTIONAL:
### I use these to adjust things like axis labels, text size, etc.
### standard_setup for a map; pv_setup for a pv plot...
def standard_setup(sp): 
  sp.set_frame_color('black')
  sp.set_tick_labels_font(size='14')
  sp.set_axis_labels_font(size='16')
  sp.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss') #axes in sexagesimal
  sp.set_tick_color('black')
  #sp.set_system_latex(True) #use latex math symbols in plots
def pv_setup(sp):
  sp.set_frame_color('black')
  sp.set_tick_labels_font(size='14')
  sp.set_axis_labels_font(size='16')
  #sp.set_tick_labels_format(xformat='mm:ss',yformat='dd:mm:ss') #axes in sexagesimal
  sp.set_tick_color('black')
  #sp.set_system_latex(True) #use latex math symbols in plots
 
### Set some source parameters
vc=8.0*u.km/u.s #central velocity
 
### Need to make moments in CASA (or using Spectral Cube)
 
filename='gasC18O_FullMapNew.corrected.fits' ## Works for this file specifically, or you can provide your own *.fits file
## This is based on the data for CARMA-7 (serps45 from Plunkett et al. 2018)
## Specifically taking a cut along the disk axis (perpendicular to the outflow)
cube0 = SpectralCube.read(filename)
 
outname='pv_testS'+str(source)+'.eps'

fig3 = plt.figure(figsize=(9,6))
fig3.clf() #(.clf() ?) clear figure (i case you run it again)
#Plot over figures
leftpanel = [0.12,0.1,0.4,0.8] # (xo,yo,width,height)
rightpanel = [0.59,0.1,0.4,0.8] 
## Set up your map
ch0='' #or you can specify
chf='' #or you can specify
if ch0 != '' and chf !='': 
  cube = cube0.subcube(zlo=ch0,zhi=chf)
else: 
  cube=cube0
  ch0 = 0
  chf = cube.header['NAXIS3']
 
c = SkyCoord(scoordra,scoorddec, frame='icrs',unit=(u.hourangle,u.deg)) ##JP## to localize the source 
RAcoord=c.ra
Deccoord=c.dec
xcoord=[0,cube.shape[2]] # Limits of each axis. cube.shape() -> (chan,y,x)
ycoord=[0,cube.shape[1]]
x_w,y_w=[RAcoord.value],[Deccoord.value] # x and y source coord in word coordinates (deg)
 
#### Set the HALFlength and angle of the PV cut
halflength = majorax/2.*u.arcsec  # length of PV cut 
npix = round(halflength.to(u.deg) / (cube.header['cdelt2']*u.deg)) #calculate equivalent number of pixels. cdelt2->interval along yaxis. i.e. value of each pixel length in degrees.
 
pa = angle*u.deg #Deg: degrees away from vertical.
#pa0 = pa-180.*u.deg 
 
# Calculate the x,y coordinates on either end of the PV cut
upperx = RAcoord+halflength.to(u.deg)*np.sin(pa.to(u.rad)) #horizontal and vertical projections obtained with sin() and cos(). I think probably upperx is lowerx in case we use pa instead of pa0: sin(a+b)=sin(a)cos(b)+cos(a)sin(b)-> sin(180+b)=-sin(b)
lowerx = RAcoord-halflength.to(u.deg)*np.sin(pa.to(u.rad))
uppery = Deccoord+halflength.to(u.deg)*np.cos(pa.to(u.rad))
lowery = Deccoord-halflength.to(u.deg)*np.cos(pa.to(u.rad))
endpoints_deg = np.array([[upperx.value,uppery.value],[lowerx.value,lowery.value]]) #pv cut endpoints in deg
 
 
################################
#### MAKE A PV CUBE
################################
 
pvcube = cube.subcube(xlo=xcoord[0],xhi=xcoord[1],ylo=ycoord[0],yhi=ycoord[1],zlo=ch0,zhi=chf+1) 
hdu = pvcube
w = wcs.WCS(hdu.header) #world coordinate system. Translating info to tell python to look at the header and translate info from pixel to world system coord (degrees).
endpoints_pix = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix([[x,y] for x,y in endpoints_deg], 0) #from the image in world coord system (w) translate the endpoints cordinates in degrees (endpoints_deg) to pixels (endpoints_pix). For that we need to know how pixel<->world coordinates are related (w.sub([wcs.WCSSUB_CELESTIAL]))
endpoints_list=endpoints_pix.tolist() #listed as pixels
path=Path(endpoints_list,width=6) #3 pixels, is approximately 1 beam. Averaging along 6 pixels perpendicular to path (pv-cut). Maybe try to include as much of the source as possible. Adele used one beamsize at each side of the pv-cut.
pv = extract_pv_slice(hdu,path,wcs=w) 
 
pv_vel=pvcube.spectral_axis #chanel velocity values
minvel = pv_vel.min().value
maxvel = pv_vel.max().value
 
pvstd = pv.data[120:160,90:120].std() #calculate the RMS of the PV plot, to later make contours. Information is in pixels
pvstd2 = pv.data[620:680,90:120].std() #fix this to get region inside image
 
pvrms = np.nanmean([pvstd,pvstd2])
 
################################
#### PANEL (b): PV PLOT
################################

#fig3.show() # blank figure. We plot over this figure.
F4 = aplpy.FITSFigure(pv,dimensions=[1,0],figure=fig3,subplot=(rightpanel)) #pvplot. [1,0] flips the axis when ploting.
pv_setup(F4) #this makes the axes nicer (see the definition above)
F4.show_grayscale(aspect='auto',pmin=.5,pmax=99.75,stretch='linear',invert=True)
#F4.show_contour(colors='cyan',levels=10,linewidths=1.5,zorder='8') # does not work as intended

F4.refresh() #Adopt my new axis that I just created
print([x.get_text()[1:-1] for x in F4._ax1.get_yticklabels()]) ##?## ._ax1? axis 1 of the figure F4
 
#for intervals that are in simpler increments of arcseconds
#this is a ridiculous workaround for prettier axes...
 
xlim=F4._ax1.get_xlim() 
ylim=F4._ax1.get_ylim()
xlim_lo = minvel #Velocity in x-axis
xlim_hi = maxvel
ylim_lo=(ylim[0])*pv.header['CDELT1']+pv.header['CRVAL1']-pv.header['CDELT1'] ##?## range of y What are these keys and why to add them like that
ylim_hi=(ylim[1])*pv.header['CDELT1']+pv.header['CRVAL1']-pv.header['CDELT1']
 
inter = 1*u.arcsec #interval you want the arcsec to be marked on (i.e. every 2 arcsec)
inter_px = inter.to(u.deg).value/(ylim_hi-ylim_lo)*(ylim[1]-ylim[0]) #transforming the interval info in arcsecs to pixels. (regla de tres) 
ntick = (int((ylim_hi-ylim_lo)/inter.to(u.deg).value)+2) #number of ticks corresponding to the chosen interval.
print([x.get_text() for x in F4._ax1.get_yticklabels()])
F4.refresh() 
#adjust minor ticks inn between major
F4._ax1.yaxis.set_major_locator(ticker.MultipleLocator(inter_px)) # tick and label
F4._ax1.yaxis.set_minor_locator(ticker.MultipleLocator(inter_px/inter.value)) #ticks in between labels
newticks = np.arange(0,inter.value*ntick,inter.value)-halflength.value #Substract halflength to make the coords correspond to how far it is from the center of the source.
F4._ax1.set_yticklabels(newticks) #--> all yticks go to 0 ##?## Why? And with what purpose?
F4.refresh()
print([x.get_text()[:-2] for x in F4._ax1.get_yticklabels()])

F4._ax1.set_xticklabels([str(float(x.get_text()[:])/1000.)[:-2] for x in F4._ax1.get_xticklabels()]) #setting x tick labels in km/s. ##?## Note that here I made the change from F4._ax1.set_xticklabels([str(float(x.get_text()[1:-1])/1000.)[:-2] for x in F4._ax1.get_xticklabels()])
F4._ax1.set_yticklabels([x.get_text()[:] for x in F4._ax1.get_yticklabels()]) ##NOTE## Here edit the indexing of x.get_text()[:-2]. #with [:-2] we achieve to do '12.0 -> '12'
F4._ax1.set_ylabel('Offset (arcsec)') 
F4._ax1.set_xlabel("$V_{LSR} (\mathrm{km\ s}^{-1})$")
print([x.get_text() for x in F4._ax1.get_yticklabels()])
F4.refresh()
 
F4._ax2.yaxis.set_major_locator(ticker.MultipleLocator(inter_px)) #ax2? 
F4._ax2.yaxis.set_minor_locator(ticker.MultipleLocator(inter_px/inter.value))
F4._ax1.yaxis.labelpad = -9 ##?## negative number? #Spacing in points from the axes bounding box including ticks and tick labels. Negative makes it closer
 
F4.refresh()
 
 
 
################################
#### PANEL (a): Moment 0 map
################################
 
## FIRST MAKE MOMENT MAP:
moment_0 = cube0.moment(order=0)
 
## NEXT PLOT (next to the PV plot):
f2 = aplpy.FITSFigure(moment_0.hdu,dimensions=[0,1],figure=fig3,subplot=(leftpanel))
standard_setup(f2) #this makes the axes nicer
f2._ax1.yaxis.labelpad = -20 ## Move label closer to figure
## you can set the minimum and max values for the grayscale plot (you choose these values)
mom0_vmin = np.nanmin(moment_0).value*0.4 #nanmin chooses the minimun, while ignorning NaN (Not a Number) values
mom0_vmax = np.nanmax(moment_0).value*0.6 #Shorten colorscale to make it more sensible for a higher number of pixels.
f2.show_grayscale(vmin=mom0_vmin,vmax=mom0_vmax,invert=True)
 
## you can recenter the plot to your chosen coordinate
xx=np.abs(cube0.header['NAXIS1']*cube0.header['cdelt1'])#transforming image width from pixel to degree
yy=np.abs(cube0.header['NAXIS2']*cube0.header['cdelt2'])#transforming image height from pixel to degree
mapwidth=5.*u.arcsec
mapheight=5.*u.arcsec
f2.recenter(RAcoord,Deccoord,width=mapwidth.to(u.deg).value,height=mapheight.to(u.deg).value)
f2.show_markers(RAcoord, Deccoord,marker="*",c='None',edgecolors="magenta",s=300,linewidths=2.0,zorder='9')
 
#### TO SHOW PV CUT
f2.show_lines([endpoints_deg.T],color='DarkMagenta',linewidths=2.0,linestyle='-') # T stands for transverse
 
#### TO Show Beam size (NOTE, SOEMTIMES THIS IS NOT INCLUDED IN YOUR DATA, YOU CAN DEFINE YOURSELF)
##JP## I think header is missing beamsize information. By inspecting other scripts I found that for this data set bmaj=1.043*u.arcsec.to(u.deg), bmin=0.7046*u.arcsec.to(u.deg), bpa=76.299. Got this from mom0 scripts.
bmaj=1.043*u.arcsec.to(u.deg)
bmin=0.7046*u.arcsec.to(u.deg)
bpa=76.299
f2.show_beam(major=bmaj,minor=bmin,angle=bpa,fill=True,color='DeepSkyBlue') 
f2.add_label(0.7, 0.93, '(a) C18O map', relative=True, color="Black", size=20)
 
################################
## FINAL ANNOTATIONS
################################
# HORIZONTAL LINE TO SHOW resolution in PV Plot. In PV Diagrams we have only one axis representing position and we need position to define resolution.
res_elem = 0.15*u.arcsec ##?## Where do I get this?? #resolution element, or beam size 
x1_res=res_elem.to(u.deg).value 
x2_res=x1_res+res_elem.to(u.deg).value 
y1_res=pvcube.spectral_axis.min().value
y2_res=pvcube.spectral_axis.max().value 
 
#### BEAMSIZE
F4.show_arrows(30080.,3./3600.,0.0,0.5*res_elem.to(u.deg).value,color='DeepSkyBlue',linestyles='-',head_length=0.1,head_width=5) 
F4.show_arrows(30080.,3./3600.,0.0,-0.5*res_elem.to(u.deg).value,color='DeepSkyBlue',linestyles='-',head_length=0.1,head_width=5)
F4.add_label(23000.,3./3600.,'0.14$^{\prime\prime}$ \n 60 AU',layer='label',color='DeepSkyBlue',size=16, verticalalignment='center')
 
# Set up a second Axis on top of pv plot, to show annotations (if necessary) ##?## invisible axes just to make annotations
ax2 = fig3.add_axes(rightpanel) 
ax2.set_xlim(xlim_hi,xlim_lo)
ax2.set_ylim(ylim_lo,ylim_hi)
 
ax2.set_axis_off()
 
# VERTICAL (dashed) LINE TO SHOW v_c in PV Plot
hx1 = 0.
hx2 = 2.*halflength.to(u.deg).value
hy1 = vc.to(u.m/u.s).value
ax2.vlines(hy1,hx1,hx2,colors='green',linestyles='dashed',label='$v_c$',zorder='3',lw=2) ##JP## I am getting vertical line shifted
ax2.annotate('$v_c$',xy=(hy1,x1_res*2.0),xycoords='data',color='green',size='16')
 
# HORIZONTAL (dashed) LINE TO SHOW star position in PV Plot. Separate position to left and right from the source along the position line. 
hx3 = xlim_lo
hx4 = xlim_hi
hy3 = halflength.to(u.deg).value
ax2.hlines(hy3,hx3,hx4,colors='magenta',linestyles='dashed',label='protostar',zorder='3',lw=2)
ax2.annotate('protostar',xy=(xlim_hi,hy3-2*res_elem.to(u.deg).value),xycoords='data',color='magenta',size='20',verticalalignment='top')
ax2.annotate('(b) CO, \n PV diag.',xy=(xlim_hi,hy3+8*res_elem.to(u.deg).value),xycoords='data',color='Black',size='20',verticalalignment='top')
#### NOTE, You have to move these around a bit depending on your target. ##NOTE## Feel free to modify text position.
 
 
###################################
###################################
###################################
###################################
#plt.show()
if writechan: 
  plt.savefig(outname)
if not writechan:
  print 'print to screen'
  fig3.show() ##?## 

#plt.close(fig3)
###################################
###################################
###################################
###################################

