# PYTHON 2.7 VERSION!
#USEFUL DOCUMENTATION:
#https://aplpy.readthedocs.io/en/stable/fitsfigure/quick_reference.html
#https://matplotlib.org/3.1.1/tutorials/index.html#advanced
#https://docs.astropy.org/en/stable/
#https://spectral-cube.readthedocs.io/en/latest/
#https://docs.astropy.org/en/stable/wcs/
#https://pvextractor.readthedocs.io/en/latest/programmatic.html#defining-a-path

#NOTE: FOR FULLY COMMENTED SCRIPT SEE /home/jp/Research/Data/pv_simplify_py2.py

#ALL LINES WITH '#MOD' NEXT TO IT HAVE MODIFICABLE QUANTITIES

#TO DO: 1) PV y-axis too many digits. Play with indexing to correct labels.
#	READ SECTION 4 MURILLO 2013!!!
#	2) Make map with elipses resulting from 2d-fit in casa (minorax,majax,angle)
#	3) Plot keplerian model lines on top. Plug mass estimates obtained (about three for each source: minorax, majorax, averax) and plot resulting curve on top of the pv-diagram. From this we add a criteria to show which mass value should be correct. Of course this works by first making the assumption that our source has keplerian rotation. See Skype chat Adele. Watch out for the units km<->arcs:
	#myvelocities=np.arange(4,14,0.1)
	#myradii=G*M*myvelocities**-2
	#plt.plot(myvel,myraddi(in arcs: from .units.km to.arcs))
#	4) Try to make show-contours work
#	5) FILTER SOURCES WITH NO INTERESTING (NO PROMETEDOR) PV-DIAGRAMS.
#	6) FOR THE INTERESTING SOURCES, ENHANCE THE PLOT. TRY DIFFERENT PV-CUT ANGLES AND AVERAGES.
#	7) ADD OTHER MODEL CURVES. READ THROUGH THESE LINES ADELE SENT AND MAKE THEM WORK HERE
#	8) ADD UNITS
#	9) CHANGE SHOW_MARKERS() TO SHOW_LINES(). WILL HAVE TO PLAY WITH ARRAY BRACKETS.
#	10) make a txt file with the possible masses organized in columns and automate this step to get info of M from there.
#	!!! 11) Plot keplerian model with several different masses and from that specify range of acceptable masses (expected mass values). See how it compares to the previously made calculations based on 2D-fitting. 
#	!!! 12) Change position angle to make plot look better. 
#	13) Try working other brilliant isolated source (simpler) and try out all this procedure.


#PACKAGES:
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
#plt.ion()
#plt.ioff()
from astropy.io import fits 
import numpy as np  
import aplpy ##JP## I get a WARNING here
from spectral_cube import SpectralCube 
from matplotlib import lines
 
from astropy import wcs
from astropy import coordinates
from astropy import units as u
from astropy.coordinates import SkyCoord
 
import os
import sys
 
from pvextractor.gui import PVSlicer 
from pvextractor import extract_pv_slice
from pvextractor.geometry import Path
import matplotlib.ticker as ticker

#PREVIOUS INFO SETUP:
writechan = False ## CHOOSE WHETHER TO WRITE IMAGE TO A FILE (True), OR DISPLAY DIRECTLY (False). #MOD
source=45 #MOD
vc=8.0*u.km/u.s #central velocity #MOD
#File:
filename='/home/juan/Desktop/Research/Data2.0/Serpens/Serpens_Feather_c18o.fits' #MOD
outname='pv_testS'+str(source)+'.eps'
#Plot over figures:
leftpanel = [0.12,0.1,0.38,0.8] # (xo,yo,width,height) #MOD
rightpanel = [0.61,0.1,0.38,0.8] #MOD
## Set up your map
cube0 = SpectralCube.read(filename)
ch0='' #MOD
chf='' #MOD
if ch0 != '' and chf !='': 
  cube = cube0.subcube(zlo=ch0,zhi=chf)
else: 
  cube=cube0
  ch0 = 0
  chf = cube.header['NAXIS3']

#Create figure
fig3 = plt.figure(figsize=(9,6)) #MOD
fig3.clf() 

#SOURCE COORDS:
filenametxt='/home/juan/Desktop/Research/Data/sourcecoordinates.txt'
ra,dec = np.loadtxt(filenametxt,unpack=True,usecols=(0,1),dtype='S')
coord=SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle,u.deg)) 
c = SkyCoord(coord[source-1].ra,coord[source-1].dec, frame='icrs',unit=(u.hourangle,u.deg)) 
#SOURCE PV INFO: THERE IS PROBABLY A MUCH SHORTER AND GENERAL WAY TO DO THIS
file2dfit='/home/juan/Desktop/Research/Data2.0/2DFit/'+str(source)+'/2DGssfit'+str(source)+'Log.txt'
file1=open(file2dfit,'r')
i=1
for line in file1:
	if i==27:
		line1=line
	if i==29:
		line2=line
	i+=1
file1.close()
for j in range(len(line1)):
	if line1[j]==':':
		k=j+1
		majax1=''
		while line1[k]!='+':
			majax1+=line1[k]
			k+=1
for j in range(len(line2)):
	if line2[j]==':':
		k=j+1
		angle1=''
		while line2[k]!='+':
			angle1+=line2[k]
			k+=1
dchan=1 #0.5 to shorten interval
amp=5
angle=float(angle1.strip())
majorax=amp*float(majax1.strip())

#PV-CUT:
halflength = majorax/2.*u.arcsec  #MOD
pa = angle*u.deg #Deg: degrees away from vertical.
avnumb=6 #number of pixels to average perpendicular to pv-cut #MOD
inter = 1*u.arcsec #interval you want the arcsec to be marked on #MOD

### OPTIONAL:
def standard_setup(sp): 
  sp.set_frame_color('black')
  sp.set_tick_labels_font(size='14')
  sp.set_axis_labels_font(size='16')
  sp.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss') 
  sp.set_tick_color('black')
  #sp.set_system_latex(True) #use latex math symbols in plots
def pv_setup(sp):
  sp.set_frame_color('black')
  sp.set_tick_labels_font(size='14')
  sp.set_axis_labels_font(size='16')
  #sp.set_tick_labels_format(xformat='mm:ss',yformat='dd:mm:ss') 
  sp.set_tick_color('black')
  #sp.set_system_latex(True)

#LOCATING SOURCE:
RAcoord=c.ra
Deccoord=c.dec
xcoord=[0,cube.shape[2]] 
ycoord=[0,cube.shape[1]]
x_w,y_w=[RAcoord.value],[Deccoord.value] 		
  
# Calculate the x,y coordinates on either end of the PV cut
upperx = RAcoord+halflength.to(u.deg)*np.sin(pa.to(u.rad))
lowerx = RAcoord-halflength.to(u.deg)*np.sin(pa.to(u.rad))
uppery = Deccoord+halflength.to(u.deg)*np.cos(pa.to(u.rad))
lowery = Deccoord-halflength.to(u.deg)*np.cos(pa.to(u.rad))
endpoints_deg = np.array([[upperx.value,uppery.value],[lowerx.value,lowery.value]]) 

################################
#### MAKE A PV CUBE
################################
pvcube = cube.subcube(xlo=xcoord[0],xhi=xcoord[1],ylo=ycoord[0],yhi=ycoord[1],zlo=ch0,zhi=chf+1) 
hdu = pvcube
w = wcs.WCS(hdu.header)
endpoints_pix = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix([[x,y] for x,y in endpoints_deg], 0) 
endpoints_list=endpoints_pix.tolist() 
path=Path(endpoints_list,width=avnumb) 
pv = extract_pv_slice(hdu,path,wcs=w) 
pv_vel=pvcube.spectral_axis 
minvel = pv_vel.min().value
maxvel = pv_vel.max().value
pvstd = pv.data[0:20,0:20].std() #choose a low signal region #MOD
pvstd2 = pv.data[50:70,170:190].std() #MOD
pvrms = np.nanmean([pvstd,pvstd2])
 
################################
#### PANEL (b): PV PLOT
################################
#PLOT:
F4 = aplpy.FITSFigure(pv,dimensions=[1,0],figure=fig3,subplot=(rightpanel)) 
pv_setup(F4)
F4.show_grayscale(aspect='auto',pmin=.5,pmax=99.75,stretch='linear',invert=True)
F4.show_contour(colors='orange',levels=5,dimensions=[1,0]) # does not work as intended
#MODIFY AXES:
F4.refresh() 
print([x.get_text()[1:-1] for x in F4._ax1.get_yticklabels()]) 
xlim=F4._ax1.get_xlim() 
ylim=F4._ax1.get_ylim()
xlim_lo = minvel 
xlim_hi = maxvel
ylim_lo=(ylim[0])*pv.header['CDELT1']+pv.header['CRVAL1']-pv.header['CDELT1']
ylim_hi=(ylim[1])*pv.header['CDELT1']+pv.header['CRVAL1']-pv.header['CDELT1']
npix = round(halflength.to(u.deg) / (cube.header['cdelt2']*u.deg))
inter_px = inter.to(u.deg).value/(ylim_hi-ylim_lo)*(ylim[1]-ylim[0])
ntick = (int((ylim_hi-ylim_lo)/inter.to(u.deg).value)+2) 
print([x.get_text() for x in F4._ax1.get_yticklabels()])
F4.refresh() 
F4._ax1.yaxis.set_major_locator(ticker.MultipleLocator(inter_px)) 
F4._ax1.yaxis.set_minor_locator(ticker.MultipleLocator(inter_px/inter.value)) 
newticks = np.arange(0,inter.value*ntick,inter.value)-halflength.value 
F4._ax1.set_yticklabels(newticks)
F4.refresh()
print([x.get_text()[:-2] for x in F4._ax1.get_yticklabels()])
F4._ax1.set_xticklabels([str(float(x.get_text()[:])/1000.*dchan)[:-2] for x in F4._ax1.get_xticklabels()])
F4._ax1.set_yticklabels([x.get_text()[:5] for x in F4._ax1.get_yticklabels()],size=10)##NOTE## Here edit the indexing of x.get_text()[:-2]. #with [:-2] we achieve to do '12.0 -> '12'
F4._ax1.set_ylabel('Offset (arcsec)',size=13) 
F4._ax1.set_xlabel("$V_{LSR} (\mathrm{km\ s}^{-1})$")
print([x.get_text() for x in F4._ax1.get_yticklabels()])
F4.refresh()
F4._ax2.yaxis.set_major_locator(ticker.MultipleLocator(inter_px)) 
F4._ax2.yaxis.set_minor_locator(ticker.MultipleLocator(inter_px/inter.value))
F4._ax1.yaxis.labelpad = -9 #Edit to modify spacing #MOD

#PLOT MODEL LINES
##KEPLERIAN ROTATION: (try with the different values of M obtained from the 2d-fit performed in CASA)
velvalues=np.arange(4,14,0.1)*1000*(u.m/u.s)-vc.to(u.m/u.s) #substract reference velocity
#velvalues=np.arange(4,6,0.1)*1000-vc.to(u.m/u.s).value
G=6.674e-11*u.m**3*u.s/u.kg
#Mlist= [8.3*u.Msun,10.1*u.Msun]
Mlist= [0.5*u.Msun,2.0*u.Msun]
colors=['cyan','red']
#In the future make a txt file with the possible masses organized in columns and automate this step to get info of M from there.
l=0
for M in Mlist:
	rvalues=G*M.to(u.kg)*velvalues**-2 ##?## we use gas mass, but what happens to dust mass?
	#(from http://icc.dur.ac.uk/~tt/Lectures/Galaxies/TeX/lec/node41.html)
	rvalarcs=rvalues.value*206265/(436*3.086*10**16)*u.arcsec#from km to arcsec (using smaal angle aprox). Distance to serpens south 436 parsecs
	F4.show_lines([np.array([velvalues.value+vc.to(u.m/u.s).value,rvalarcs.to(u.deg).value+halflength.to(u.deg).value])],color=colors[l],linewidths=2.0,linestyles='-')
	#since we also have negative radii we als have to add
	F4.show_lines([np.array([velvalues.value+vc.to(u.m/u.s).value,-rvalarcs.to(u.deg).value+halflength.to(u.deg).value])],color=colors[l],linewidths=2.0,linestyles='-')
	l+=1

######ADD#########
'''if plotrotation: 
  #### ROTATION LINES
  #### SEE Murillo et al. 2013

  M = 0.5*u.Msun

  R_as = np.arange(0,3,0.1)*u.arcsec
  R_pc = R_as.value*dd/206265.

  v_ff=np.sqrt(2*const.G*M/R_pc) #### EQ 2, free-fall
  v_rot=np.sqrt(const.G*M/R_pc) #### EQ 3, pure Keplerian rotation
  C = 2.
  R_cam = 20.*u.AU
  v_am=C*np.sqrt(const.G*M/R_cam)*R_cam/R_pc #### EQ 4, conserved angular momentum
  R_cinf = 200.*u.AU
  v_inf=np.sqrt(const.G*M/R_cinf)*R_cinf/R_pc #### EQ 5, infall + Keplerian

  #### PURE FREEFALL
  v_ff_blu = vc-v_ff
  v_ff_red = vc+v_ff
  pv_ff_blu = np.array([(halflength-R_as).to(u.deg),v_ff_blu.to(u.m/u.s)])
  pv_ff_red = np.array([(halflength+R_as).to(u.deg),v_ff_red.to(u.m/u.s)])

  #F4.show_lines([pv_ff_blu],color='blue',linewidths=1.0,linestyles='-')
  #F4.show_lines([pv_ff_red],color='red',linewidths=1.0,linestyles='-')

  #### PURE KEPLERIAN
  v_rot_blu = vc-v_rot
  v_rot_red = vc+v_rot
  pv_rot_blu = np.array([(halflength-R_as).to(u.deg),v_rot_blu.to(u.m/u.s)])
  pv_rot_red = np.array([(halflength+R_as).to(u.deg),v_rot_red.to(u.m/u.s)])

  F4.show_lines([pv_rot_blu],color='cyan',linewidths=2.0,linestyles='-')
  F4.show_lines([pv_rot_red],color='pink',linewidths=2.0,linestyles='-')

  #### KEPLERIAN + ANGULAR MOMENTUM
  #v_am_blu = vc-v_am
  #v_am_red = vc+v_am
  #v_am_blu[R_cam>R_pc] = vc-v_rot[R_cam>R_pc]
  #v_am_red[R_cam>R_pc] = vc+v_rot[R_cam>R_pc]
  #pv_am_blu = np.array([(R_as+halflength).to(u.deg),v_am_blu.to(u.m/u.s)])
  #pv_am_red = np.array([(halflength-R_as).to(u.deg),v_am_red.to(u.m/u.s)])

  #F4.show_lines([pv_am_blu],color='blue',linewidths=1.0,linestyles='--')
  #F4.show_lines([pv_am_red],color='red',linewidths=1.0,linestyles='--')

  #### KEPLERIAN + INFALL
  v_inf_blu = vc-v_inf
  v_inf_red = vc+v_inf
  pv_inf_blu = np.array([(halflength-R_as).to(u.deg),v_inf_blu.to(u.m/u.s)])
  pv_inf_red = np.array([(halflength+R_as).to(u.deg),v_inf_red.to(u.m/u.s)])

  F4.show_lines([pv_inf_blu],color='cyan',linewidths=2.0,linestyles='--')
  F4.show_lines([pv_inf_red],color='pink',linewidths=2.0,linestyles='--')

  pv_rot_blu = np.array([(halflength-R_as).to(u.deg),v_rot_blu.to(u.m/u.s)])
  pv_rot_red = np.array([(halflength+R_as).to(u.deg),v_rot_red.to(u.m/u.s)])

  F4.show_lines([pv_rot_blu],color='cyan',linewidths=2.0,linestyles='-')
  F4.show_lines([pv_rot_red],color='pink',linewidths=2.0,linestyles='-')

pv_rot_blu = np.array([rvalarcs.to(u.deg),velvalues*u.m/u.s+vc.to(u.m/u.s)])'''
######

F4.refresh()

##INFALL: (try with the different values of M obtained from the 2d-fit performed in CASA)

################################
#### PANEL (a): Moment 0 map
################################
## FIRST MAKE MOMENT MAP:
moment_0 = cube0.moment(order=0)
## NEXT PLOT (next to the PV plot):
f2 = aplpy.FITSFigure(moment_0.hdu,dimensions=[0,1],figure=fig3,subplot=(leftpanel))
standard_setup(f2) 
f2._ax1.yaxis.labelpad = -20 ## Move label closer to figure #MOD
mom0_vmin = np.nanmin(moment_0).value*0.4 #modify multiplicative factor to make #MOD better contrast
mom0_vmax = np.nanmax(moment_0).value*0.6 #Shorten colorscale to make it more sensible for a higher number of pixels. #MOD
f2.show_grayscale(vmin=mom0_vmin,vmax=mom0_vmax,invert=True)
xx=np.abs(cube0.header['NAXIS1']*cube0.header['cdelt1'])
yy=np.abs(cube0.header['NAXIS2']*cube0.header['cdelt2'])
mapwidth=5.*u.arcsec #MOD
mapheight=5.*u.arcsec #MOD
f2.recenter(RAcoord,Deccoord,width=mapwidth.to(u.deg).value,height=mapheight.to(u.deg).value)
f2.show_markers(RAcoord, Deccoord,marker="*",c='None',edgecolors="magenta",s=300,linewidths=2.0,zorder='9')#MOD
 
#### TO SHOW PV CUT
f2.show_lines([endpoints_deg.T],color='DarkMagenta',linewidths=2.0,linestyle='-')
#BEAMSIZE
bmaj=1.043*u.arcsec.to(u.deg) #MOD
bmin=0.7046*u.arcsec.to(u.deg) #MOD
bpa=76.299 #MOD
f2.show_beam(major=bmaj,minor=bmin,angle=bpa,fill=True,color='DeepSkyBlue') #MOD
f2.add_label(0.22, 0.9, '(a) C18O map', relative=True, color="Black", size=15) #MOD
 
################################
## FINAL ANNOTATIONS
################################
# HORIZONTAL LINE TO SHOW resolution in PV Plot. 
res_elem = 0.15*u.arcsec ##?## Where do I get this?? #resolution element, or beam size #MOD
x1_res=res_elem.to(u.deg).value 
x2_res=x1_res+res_elem.to(u.deg).value 
y1_res=pvcube.spectral_axis.min().value
y2_res=pvcube.spectral_axis.max().value 
 
#### BEAMSIZE #MOD
F4.show_arrows(30080.,3./3600.,0.0,0.5*res_elem.to(u.deg).value,color='DeepSkyBlue',linestyles='-',head_length=0.1,head_width=5) 
F4.show_arrows(30080.,3./3600.,0.0,-0.5*res_elem.to(u.deg).value,color='DeepSkyBlue',linestyles='-',head_length=0.1,head_width=5)
F4.add_label(23000.,3./3600.,'0.14$^{\prime\prime}$ \n 60 AU',layer='label',color='DeepSkyBlue',size=16, verticalalignment='center')
 
# Set up a second Axis on top of pv plot, to show annotations (if necessary) ##?## invisible axes just to make annotations
ax2 = fig3.add_axes(rightpanel) 
ax2.set_xlim(xlim_lo,xlim_hi)
ax2.set_ylim(ylim_lo,ylim_hi)
 
ax2.set_axis_off()
 
# VERTICAL (dashed) LINE TO SHOW v_c in PV Plot
hx1 = 0.
hx2 = 2.*halflength.to(u.deg).value
hy1 = vc.to(u.m/u.s).value
ax2.vlines(hy1,hx1,hx2,colors='green',linestyles='dashed',label='$v_c$',zorder='3',lw=2) ##?## ##JP## I am getting vertical line shifted #MOD
#ax2.annotate('$v_c$',xy=(hy1,x1_res*2.0),xycoords='data',color='green',size='16')
 
# HORIZONTAL (dashed) LINE TO SHOW star position in PV Plot. Separate position to left and right from the source along the position line. 
hx3 = xlim_lo
hx4 = xlim_hi
hy3 = halflength.to(u.deg).value
#### NOTE, You have to move these around a bit depending on your target. #MOD
ax2.hlines(hy3,hx3,hx4,colors='magenta',linestyles='dashed',label='protostar',zorder='3',lw=2)
#ax2.annotate('protostar',xy=(xlim_hi,hy3-2*res_elem.to(u.deg).value),xycoords='data',color='magenta',size='20',verticalalignment='top')
#ax2.annotate('(b) CO, \n PV diag.',xy=(xlim_hi,hy3+8*res_elem.to(u.deg).value),xycoords='data',color='Black',size='20',verticalalignment='top')
F4.add_label(0.18,0.9,'(b) C18O \n PV diag.',relative=True,color='Black',size=15)
F4.add_label(0.8,0.55,'protostar',relative=True,color='magenta',size=13)
F4.add_label(0.5,0.1,'Vc',relative=True,color='green',size=13)

##SHOW/SAVE FIGURE:
#plt.show()
if writechan: 
  plt.savefig(outname)
if not writechan:
  print 'print to screen'
  fig3.show() 
#plt.close(fig3)


