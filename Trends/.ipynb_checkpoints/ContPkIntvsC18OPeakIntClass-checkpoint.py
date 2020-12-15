import matplotlib.pyplot as plt
import numpy as np

filenametxt='/home/jp/Research/Data/PeakFluxMassIntGasClass.txt'
contpkint,pkintC18O,sclass= np.loadtxt(filenametxt,unpack=True,usecols=(1,4,5),dtype='S')

lcontpkint=[]
lpkintC18O=[]
for i in range(len(pkintC18O)):
	lpkintC18O.append(float(pkintC18O[i]))
	lcontpkint.append(float(contpkint[i]))

lcontpkintN=[]
lcontpkint0_I=[]
lcontpkintF=[]
lpkintC18ON=[]
lpkintC18O0_I=[]
lpkintC18OF=[]
for i in range(len(pkintC18O)):
	if sclass[i]=='N':
		lpkintC18ON.append(float(pkintC18O[i]))
		lcontpkintN.append(float(contpkint[i]))
	elif sclass[i]=='F':
		lpkintC18OF.append(float(pkintC18O[i]))
		lcontpkintF.append(float(contpkint[i]))
	else:
		lpkintC18O0_I.append(float(pkintC18O[i]))
		lcontpkint0_I.append(float(contpkint[i]))

figname='Continuum Peak Intensity vs. C18O Peak Intensity'
fig=plt.figure('Continuum Peak Intensity vs. C18O Peak Intensity',figsize=(18,12))
plt.title('Continuum Mass vs. C18O Peak Intensity')
plt.xlabel('Intensity [mJy/bm]')
plt.ylabel('Intensity [Jy/bm]')
plt.rc('axes', titlesize=18)
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
plt.xticks(np.arange(min(lcontpkint),max(lcontpkint),10.))
plt.yticks(np.arange(min(lpkintC18O),max(lpkintC18O),0.05))
plt.scatter(lcontpkintF,lpkintC18OF,color='blue',label='F')
plt.scatter(lcontpkintN,lpkintC18ON,color='gray',label='None')
plt.scatter(lcontpkint0_I,lpkintC18O0_I,color='red',label='Class 0-I')
plt.legend()
fig.savefig(figname+'.png')
plt.show()
plt.close(fig)
