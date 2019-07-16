'''
Created on Jun 27, 2018

@author: Keisuke
'''

from pylab import *
from Labdata02 import *
from Savefig3 import Savefig3
from Savetxt import *

# rc('text', usetex=True)
# rc('font', family='serif')
# rc('font', serif='Times New Roman')

rcParams.update({'mathtext.default': 'regular' })

#Plot setting##############
rcParams.update({'font.size': 22})
rcParams.update({'lines.markersize': 10})
rcParams.update({'lines.markeredgewidth': 0.5})
rcParams.update({'font.size': 15})
rcParams.update({'lines.linewidth': 2.5})
rcParams.update({'figure.autolayout': True})
rcParams['figure.figsize']=5.5,4
rcParams.update({'figure.facecolor':'W'})
rcParams.update({'axes.linewidth':1.2}) 
rcParams.update({'xtick.major.width':1.2})  #
rcParams.update({'ytick.major.width':1.2})  #
rcParams.update({'lines.markeredgewidth': 1}) #default = 1
rcParams.update({'text.latex.preamble': ['\\usepackage[greek,english]{babel}']})
rcParams.update({'xtick.major.pad': 10})
rcParams.update({'ytick.major.pad': 8})
###########################

savefolder='Croco\\Light06'
dpi=600

Xlabel='$I$ ($\mu$mol m$^{-2}$ s$^{-1}$)'

def gft(name):
    return genfromtxt('Light06\\'+name+'.csv',delimiter=',')

def FigSetting():
    xticks([0,500,1000,1500])

Iarray = gft('IarrayN2fix')
N2fix = gft('N2fix')
Iarray2 = gft('IarrayPchl')
PI = gft('Pchl')

C2H4Nratio=1.5    #(molC2H4 molN-1) conversion
Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll (from a 704)

figure(11)
ax1 = gca(); ax2 = twinx(ax1)
ln1=ax2.plot([25,50,100,140,190,250,375,500,750,1000],array([0,0,0,4.17,3.5,3.36,4.85,2.86,2.58,2.37])*2,'o',color='cyan',label='Data')    #Data from Masuda
ln2=ax1.plot(Iarray,N2fix/C2H4Nratio,label='Model',color='b')
lns = ln1 + ln2
ax1.set_xlabel(Xlabel)
ax1.set_ylim(ymin=0,ymax=15)
ax2.set_ylim(ymin=0,ymax=15/1.8*1.65)
ax2.set_ylabel('N$_{2}$fix (fmolC$_{2}$H$_{4}$ cell$^{-1}$ d$^{-1}$)')
ax1.set_ylabel('N$_{2}$fix (fmolN cell$^{-1}$ d$^{-1}$)')
labs = [l.get_label() for l in lns]
legend(lns, labs, loc='lower right',fontsize=18)
FigSetting()
xlim(0,1000)
Savefig3(savefolder,11,dpi)

figure(12)
ax1 = gca(); ax2 = twinx(ax1)
ln2 = ax1.plot(Iarray2,PI*3600,'g',label='Model')
ln1 = ax2.plot([20,60,100,250,450,700,900,1200,1600],array([51,192,327,765,1134,1279,1203,847,578])/1000*Mchl,'o',color='#00FF00',label='Data')    #from "XXXXXE_ETR Crocosphaera.xlsx"

ax1.set_xlabel(Xlabel)
ax1.set_ylabel('$P_{Chl}$ (molC mol Chl$^{-1}$ h$^{-1}$)')
ax2.set_ylabel('$ETR$ (mol e$^-$ mol Chl$^{-1}$ h$^{-1}$)')
ax2.set_ylim(0,1280)
lns = ln1 + ln2
labs = [l.get_label() for l in lns]
legend(lns, labs, loc='lower right',fontsize=18)
FigSetting()
xlim(0,1700)
Savefig3(savefolder,12,dpi)

show()