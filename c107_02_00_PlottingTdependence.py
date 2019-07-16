'''
Created on Mar 13, 2019

@author: keiin
'''


from pylab import *
from Savefig3 import *

rcParams.update({'font.size': 25,
                 'lines.markersize':12,
                 'lines.markeredgewidth':1})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'font.serif': 'Times New Roman'})
rcParams.update({'figure.autolayout': True})
rcParams['figure.figsize']=8,6.5
rcParams.update({'figure.facecolor':'W'})
rcParams.update({'lines.linewidth':3})   
rcParams.update({'patch.edgecolor':'none'})

rcParams.update({'axes.linewidth':1.5})
rcParams.update({'xtick.major.width':1})
rcParams.update({'ytick.major.width':1})
rcParams.update({'mathtext.default': 'regular' })

def gft(name):
    return genfromtxt('C:\\Users\\keiin\\OneDrive\\Desktop\\figures\\01\\Croco\\Tdepend02\\' + name, delimiter=',')

def sf(name):
    Savefig3('Croco\\Tdepend02\\Paper2',name,600)

C2H4Nratio = 1.5

T = gft('Tall\\IarrayN2fix04.csv')
Photo = gft('Tall\\Photo.csv')
Resp = gft('Tall\\Resp.csv')
N2fix0 = gft('Tall\\N2fix.csv')
N2fix1 = gft('TnfixOnly\\N2fix.csv')
N2fix2 = gft('TphotoONly\\N2fix.csv')
N2fix3 = gft('TrespOnly\\N2fix.csv')

Xlabel='$T$ ($^\circ$C)'

figure(1,figsize=(9,6.5))
ax1=gca(); ax2=ax1.twinx()
ln1=ax1.plot(T,N2fix0/C2H4Nratio,label='N$_2$fix')
ln2=ax2.plot(T,Resp,color='red',label='Resp')
ln3=ax2.plot(T,Photo,color='green',label='Photo')
ax1.set_xlabel(Xlabel)
lns=ln1+ln2+ln3
labs = [l.get_label() for l in lns]
legend(lns, labs, loc='upper left',fontsize=25)
title('$T$ dependence',y=1.03)
ax1.set_ylabel('N$_{2}$fix (fmolN cell$^{-1}$ d$^{-1}$)')
ax2.set_ylabel('O$_2$ fluxes (fmol cell$^{-1}$ h$^{-1}$)')
ax1.set_ylim(0,10)
ax2.set_ylim(0,160)
xlim(14,28)
sf('Tdepend')

figure(2)
plot(T,N2fix0/C2H4Nratio,'b',label='All')
plot(T,N2fix1/C2H4Nratio, color='b',dashes=[3,5],label='N$_{2}$fix')
plot(T,N2fix2/C2H4Nratio, color='b',dashes=[11,7],label='Photo')
plot(T,N2fix3/C2H4Nratio, color='b',dashes=[3,3,10,3],label='Resp')
xlim(14,28)
ylim(0,10)
title('T effects on N$_2$fix',y=1.03)
ylabel('N$_{2}$fix (fmolN cell$^{-1}$ d$^{-1}$)')
xlabel(Xlabel)
legend(loc='upper left',fontsize=25)
sf('N2fixTdepend')

show()


