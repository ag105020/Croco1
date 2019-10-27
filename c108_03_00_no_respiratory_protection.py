'''
Created on Mar 13, 2019

@author: keiin
'''

from pylab import *
from Savefig3 import *

#======================================================
# Figure Setting
#======================================================
rcParams.update({'lines.markersize': 10})
rcParams.update({'lines.markeredgewidth': 0.5})
rcParams.update({'font.size': 15})
rcParams.update({'lines.linewidth': 2.5})
rcParams.update({'figure.autolayout': True})
rcParams['figure.figsize']=5,4
rcParams.update({'figure.facecolor':'W'})
rcParams.update({'axes.linewidth':1.2}) 
rcParams.update({'xtick.major.width':1.2})  #
rcParams.update({'ytick.major.width':1.2})  #
rcParams.update({'lines.markeredgewidth': 1}) #default = 1
rcParams.update({'mathtext.default': 'regular' })

#==========================================
# Pre-defining parameters for defs
#==========================================
n = "00"
n1 = "03"
dpi = 600
Tmax = 86400 #(s) maximum time step (177-32)
dT = 50    #(s) time step    (177-32)
Days = 1
T = arange(0.,Tmax*Days+10**(-10),dT)
Th = T/3600
Days = 2
T2 = arange(0.,Tmax*Days+10**(-10),dT)
Th2 = T2/3600
C2H4Nratio = 1.5
color05 = 'k'
color20 = 'r'
color051 = 'gray' 
color201 = 'pink'
Xlabel = '$t$ (h)'
xedge = 0.2

#==========================================
# Defining functions
#==========================================
def gft(name):
    return genfromtxt('C:\\Users\\keiin\\OneDrive\\Desktop\\figures\\01\\Croco\\Output01\\' + n + name+ '.csv', delimiter=',')

def g05(name):
    return genfromtxt('C:\\Users\\keiin\\OneDrive\\Desktop\\figures\\01\\Croco\\Output01\\' + n + name + '05.csv', delimiter=',')

def g20(name):
    return genfromtxt('C:\\Users\\keiin\\OneDrive\\Desktop\\figures\\01\\Croco\\Output01\\' + n + name + '20.csv', delimiter=',')

def g05b(name):
    return genfromtxt('C:\\Users\\keiin\\OneDrive\\Desktop\\figures\\01\\Croco\\Output01\\' + n1 + name + '05.csv', delimiter=',')

def g20b(name):
    return genfromtxt('C:\\Users\\keiin\\OneDrive\\Desktop\\figures\\01\\Croco\\Output01\\' + n1 + name + '20.csv', delimiter=',')

def sf00(name):
    Savefig3('Croco\\Output01\\Figs\\03NoDiffPro',name,dpi)

def dp05(filename):
    if filename=='N2fix05all':
        da=genfromtxt('data02\\'+filename+'.txt',)
        errorbar(da[:,0],da[:,1]/C2H4Nratio,(da[:,5]/C2H4Nratio,da[:,4]/C2H4Nratio),fmt='o',color=color05,elinewidth=1.5,ecolor=color05,capthick=1.5,capsize=5)

    else:
        da=genfromtxt('data02\\'+filename+'.txt',)
        errorbar(da[:,0],da[:,1],(da[:,5],da[:,4]),fmt='o',color=color05,elinewidth=1.5,ecolor=color05,capthick=1.5,capsize=5)
    
def dp20(filename):
    if filename=='N2fix20all':
        da=genfromtxt('data02\\'+filename+'.txt',)
        errorbar(da[:,0],da[:,1]/C2H4Nratio,(da[:,5]/C2H4Nratio,da[:,4]/C2H4Nratio),fmt='o',color=color20,elinewidth=1,ecolor=color20,capthick=1,capsize=5)
    else:
        da=genfromtxt('data02\\'+filename+'.txt',)
        errorbar(da[:,0],da[:,1],(da[:,5],da[:,4]),fmt='o',color=color20,elinewidth=1,ecolor=color20,capthick=1,capsize=5)      

def dp051(filename):
    da=genfromtxt('data02\\'+filename+'.txt',)
    errorbar(da[:5,0],da[:5,1],(da[:5,5],da[:5,4]),fmt='o',color=color051,elinewidth=1.5,ecolor=color051,capthick=1.5,capsize=5)

def dp201(filename):
    da=genfromtxt('data02\\'+filename+'.txt',)
    errorbar(da[:5,0],da[:5,1],(da[:5,5],da[:5,4]),fmt='o',color=color201,elinewidth=1,ecolor=color201,capthick=1,capsize=5)

def pl(fignun,name,m05,m20,m05b,m20b,d05,d20,Ylabel):
    figure(fignun)
    plot(Th,m05,color=color05)
    if name=='N2fix':
        plot(Th,m05b,dashes=(5,2),color=color05)
    else:
        plot(Th,m05b,'--',color=color05)
    dp05(d05)
    plot(Th,m20,color=color20)
    plot(Th,m20b,'--',color=color20)
    dp20(d20)
    xlabel(Xlabel)
    ylabel(Ylabel)
    xlim([0,24+xedge])
    xticks(arange(0,24+1,2))
    sf00(name)

#====================================
# Generating from files
#====================================
#----Default-------------
N2fix05 = g05('N2fix')
Pho05 = g05('Pho')
Res05 = g05('Res')
C05 = g05('C')
N05 = g05('N')
CN05 = g05('CN')
O2_05 = g05('_O2_')

N2fix20 = g20('N2fix')
Pho20 = g20('Pho')
Res20 = g20('Res')
C20 = g20('C')
N20 = g20('N')
CN20 = g20('CN')
O2_20 = g20('_O2_')

#----Altered simulation---------------
N2fix05b = g05b('N2fix')
Pho05b = g05b('Pho')
Res05b = g05b('Res')
C05b = g05b('C')
N05b = g05b('N')
CN05b = g05b('CN')
O2_05b = g05b('_O2_')

N2fix20b = g20b('N2fix')
Pho20b = g20b('Pho')
Res20b = g20b('Res')
C20b = g20b('C')
N20b = g20b('N')
CN20b = g20b('CN')
O2_20b = g20b('_O2_')

#==================================
# Plotting
#==================================
figure(3)
# plot(Th,Pho05,color=color051)
# plot(Th,Pho05b,'--',color=color051)
# dp051('O2fluxPho05all')
plot(Th,Res05,color=color05)
plot(Th,Res05b,'--',color=color05)
dp05('O2fluxRes05all')
# plot(Th,Pho20,color=color201)
# plot(Th,Pho20b,'--',color=color201)
# dp201('O2fluxPho20all')
plot(Th,Res20,color=color20)
plot(Th,Res20b,'--',color=color20)
dp20('O2fluxRes20all')
xlabel(Xlabel)
ylabel('O$_2$ flux (fmol O$_2$ cell$^{-1}$ h$^{-1}$)')
xlim([0,24+xedge])
xticks(arange(0,24+1,2))
ylim(bottom=-0.72,top=18);yticks(arange(0,19,2))
sf00('O2flux')

pl(1,'N2fix',N2fix05,N2fix20,N2fix05b,N2fix20b,'N2fix05all','N2fix20all','N$_{2}$fix (fmolN cell$^{-1}$ h$^{-1}$)')
ylim(bottom=-0.08); sf00('N2fix')
pl(4,'C',C05,C20,C05b,C20b,'Ccell05all','Ccell20all','C per cell (fmol C cell$^{-1}$)')
pl(5,'N',N05,N20,N05b,N20b,'Ncell05all','Ncell20all','N per cell (fmol N cell$^{-1}$)')
pl(6,'CN',CN05,CN20,CN05b,CN20b,'CNcell05all','CNcell20all','Molar C:N ratio (molC molN$^{-1}$)')

figure(7)
plot(Th,O2_05,color=color05)
plot(Th,O2_05b,'--',color=color05)
plot(Th,O2_20,color=color20)
plot(Th,O2_20b,'--',color=color20)
xlabel(Xlabel)
ylabel('O$_2$ cell (molO$_2$ m$^{-3}$)')
xlim([0,24+xedge])
xticks(arange(0,24+1,2))
ylim([-0.04,1])
yticks(arange(0,1.01,0.2)) 
sf00('O2')

show()


