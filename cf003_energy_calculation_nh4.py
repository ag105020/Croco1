'''
Created on Jan 24, 2016

@author: Keisuke
'''




    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Computation of fe0, fpr and fn considering material, redox and energy balance
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
from pylab import *

class Evalue_nh4:
    def __init__(self,E,O2_bioC_ratio):
        self.E=E
        self.O2_bioC_ratio=O2_bioC_ratio
        


def evalue_nh4(ep,CNratio):

    n=5         #number of carbon in biomass (167-6)
    a=7         #number of hydrogen in biomass (167-6)
    b=2         #number of oxygen in biomass (167-6)
    c=n/CNratio         #number of nitrogen in biomass (167-6)
#     c=1
#     d=1/30      #number of phosphorus in biomass (167-6)
#     d=0
    C=6         #number of carbon in sugra substrate (glucose->6)
    
 #   f=4*n+5*c+5*d-2*b+a             #(from 167-6: Nitrate case)
    d=4*n+a-2*b-3*c             #(from Rittman and Mccarty: NH4+ case)
    g=12*n+a+16*b+14*c     #(from 167-6)
    y=c/d           #coefficient of NH4+ in the half reaction of BB production
    z=1/4           #coefficient of NH4+ in the half reaction of nitrogen-fixation
    
#    fs00=(1/y)/(1/y+1/z)        #The ratio of electron used for protein synthesis
    fs00=1
#    fn00=(1/z)/(1/y+1/z)        #The ration of electron used for nitrogen fixation
    fn00=0

    dgc0=41.35      #The free energy necessary for the half reaction of glucose production (kJ/e-mol)
    dgATP=50      #(kJ/ATP): energy produced by the reaction of ATP -> ADP (147-19)
    dgn=2*dgATP-dgc0*ep   #The free energy necessary (dg) for the half reaction of nitrogen fixation (kJ/e-mol)
    dgp=35.09-dgc0  #dg for production of pyruvate from glucose (kJ/e-mol))
    dgpc=3.33*g/d  #dg for the production of BB (bacterial biomass) from pyruvate) (147-17)
    dgr=-120.07     #-dg for the energy production pathway (kJ/e-mol)
#    ep=0.6          #energy efficiency for the production of energy and the consumption of energy
    if dgp<0:       
        ep1=1/ep    #change ep1 depending on the sign of dgp    
    else:
        ep1=ep
    A=(fn00*dgn+fs00*(dgp/ep1+dgpc/ep))/(-ep*dgr) #A is related to fs0 and fe0
    fe0=A/(1+A)     #the ratio of electron used for energy production
    
    fs0=1-fe0       #the ratio of electron used for biomass synthesis+nitrogen fixation
    fpr=fs0*fs00    #the ratio of electron used for biomass synthesis
    fn=fs0*fn00     #the ratio of electron used for nitrogen fixation

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #2.--Stoichiometry (to get E1(E for the case O2cri>O2in))------------------
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    S1=array([["","CH","H2O","CO2","O2","HCO3-","NH4+","N2","H2","BB","H+","e-"],
              ["'-Rd",1/24,0.25,-0.25,0.,0.,0.,0.,0.,0.,-1.,-1.],
              ["Ra",0.,-0.5,0.,0.25,0.,0.,0.,0.,0.,1.,1.],
              ["Rpr",0.,-(2*n-b+c)/d,(n-c)/d,0.,c/d,c/d,0.,0.,-1/d,1.,1.],
              ["Rn",0.,0.,0.,0.,0.,-0.25,0.125,-0.125,0.,1.25,1.]])

    #=====================================================================================
    #NH4+ absorption case
    #=====================================================================================
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S2 (f*R for electron acceptance)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    width=12
    depth=5    
    Mua=array([[1],[fe0],[fpr],[0.0]])  #column of f
    S2a=copy(S1[1:,1:])               #use copy so S2 does not respond to the change in S1
    S2a=S2a.astype(float64)
    #add numbers for columns and raws for counting 
    S21a=arange(1,width,1)
    S22a=arange(0,depth,1)
    S22a=S22a.reshape(depth,1)
    S2a=vstack((S21a,S2a))
    S2a=hstack((S22a,S2a))
    #calculate f*R
    S2a[1:,1:]=Mua*S2a[1:,1:]

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S3 (f*R for electron donation)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    S3a=vstack((S2a[1,1:],S2a[1,1:],S2a[1,1:]))
    S3a=Mua[1:]*S3a
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S4 (f*R for "electron donation + electron acceptance")
    # and RR, which is the entire reaction
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    S4a=S2a[2:,1:]+S3a     #S4 is f*R for "electron donation + electron acceptance"
    RRa=S4a[0]+S4a[1]+S4a[2]
    RR1a=copy(RRa)
    RR1a=vstack((S1[0,1:],RR1a))    
    pa=-S4a[1,8]*n                        #(molC/e-mol) CH consumption for protein production
    ha=(S4a[1,0]+S4a[2,0])*C               #(molC/e-mol) CH consumption for other than energy production
    alpa=pa/ha                             #(CO2 from (Ra-Rd))/(CH for protein production) (see 72-5)
    betaa=-(S4a[1,2]+S4a[1,4]+S4a[2,2])/ha   #(CO2 from (Rpr-Rd) + CO2 from (RN-Rd))/(CH for protein production) (see 72-5)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #output of each array into CSV files
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

 
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Getting yield (Y) and the ratio of CO2 production rate to CH consumption (E)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Y1=-(RRa[8]*n)/(C*RRa[0])              #Yield
    E3=1/Y1-1 
    E=E3
    O2_bioC_ratio=-RRa[3]/(RRa[8]*n)
    
    E1=Evalue_nh4(E,O2_bioC_ratio)
    
    return(E1)
    