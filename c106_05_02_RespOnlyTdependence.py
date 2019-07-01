
'''
Here I run the model for various temperatuer with O2 variation
'''

from pylab import *
from Labdata02 import *
from cf002_energy_calculation_n2 import *
from cf003_energy_calculation_nh4 import *
from Savefig3 import Savefig3
from Savetxt import *
import time
from TemperatureViscosity import *   #in the source folder "Functions"
from DissolvedO2Saturation import *

t0=time.time()

rcParams.update({'mathtext.default': 'regular' })

###############################################################################
def croco(data,Temperature):
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Getting values from data
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Datanumber=data.Tag
    Da=copy(Datanumber)
    I_const=data.Irradiance
    Days=data.Days
    light_period_length=data.Daytime  #(s)length of light period
    dark_period_length=data.Nighttime    #(s)length of dark period
    #Temperature=data.Temperature    #(K) temperature of the environment
    Kelvin = 273.15
    Temperature = Temperature + Kelvin
    O2_v=data.O2               #(molO2 m-3) environmental oxygen concentration (mol m-3) (post 1982)
    O2_v = DissolvedO2Saturation(Temperature-Kelvin)/1000
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Choosing the data to work on 
    #0:Mohr 2010, 1:GroBkopf 2012 5% O2, 
    #2:GroBkopf 2012 20% O2, 3:Dron 2012, 4: Satio 2011
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    if (Da==0) or (Da==1) or (Da==3) or (Da==4):                    #Plotting GroBkopf 2012 O2 20%
        return
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Time setting
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Tmax=86400 #(s) maximum time step (177-32)
    dT=100
    T=arange(0.,Tmax+10**(-10),dT)  #(s) time array (177-32)
    
    #------------------------------
    #Size related values
    #------------------------------
    RLg0=1.32e-6
    RLg=RLg0    #(m) radius of the cell including the cell membrane
    if Datanumber==1:
        RLg=1.28e-6
    
    Rg=1/27         #(dimensionless) ratio of glycolipid layer to the cell radius
    R=RLg/(1+Rg)    #(m) radius of the cell without including the cell membrane
    V=4/3*pi*RLg**3   #(m3) volume of the cell
    Lg=R*Rg         #(m3) thickness of the cell membrane (g=glycolipid)
    RLgoriginal=1.5e-6    #(m) originally set radius of the cell including the cell membrane
    Voriginal=4/3*pi*RLgoriginal**3   #(m3) originally set volume of the cell 
    Sizeadjustment=V/Voriginal   #(dimensionless) size adjustment term
    Sizeadjustment2=RLg0**3/RLg**3  #Counteracting size effect
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Values to play with
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #*******************************************
    #play
    #*******************************************
    RFeChl=5.236e-3     #(molFe molC chl-1) Fe in the photosystem to chol carbon ratio
    N2fixperFe=0.1            #(molN s-1 molFe-1) N2 fixation rate per nitrogenase iron (186-28) (value determination 186-36)
    Rpgly=6.32e-5          #(s-1) glycogen production rate (186-43)
    glymax_v=12383         #(molC m-3) maximum glycogen storage (per volume) (177-26)
    Qc_v=15000     #(molC m-3) biomass carbon content (storage not included) (per volume) (177-39)
    Qnc=1/6             #(molN molC-1) N to C molar ratio in the biomass except for storage (177-26)
    lsmaxmax_v=0.2607    #(molC m-3 s-1) maximum biomass productioin rate (per volume)
    Kcy_v=91.665            #(Nmol m-3) Half saturation of cyanophysin for biomass production (per volume)
    Kcy=Kcy_v*V          #(Nmol cell-1) Half saturation of cyanophysin for biomass production (per cell)
    cymax_v=47665.8   #(molN m-3) maximum storage of cyanophysin (per volume)(177-50)
    if Datanumber==1:
        cymax_v=cymax_v*Sizeadjustment2
    
    NFe=4869633         #(Fe cell-1) Number of Fe molecules in the cell
    Rchldecomposition=5e-4   #(dimensionless) Chlorophyll adjusting speed in degradation of chlorophyll (Rspeed2 in (186-35)
    epd=1.55e-5
    C2H4Nratio=1.5    #(molC2H4 molN-1) conversion
    erecycle=4/3        #(dimensionless) how much more N2 fixaction can happn by recycling the electron (100% recycle leads to 4/3 times more N2 fixation)
                        #read GroBkopf 2012 for more details
    N2fixidealpowerfactor=1    #(dimensionless) powerfactor on (cymax-cy)/cymax
    startingpoint=5*3600         #(s) counting starting point
    #**********************************************
    
    glymax=glymax_v*V   #(molC cell-1) maximum glycogen storage (per cell) (177-26)
    Qc=Qc_v*V           #(molC cell-1) biomass carbon content (storage not included) (per cell)
    lsmaxmax=lsmaxmax_v*V    #(molC cell-1 s-1) maximum biomass production rate (per cell)
    cymax=cymax_v*V     #(molN cell-1) maximum storage of cyanophycin (per cell)(177-50)
    Chlmin=7.13e-16*Sizeadjustment      #(molCchl cell-1) Chlorophyll min per cell (186-40)
    N2fixmax=7.5e-19*Sizeadjustment*Sizeadjustment2          #(molN s-1 cell-1) maximum nitrogen fixation rate when cyanophycin is depleated
    Rchlproductionmax=4.63e-5*Sizeadjustment  #(dimensionless)
    Kgly=9.25e-15*Sizeadjustment #(glyC cell-1) half saturation constant for glycogen consumption (177-28)
    #+++++++++++++++++++++++++++++++++++++
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Light dark period setting
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    light_period_step=light_period_length/dT    #light period in step
    light_or_dark=copy(T)*0     #array for light or dark 1=light, 0=dark
    light_or_dark[0:(light_period_step)]=1    #
    Oneday_light_or_dark=copy(light_or_dark)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Adjustment for multiple days
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    T=arange(0.,Tmax*Days+1e-10,dT)
    Th=T/3600
    U=T/dT
    Daysforloop=arange(0,Days-1,1)
    for i in Daysforloop:
        light_or_dark=hstack((light_or_dark,Oneday_light_or_dark))
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #light intensity setting
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    I=copy(light_or_dark)*I_const     #(uE m-2 s-1) light intensity (unit from Healey)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Parameter setting
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #==============================
    #constant values
    #==============================
    #-------------------------------
    #Diffusion related values
    #-------------------------------
    x0=1/epd*(1/R-1/(R+Lg))+1/(R+Lg)  #effect of glycolipid layer (Kei 103-16)
    r5=1/(x0*R)         #diffusivity efficiency of cell membrane
    Do2_25=2.12e-9 #Diffusion coefficient of O2 in the water at 25C (m2/s)

    #======================================================
    #Temperature adjustment of the diffusion coefficient
    #======================================================

    mu25 = 890.3
    Kelvin = 273.15
    mu = TemperatureViscosity(Temperature)
    Tfactor = Temperature/(25+Kelvin)/mu*mu25  #Lerned from 633 00 07
    #Conversion term is based on "15 Computation of diffusion coefficient for different temperature.xlsx"
    Do20 = Do2_25*Tfactor                  
    Do2 = Do20*r5         #Diffusion coefficient of O2 in the water (m2/s)
    #======================================================
    
    #--------------------------------
    #physiology related parameters
    #--------------------------------
    O2CH=1              #(O2mol CHmol) O2 production to Carbohydrate production ratio in photosynthtesis 
    
    #'''''''''''''''''''''''''''''''''''
    #night time specific parameters
    #'''''''''''''''''''''''''''''''''''
    Dglymaxmax=lsmaxmax*2.0745       #(molC m-3 s-1) maximum glycogen decomposition rate (tunable) (177-28)
    O2CH=O2CH           #(molO2 / molC CH) oxygen consumption to CH consumption ratio (177-50)(value: 72-04)
    
    #Temperature factor
    #Based on Geider 1997 and Deutsch 2015
    E0kB = 0.8e4 
    TempModifyFactor = exp(-E0kB*(1/Temperature-1/(28+Kelvin)))
    
    #-----------------------------------
    #photosynthesis related parameters
    #-----------------------------------
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll (from a 704)
    Pmax0=7                     #(g C /(g Chl h) Maximum production rate per chlorophyll (around 6 by Cullen 1990) (from a704-16 (in Healey 1985 model))
    Pmax=Pmax0*Mchl/12/3600/55       #(mol C s-1 mol chl-1) carbon fixing rate (156-10) (156-15) for unit conversion) (from a704-16)
    Pmax=0.004738
    O=0.01  #(m2 / umol photon) absorption cross section (from a704_16)
    Tau=0.5    #(s) handling time (from a704_16) (Irina: half life of RNA is less than 1 minute (30 (s)) or so)
    PI=Pmax*(1-exp(-O*I*Tau))#*TempModifyFactor     #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)
    
    
    
    #-----------------------------------------------------------
    #Iron introduction related parameter (mostly from Kei 186)
    #-----------------------------------------------------------
    Avogadro=6.02214e23    #(atoms mol-1) Avogadro constant  
    KFepho=4e6/Avogadro/10*Sizeadjustment    #(molFe cell-1) Half satulation of Fe in photosynthetic aparatus its degradation (186-36) based on Saito 2011.
    KFepho_v=KFepho/V       #(molFe m-3) Half satulation of Fe in photosynthetic aparatus in its degradation
    KFebuffer=copy(KFepho)  #(molFe cell-1) Half satulation of Fe in the buffer in its degradation (186-36)
    KFebuffer_v=KFebuffer/V   #(molFe m-3) Half satulation of Fe in the buffer in its degradation (186-36)
    KFenitrogenase=copy(KFepho)  #(molFe cell-1) Half satulation of Fe in nitrogenase in its degradation (186-36) 
    KFenitrogenase_v=KFenitrogenase/V #(molFe m-3) Half satulation of Fe in nitrogenase in its degradation (186-36)
    #++++++++++++++++++++++++++++++
    #importatnt value
    #++++++++++++++++++++++++++++++

    Rnitrogenasedegradation=1/300 #(dimensionless) Rate constant for nitrogenase degradation (Rspeed4 in 186-36)

    #---------------------------------------
    #Array preparation for each prarameter
    #---------------------------------------
    o=copy(T)*0             # this creates zero array for the right size for the time steps
    Ycnbio=copy(o)          #(mol C molN-1) real C:N ratio of biomass
    E=copy(o)               #(dimensionless) ratio of CO2 production to biomass production
    chl_v=copy(o)         #(molCchol m-3) chrolophyll concentration (per volume)
    cy_v=copy(o)        #(Nmol m-3) cyanophycin concentration (per volume)
    gly_v=copy(o)       #(Cmol m-3) glycogen concentration    (per volume)
    O2c_v=copy(o)         #(O2mol m-3) oxygen concentration in the cell (per volume)
    chl=copy(o)   #(molCchl cell-1) chlorophyll concentration (per cell)
    cy=copy(o)    #(molN cell-1) cyanophycin concentration (per cell)
    gly=copy(o)   #(molC cell-1) glycogen concentraiton (per cell)
    lsmax=copy(o) #(molC cell-1 s-1) biomass production rate given sufficient chlorophyll (per cell)
    chl_ratio=(o)   #(dimensionless) chl/chlideal    (177-26)
    chlideal=copy(o)   #(molCchl cell-1) ideal chlorophyll concentration (per cell) (177-27) 
    ls=copy(o)      #(molC cell-1 s-1) actual biomass production rate (per cell)  (177-26)
    Dcy=copy(o)     #(molN cell-1 s-1) cyanophycin decomposition rate (per cell) (177-26)
    Pglymax=copy(o) #(molC cell-1 s-1) glycogen production rate given sufficient chlorophyll (per cell) (177-26)
    Pgly=copy(o)    #(molC cell-1 s-1) actual glycogen production rate(per cell) (177-26)
    Ex=copy(o)      #(molC cell-1 s-1) excretion of carbohydrate (per cell) (177-26)
    Po2=copy(o)     #(molO2 cell-1 s-1) oxygen evolution rate (per cell) (177-42)
    Vo2=copy(o)     #(molO2 cell-1 s-1) oxygen uptake rate (per cell) (177-42)
    Dglymax=copy(o)   #(molC cell-1 s-1) 
    Pcy=copy(o)      #(molN cell-1 s-1) cyanophycin production rate (per cell) (177-28)
    Dgly=copy(o)    #(molC cell-1 s-1) decomposition of glycogen (per cell) (177-51)
    res=copy(o)     #(molO2 cell s-1)  Oxygen consumption rate (per cell) (177-51)
    N2fix=copy(o)   #(molN cell s-1) Nitrogen fixation rate (per cell) (177-51)
    N2fixideal=copy(o) #(for calculation (177-59)
    daytime=copy(o)     #(186-42)
    nighttime=copy(o)   #(186-42)
    Rnitrogenaseproduction=copy(o) #(dimensionless) Rate constant for nitrogenase production from ferritin (Rspeed3 in 186-36)

    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    #Iron introduction related parameter (mostly from Kei 186)
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    Fepho=copy(o)   #(molFe cell-1) Fe concentration in photosynthetic aparatus (186-35)
    Febuffer=copy(o)  #(molFe cell-1) Fe concentration in the buffer (186-35)
    Fenitrogenase=copy(o) #molFe cell-1) Fe concentration in nitrogenase (186-29)
    Fepho_v=copy(o) #(molFe m-3) Fe concentration in phototynthetic aparatus (186-35)
    Febuffer_v=copy(o) #(molFe m-3) Fe concentration in the buffer (186-35)
    Fenitrogenase_v=copy(o)  #(molFe m-3) Fe concentration in nitrogenase (186-29)
    Fenitrogenaseideal=copy(o)  #(molN s-1 cell-1) ideal nitrogenase amount per cell (186-29)
    Rchlproduction=copy(o)   #(dimensionless) Chlorophyll adjusting speed in production of chlorophyll (R=rapidity) (Rspped1 in (186-35)(177-27)
    resmax=copy(o)          #(molO2 s-1 cell-1) maximum respiration rate during the night time
    Vo2max=copy(o)          #(molO2 s-1 cell-1) maximum O2 uptake rate during the night time
    N2fixpotential=copy(o)      #(molN s-1 cell-1) potential nitrogen fixation rate based on the availability of nitrogenase
    Dglyideal=copy(o)       #(molC s-1 cell-1) the sum of N2 fixpotential and maximum oxygen uptake
    N2fixideal0=copy(o)     #(molN s-1 cell-1) ideal nitrogen fixation rate based on the availability of cyanophycin (187-4)
    Dglyideal0=copy(o)      #(molC s-1 cell-1) ideal glycogen decompositioin rate based on N2fixideal0 (187-4)
    which=copy(o)       #(dimensionless) To check which rout it takes
    which2=copy(o)  #(dimensionless) To check which rout it takes
    which3=copy(o)    #(dimensionless) To check which rout it takes
    nighttimeafter=copy(o)      #(s) time after certain time in the night time
    N2fixidealgly=copy(o)      #(molN s-1 cell-1) Ideal N2fixation rate in case of glycogen limiting
    Resmin=copy(o)          #(molO2 s-1 cell-1) minimum respiration rate for nitrogen fixation
    Resoriginal=copy(o)         #(mol O2 s-1 cell-1) original respiration before modification by res min (187-10)
    #==================================
    #Initial conditions
    #==================================
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Key initial conditions
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    cy_v[0]=733.32     #(molN m-3) original cyanophysin concentration (per volume)
    if Datanumber==0:
        cy_v[0]=229.26
    if Datanumber==1:
        cy_v[0]=229.26
    if Datanumber==2:
        cy_v[0]=229.26    #(molN m-3) original cyanophysin concentration (per volume)
    if Datanumber==3:
        cy_v[0]=229.26     #(molN m-3) original cyanophysin concentration (per volume)
    chl_v[0]=Chlmin/V    #(molCchl m-3) original chlorophyll concentration (per volume)

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #----------------------------
    #each parameter per volume
    #----------------------------
    gly_v[0]=0    #(molC m-3) original glycogen concentration    (per volume)
    if Datanumber==4:
        gly_v[0]=2200
    
    O2c_v[0]=O2_v   #(molO2 m-3) original oxygen concentration (per volume)

    #''''''''''''''''''''''''''''''''''''''''''''
    #Fe related parameters (mostly from kei 186)
    #''''''''''''''''''''''''''''''''''''''''''''
    Fepho_v[0]=chl_v[0]*RFeChl   #(molFe cell-1) original Fe concentration in photosynthetic aparatus (186-25)
    Febuffer_v[0]=(NFe-(Fepho_v[0]*V*Avogadro))/Avogadro/V   #(molFe cell-1) original Fe concentration in the buffer (186-25)
    #--------------------------------------------
    #each parameter per cell (_c means per cell)
    #--------------------------------------------
    chl[0]=chl_v[0]*V       #(molCchl cell-1) original chlorophyll concentration (per cell)
    cy[0]=cy_v[0]*V         #(molN cell-1) original cyanophysin concentration (per cell)
    gly[0]=gly_v[0]*V       #(molC cell-1) original glycogen concentration    (per cell)
    
    #''''''''''''''''''''''''''''''''''''''''''''
    #Fe related parameters (mostly from kei 186)
    #''''''''''''''''''''''''''''''''''''''''''''
    Fepho[0]=Fepho_v[0]*V  #(molFe cell-1) original Fe concentration in photosynthetic aparatus (per cell) (186-25)
    Febuffer[0]=Febuffer_v[0]*V   #(molFe cell-1) original Fe concentration in the buffer (per cell) (186-25)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for sentence for time steps
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    for i in U:
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Light period equations
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #======================================================================
    #Energetics and stoichiometry based on different C:N ratio in biomass
    #======================================================================
        Ycnbio[i]=(Qc+gly[i])/(Qc*Qnc+cy[i])
        ep=0.41          #energy efficiency for the production of energy and the consumption of energy
        RM_n2=evalue_n2(ep,Ycnbio[i])                 #RM means Rittmann and McCarty              
        RM_nh4=evalue_nh4(ep,Ycnbio[i])
        E_n2=RM_n2.E                #(dimensionless) the ratio of CO2 production to biomass production when the nitrogen source is N2 (nitrogen fixation)
        E_nh4=RM_nh4.E              #(dimensionless) the ratio of CO2 production to biomass production when the nitrogen sources is NH4+
        O2_bio_ratio_nh4_balanced=RM_nh4.O2_bioC_ratio
        
        #Temperature factors
        E0kB = 0.8e4 
        TempModifyFactor = exp(-E0kB*(1/Temperature-1/(28+Kelvin)))
        
    #=====================================================================   
     
        if light_or_dark[i]==1:
            if i>0:
                if light_or_dark[i-1]==0:
                    daytime[i]=0
                    
                else:
                    daytime[i]=daytime[i-1]+dT
            #-------------------------------------------------------------
            #obtaining current values    (overall solution map at 177-38)
            #-------------------------------------------------------------
            Rchlproduction[i]=(daytime[i]+dT/2)*1.1574e-9    #(186-42) learned from Rnitrogenaseproduction
            Rchlproduction[i]=Rchlproductionmax
            E[i]=E_nh4
            lsmax[i]=lsmaxmax*cy[i]/(cy[i]+Kcy/Ycnbio[i]/Qnc)#*TempModifyFactor     #(177-26)(167-1)(167-2)

            Pglymax[i]=(glymax-gly[i])*Rpgly#*TempModifyFactor   #(177-26)
            chlideal[i]=(lsmax[i]*(1+E[i])+Pglymax[i])/PI[i]   #(186-35)
            chl_ratio[i]=chl[i]/chlideal[i]        #(177-26)
            if chl_ratio[i]>1:
                chl_ratio[i]=1
            ls[i]=lsmax[i]*chl_ratio[i]       #(177-37)
            res[i]=ls[i]*O2_bio_ratio_nh4_balanced  #(177-26)
            Dcy[i]=ls[i]/Ycnbio[i]                        #(177-26)
            Pgly[i]=Pglymax[i]*chl_ratio[i]           #(177-26)
            Ex[i]=chl[i]*PI[i]-ls[i]*(1+E[i])-Pgly[i]             #(177-26)
            Po2[i]=PI[i]*chl[i]*O2CH                #(177-53)
            Vo2[i]=-Po2[i]+res[i]           #(177-42)
            O2c_v[i]=((-1)*Vo2[i])/(4*pi*R*Do2)+O2_v   #(177-54)
            Fenitrogenaseideal[i]=0         #(186-42)
            
            #---------------------------------------------------------------------
            #obtaining next time step values    (overall solution map at 177-38)
            #---------------------------------------------------------------------
            if i<(size(U)-1):
                cy[i+1]=cy[i]-Dcy[i]*dT                 #(177-36)
                gly[i+1]=gly[i]+Pgly[i]*dT              #(177-36)
            #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            #Fe related part (based mostly on Kei 186)
            #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                if chlideal[i]>chl[i]:
                    chl[i+1]=chl[i]+(chlideal[i]-chl[i])*Rchlproduction[i]*dT*(Febuffer_v[i]/(Febuffer_v[i]+KFebuffer_v))    #(186-35)(177-36)
                else:
                    if chlideal[i]<Chlmin:
                        chlideal[i]=Chlmin
                    chl[i+1]=chl[i]+(chlideal[i]-chl[i])*Rchldecomposition*dT*(Fepho_v[i]/(Fepho_v[i]+KFepho_v))     #(186-35)
                Fepho[i+1]=chl[i+1]*RFeChl       #(186-35)
                Fepho_v[i+1]=Fepho[i+1]/V   #(186-35)
                Fenitrogenase[i+1]=Fenitrogenase[i]+dT*(Fenitrogenaseideal[i]-Fenitrogenase[i])*Rnitrogenasedegradation*(Fenitrogenase_v[i]/(Fenitrogenase_v[i]+KFenitrogenase_v))     #(186-42)(186-39)
                Fenitrogenase_v[i+1]=Fenitrogenase[i+1]/V  
                Febuffer[i+1]=Febuffer[i]-(Fepho[i+1]-Fepho[i])-(Fenitrogenase[i+1]-Fenitrogenase[i])   #(186-42)(186-39)(186-35)
                Febuffer_v[i+1]=Febuffer[i+1]/V   #(186-35)
                      
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Dark period equations
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOO        
        elif light_or_dark[i]==0:
            if i>0:
                if light_or_dark[i-1]==1:
                    nighttime[i]=0
                else:
                    nighttime[i]=nighttime[i-1]+dT
            nighttimeafter[i]=nighttime[i]-(dark_period_length-startingpoint)
            if nighttimeafter[i]<0:
                nighttimeafter[i]=0
            
            N2fixmodification=(startingpoint-nighttimeafter[i])/startingpoint
            Vo2maxmodification=(startingpoint-nighttimeafter[i]/2)/startingpoint
            
            #=====================
            #something to play with
            #=====================
            Rnitrogenaseproduction[i]=6.7e-16*nighttime[i]**3*(-3.57*(O2c_v[i]-0.046)+1)*(gly[i]/(gly[i]+Kgly*5))**3     #(186-42)(187-7)
            #=====================
            
            resmax[i]=res[i-1]+nighttime[i]*3.125e-24        #(186-49)
            resmax[i]=nighttime[i]**2*2e-26
            
            resmaxMax = 1e-14/3600 #(mol O2 cell-1 s-1)
            
            if resmax[i] > resmaxMax:
                resmax[i] = resmaxMax
            
            E0kB = 0.8e4 
            TempModifyFactor = exp(-E0kB*(1/Temperature-1/(28+Kelvin)))
            
            resmax[i] = resmax[i]*TempModifyFactor
            E[i]=E_n2               #(177-50)
            Dglymax[i]=Dglymaxmax*gly[i]/(gly[i]+Kgly)     #(177-28)
            Vo2max[i]=4*pi*R*Do2*O2_v       #(177-58)(186-29)
            
            #==========================
            #Vo2max modification
            #==========================
            Vo2max[i]=Vo2max[i]*Vo2maxmodification
            #==========================
            if Vo2max[i]>resmax[i]:     #(186-49)
                Vo2[i]=resmax[i]
                O2c_v[i]=O2_v-resmax[i]/(4*pi*R*Do2)
                which3[i]=1
            else:
                Vo2[i]=Vo2max[i]
                O2c_v[i]=O2_v-Vo2[i]/(4*pi*R*Do2)
                which3[i]=2
            #=========================
            #revising based on 187-4
            #=========================
            Ocri=0.1
            if O2c_v[i]>Ocri: 
                N2fixpotential[i]=0
                which2[i]=1
            else:
                N2fixpotential[i]=Fenitrogenase[i]*N2fixperFe*(Ocri-O2c_v[i])/Ocri#*TempModifyFactor    #(187-4)(186-29)
                which2[i]=2
            Dglyideal[i]=N2fixpotential[i]/erecycle+Vo2[i]        #(187-4)
            if Dglymax[i]<Vo2[i]:               #(186-30)
                N2fixideal[i]=0                 #(186-30)
                N2fix[i]=0
                Vo2[i]=Dglymax[i]               #(186-30)
                Dgly[i]=Dglymax[i]
                
            else:
                if Dglyideal[i]>Dglymax[i]:
                    Dgly[i]=Dglymax[i]-(Vo2max[i]-Vo2[i])         #(check)
                    N2fix[i]=(Dglymax[i]-Vo2max[i])*erecycle
                    N2fixidealgly[i]=(Dglymax[i]-Vo2max[i])*erecycle          #maybe we should calculate N2fixideal in the next step part
                    which[i]=1
                    
                else:
                    Dgly[i]=N2fixpotential[i]/erecycle+Vo2[i]
                    N2fix[i]=N2fixpotential[i]
                    N2fixideal0[i]=((cymax-cy[i])/cymax)**N2fixidealpowerfactor*N2fixmax*gly[i]/(gly[i]+Kgly)#*(Ocri-O2c_v[i])/Ocri
                    Dglyideal0[i]=N2fixideal0[i]/erecycle+Vo2max[i]
                    if Dglyideal0[i]>Dglymax[i]:
                        N2fixideal[i]=(Dglymax[i]-Vo2max[i])*erecycle
                        which[i]=2
                    else:
                        N2fixideal[i]=N2fixideal0[i]
                        which[i]=3
            
            #==================================================
            #Modifying nitrogen fixation by the remaining time
            #==================================================
            N2fixideal[i]=N2fixideal[i]*N2fixmodification

            #==================================================            
            Pcy[i]=N2fix[i]   
            Fenitrogenaseideal[i]=N2fixideal[i]/N2fixperFe
            res[i]=Vo2[i]
            O2c_v[i]=O2_v-res[i]/(4*pi*R*Do2)

            #================================================
            #Modifying respiration with minimum respiration
            #================================================
            #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            #Nitrogen fixation energy balance
            #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            dgc0=41.35      #The free energy necessary for the half reaction of glucose production (kJ/e-mol)
            dgATP=50      #(kJ/ATP): energy produced by the reaction of ATP -> ADP (147-19)
            dgn=2*dgATP-dgc0*ep   #The free energy necessary (dg) for the half reaction of nitrogen fixation (kJ/e-mol)
            dgr=-120.07     #-dg for the energy production pathway (kJ/e-mol)
            An2fix=dgn/(-ep*dgr)    #(dimensionless) per-electron ratio of respiration that supports nitrogen fixtation to nitrogen fiation
            
            #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            Resmin[i]=N2fix[i]*8/4.28   #(187-10)    #value is based on GroBkopf (2012) reffering to Rarven (2009)
            Resmin[i]=N2fix[i]*An2fix   #from energetic balance
            if res[i]<Resmin[i]:
                Resoriginal[i]=copy(res[i])
                res[i]=Resmin[i]
                Dgly[i]=Dgly[i]+res[i]-Resoriginal[i]
            
            if Dgly[i]>Dglymax[i]:
                print("issue with respiration")
            #here I should update other values especially glucose consumption and think about whether resmin is possible
            #fix it when this becomes an issue
            
            #---------------------------------------------------------------------
            #obtaining next time step values    (overall solution map at 177-38)
            #---------------------------------------------------------------------
            if i<(size(U)-1):
                cy[i+1]=cy[i]+Pcy[i]*dT     #(177-52)
                if cy[i+1]>cymax:
                    cy[i+1]=cymax
                gly[i+1]=gly[i]-Dgly[i]*dT  #(177-52)
                Fepho[i+1]=Fepho[i]+dT*(Chlmin-chl[i])*Rchldecomposition*(Fepho_v[i]/(Fepho_v[i]+KFepho_v))*RFeChl  #(186-39)
                Fepho_v[i+1]=Fepho[i+1]/V       
                if Fenitrogenaseideal[i]>Fenitrogenase[i]:
                    Fenitrogenase[i+1]=Fenitrogenase[i]+dT*(Fenitrogenaseideal[i]-Fenitrogenase[i])*Rnitrogenaseproduction[i]*(Febuffer_v[i]/(Febuffer_v[i]+KFebuffer_v))     #(186-39)
                else:
                    Fenitrogenase[i+1]=Fenitrogenase[i]+dT*(Fenitrogenaseideal[i]-Fenitrogenase[i])*Rnitrogenasedegradation*(Fenitrogenase_v[i]/(Fenitrogenase_v[i]+KFenitrogenase_v))     #(186-39)
                Fenitrogenase_v[i+1]=Fenitrogenase[i+1]/V
                Febuffer[i+1]=Febuffer[i]-(Fepho[i+1]-Fepho[i])-(Fenitrogenase[i+1]-Fenitrogenase[i])   #(186-39)
                Febuffer_v[i+1]=Febuffer[i+1]/V
                chl[i+1]=Fepho[i+1]/RFeChl      #updating chlorophill

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Print
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Conversion for plot  (C: compare)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Cres=res*3600*1e15  #(fmolO2 cell-1 h-1)
    CPo2=Po2*3600*1e15  #(fmolO2 cell-1 h-1)
    #========================
    #molar C:N
    #========================
    TrueQc=Qc*ones(size(Th))+gly    #(molC cell-1)
    TrueQn=Qc*Qnc*ones(size(Th))+cy   #(molN cell-1)
    MolarCN=TrueQc/TrueQn           #(molC molN-1) molar C:N ratio of cells
    #========================
    CN2fix=N2fix*3600*1e15  #(fmolN cell-1 h-1)
    PC2H4=CN2fix*C2H4Nratio  #(fmolC2H4 cell-1 h-1)
    fmolC=TrueQc*1e15   #(fmolC cell-1)
    fmolN=TrueQn*1e15   #(fmolN cell-1)
    
    #============
    #Fe
    #============
    NFebuffer=Febuffer*Avogadro  #(Fe cell-1)
    NFepho=Fepho*Avogadro   #(Fe cell-1)
    NFenitrogenase=Fenitrogenase*Avogadro   #(Fe cell-1)
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Print parameters
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Free parameters
    #OOOOOOOOOOOOOOOOOOOO
#     print("Ychl-Fe_photo",1/RFeChl)
#     print("lmax",lsmaxmax)
#     print("Nstore_initial",cy[0])
#     print("Knstore", Kcy)
#     print("Cmax_store",glymax)
#     print("Rcstore",Rpgly)
#     print("PmaxI",Pmax)
#     print("Apho",O*Tau)
#     print("DmaxCtore",Dglymaxmax)
#     print("Kdecom_Cstore",Kgly)
#     print("P3",2)
#     print("Cpotential_O2",2e-26)
#     print("CN2fix_Fenitroge",N2fixperFe)
#     print("CPhoto-Fe_Buffer-Fe",(4*60*60)*(1/1500)*20/800)
#     print("Kfe",KFepho)
#     print("RBuffer-Fe_photo-Fe",Rchldecomposition)
#     print("RBuffer-Fe_Nitrogen_Fe",1/3000*10)
#     print("CNitroge-Fe_Buffer-Fe",(1/(4*60*60))**3*(1/1500)*3*((0.5-1)/(0.186-0.046))*(-1))
#     print("P4",3)
#     print("O2Cell_nitroge_cri",(-0.046+1/((0.5-1)/(0.186-0.046)))*(-1))
#     print("Knitroge_Cstore",Kgly*5)
#     print("Nmax_store",cymax)
#     print("N2max_fix",N2fixmax)
#     print("KN2fix_Ctore",Kgly)
    
    return nansum(PC2H4*dT/3600),nansum(Cres*dT/3600),nansum(CPo2*dT/3600)
#############################################################################

Labdata=labdata()
i=0

Tmin = 14
Tmax = 30
Tstep = 1
Tarray = arange(Tmin,Tstep+Tmax,Tstep)
U = arange(size(Tarray))
N2fix = zeros(size(Tarray))
ResArray = zeros(size(Tarray))
PhotoArray = zeros(size(Tarray))

for aa in Labdata:
   if i==2:
       count=0
       tl=time.time()
       for j in U:
           N2fix[j],ResArray[j],PhotoArray[j] = croco(aa,Tarray[j])
           print(count,'looptime =',round(time.time()-tl,2),'(s)')
           tl=time.time()   
           count=count+1
              
   i=i+1

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Plot
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rcParams.update({'font.size': 22})
rcParams.update({'lines.markersize': 10})
rcParams.update({'lines.markeredgewidth': 0.5})
rcParams.update({'font.size': 18})
rcParams.update({'lines.linewidth': 2.5})
rcParams.update({'figure.autolayout': True})
rcParams['figure.figsize']=5,4
rcParams.update({'figure.facecolor':'W'})
rcParams.update({'axes.linewidth':1.2}) 
rcParams.update({'xtick.major.width':1.2})  #
rcParams.update({'ytick.major.width':1.2})  #
rcParams.update({'lines.markeredgewidth': 1}) #default = 1

Xlabel='$T$ ($^\circ$C)'

figure(11)
plot(Tarray,N2fix)
xlabel(Xlabel)
ylabel('N$_{2}$fix (fmolC$_{2}$H$_{4}$ cell$^{-1}$ d$^{-1}$)')

savefolder='Croco\\Tdepend02'
dpi=300

figure(13)
ax1=gca(); ax2=ax1.twinx()
ln1=ax1.plot(Tarray,N2fix,label='N$_2$fix')
ln2=ax2.plot(Tarray,ResArray,color='red',label='Resp')
ln3=ax2.plot(Tarray,PhotoArray,color='green',label='Photo')
TarrayGen = genfromtxt("C:\\Users\\Keiin\\OneDrive\\Desktop\\figures\\01\\Croco\\Tdepend02\\Tall\\IarrayN2fix04.csv",delimiter=',')
NfixGen = genfromtxt("C:\\Users\\Keiin\\OneDrive\\Desktop\\figures\\01\\Croco\\Tdepend02\\Tall\\N2fix.csv",delimiter=',')
ln4=ax1.plot(TarrayGen,NfixGen,'--',dashes=(6, 5),color='b',label='N$_2$fix')
ax1.set_xlabel(Xlabel)
lns=ln1+ln2+ln3
labs = [l.get_label() for l in lns]
title('Resp only $T$ dependent')
#legend(lns, labs, loc='upper left',fontsize=21)

ax1.set_ylabel('N$_{2}$fix (fmolC$_{2}$H$_{4}$ cell$^{-1}$ d$^{-1}$)')
ax2.set_ylabel('O$_2$ fluxes (fmol cell$^{-1}$ h$^{-1}$)')
ax1.set_ylim(0,14)
ax2.set_ylim(0,160)
xlim(14,28)
Savefig3(savefolder,'ResTonly',dpi)
#Save txt"
Condition = 'TrespOnly\\'
Savetxt(Tarray,savefolder,Condition+'IarrayN2fix04')
Savetxt(N2fix,savefolder,Condition+'N2fix')
Savetxt(PhotoArray,savefolder,Condition+'Photo')
Savetxt(ResArray,savefolder,Condition+'Resp')

t1=time.time()
print('Total',round((t1-t0)/60,2),'(min)')

show()
