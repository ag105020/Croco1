'''
Created on Feb 1, 2015
About c103
Here we introduce the influence of iron.
The update is based mostly on Kei 186
Comments----------------------------------------------------------------
Here I modify dt/2 part in nitrogenase iron production
------------------------------------------------------------------------
@author: Keisuke
'''

from pylab import *
from Labdata02 import *
from cf002_energy_calculation_n2 import *
from cf003_energy_calculation_nh4 import *
from Savefig2 import Savefig2

#font setting as always
# rc('text', usetex=True)
# rc('font', family='serif')
# rc('font', serif='Times New Roman')

###############################################################################
def croco(data):
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Getting values from data
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Datanumber=data.Tag
    Da=copy(Datanumber)
    I_const=data.Irradiance
    Days=data.Days
    light_period_length=data.Daytime  #(s)length of light period
    dark_period_length=data.Nighttime    #(s)length of dark period
    Temperature=data.Temperature    #(K) temperature of the environment
    O2_v=data.O2               #(molO2 m-3) environmental oxygen concentration (mol m-3) (post 1982)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Choosing the data to work on 
    #0:Mohr 2010, 1:GroBkopf 2012 5% O2, 
    #2:GroBkopf 2012 20% O2, 3:Dron 2012, 4: Satio 2011
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#    if (Da==4) or (Da==0) or (Da==2) or (Da==3):        #Plotting GroBkopf 2012 O2 5%
    if (Da==4) or (Da==0) or (Da==3):                    #Plotting GroBkopf 2012 O2 5% and 20%
  #  if (Da==6):                                        #Plotting all
   # if (Da==1) or (Da==2) or (Da==3) or (Da==4):         #Plotting Mohr 2010
   # if (Da==3) or (Da==4):
        return
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Time setting
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    print(I_const)
    Tmax=86400 #(s) maximum time step (177-32)
    dT=1800    #(s) time step    (177-32)
    dT=50
    T=arange(0.,Tmax+10**(-10),dT)  #(s) time array (177-32)
    U=T/dT                   #(s) U array for for loop (from a704_16)
    Th=T/3600       #(h) time 
    
    #------------------------------
    #Size related values
    #------------------------------
    RLg=1.32*10**(-6)    #(m) radius of the cell including the cell membrane
    if Datanumber==1:
        RLg=1.28*10**(-6)
    
    Rg=1/27         #(dimensionless) ratio of glycolipid layer to the cell radius
    R=RLg/(1+Rg)    #(m) radius of the cell without including the cell membrane
    V=4/3*pi*RLg**3   #(m3) volume of the cell
    Lg=R*Rg         #(m3) thickness of the cell membrane (g=glycolipid)
    RLgoriginal=1.5*10**(-6)    #(m) originally set radius of the cell including the cell membrane
    Voriginal=4/3*pi*RLgoriginal**3   #(m3) originally set volume of the cell 
    Sizeadjustment=V/Voriginal   #(dimensionless) size adjustment term
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Values to play with
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #*****************
    #play
    #*****************
    RFeChl=12/(250*55)*6     #(molFe molC chl-1) Fe in the photosystem to chol carbon ratio
    N2fixperFe=0.1            #(molN s-1 molFe-1) N2 fixation rate per nitrogenase iron (186-28) (value determination 186-36)
    Rpgly=8e-5*0.79          #(s-1) glycogen production rate (186-43)
    glymax_v=18333.333333*0.95*0.9*0.79         #(molC m-3) maximum glycogen storage (per volume) (177-26)
    Qc_v=0.18*10**6/12     #(molC m-3) biomass carbon content (storage not included) (per volume) (177-39)
    Qnc=1/6             #(molN molC-1) N to C molar ratio in the biomass except for storage (177-26)
    lsmaxmax_v=0.33*0.79    #(molC m-3 s-1) maximum biomass productioin rate (per volume)
    Kcy_v=18333/5/2/2/10            #(Nmol m-3) Half saturation of cyanophysin for biomass production (per volume)
    Kcy=Kcy_v*V          #(Nmol cell-1) Half saturation of cyanophysin for biomass production (per cell)
    cymax_v=18333*0.5*1/5*10*2*1.3   #(molN m-3) maximum storage of cyanophysin (per volume)(177-50)
    if Datanumber==1:
        cymax_v=cymax_v*(1.3*10**(-6))**3/RLg**3
    
    
    NFe=4869633         #(Fe cell-1) Number of Fe molecules in the cell
    Rchldecomposition=1/3000*1.5   #(dimensionless) Chlorophyll adjusting speed in degradation of chlorophyll (Rspeed2 in (186-35)
    epd=0.000015         #diffusivity of glycolipid layer
    epd=0.0000155
    C2H4Nratio=1.5    #(molC2H4 molN-1) conversion
    erecycle=4/3        #(dimensionless) how much more N2 fixaction can happn by recycling the electron (100% recycle leads to 4/3 times more N2 fixation)
                        #read GroBkopf 2012 for more details
    N2fixidealpowerfactor=1    #(dimensionless) powerfactor on (cymax-cy)/cymax
    startingpoint=5*3600         #(s) counting starting point
    
    #*****************
    glymax=glymax_v*V   #(molC cell-1) maximum glycogen storage (per cell) (177-26)
    Qc=Qc_v*V           #(molC cell-1) biomass carbon content (storage not included) (per cell)
    lsmaxmax=lsmaxmax_v*V    #(molC cell-1 s-1) maximum biomass production rate (per cell)
    cymax=cymax_v*V     #(molN cell-1) maximum storage of cyanophycin (per cell)(177-50)
    Chlmin=2e-19/(12/(250*55)*3.6)*10*1.6*Sizeadjustment*0.7      #(molCchl cell-1) Chlorophyll min per cell (186-40)
    N2fixmax=1e-18*Sizeadjustment*(1.3*10**(-6))**3/RLg**3*0.8          #(molN s-1 cell-1) maximum nitrogen fixation rate when cyanophycin is depleated
    Rchlproductionmax=50/(4*60*60)*(1/1500)*20/2*Sizeadjustment*2  #(dimensionless)
    Kgly=9.25e-15*Sizeadjustment #(glyC cell-1) half saturation constant for glycogen consumption (177-28)
    #+++++++++++++++++++++++++++++++++++++
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Light dark period setting
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    light_period_step=light_period_length/dT    #light period in step
    dark_period_step=dark_period_length/dT      #light period in step
    light_or_dark=copy(T)*0     #array for light or dark 1=light, 0=dark
    light_or_dark[0:(light_period_step)]=1    #
    Oneday_light_or_dark=copy(light_or_dark)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Adjustment for multiple days
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    T=arange(0.,Tmax*Days+10**(-10),dT)
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
    Do2_25=2.12*10**(-9)  #Diffusion coefficient of O2 in the water at 25C (m2/s)

    #======================================================
    #Temperature adjustment of the diffusion coefficient
    #======================================================
    if Datanumber==4:
        Do20=Do2_25*1.050507703     #(m2/s) Temperature adjusted O2 diffusivity in water at 27C 
                                    #Conversion term is based on "15 Computation of diffusion coefficient for different temperature.xlsx"
    else:
        Do20=Do2_25*1.077446314     #(m2/s) Temperature adjusted O2 diffusivity in water at 28C
                                    #Conversion term is based on "15 Computation of diffusion coefficient for different temperature.xlsx"
    Do2=Do20*r5         #Diffusion coefficient of O2 in the water (m2/s)
    #======================================================
    
    #--------------------------------
    #physiology related parameters
    #--------------------------------
    G=2                 #(dimensionless) glycogen/biomass production ratio when zero glycogen (177-26)
    E_nh4=0.383              #("molC CO2" "molC biomass"-1) CO2 production to biomass production ratio for NH4+consumption case (from 314_64)
    E_n2=0.683          #("molC CO2" "molC biomass"-1) CO2 production to biomass production ratio for N2 fixing case (from m02)
                        #This value tells how fast the difference between ideal chlorophyll and chlorophyll diminishes
                        #if 1/100, the differnce will be zero in 100 seconds if the chl concentration is not updated  
    O2_bioC_ratio=0.383     #(molO2 molCbio) the ratio of O2 production and biomass C production (from m02)
    O2_bio_ratio_nh4_balanced=0.319       #(dimensionless) the ratio of oxygn evolution to biomass production for daytime (balanced respiration +NH4+ useage) (from m02)
    O2CH=1              #(O2mol CHmol) O2 production to Carbohydrate production ratio in photosynthtesis 
    KIn=6               #(uE m-2 s-1) Half satulation of light on biomass production(167-1)(167-17)
    lsin=6              #(dimentionless) Power factoer of light influence on biomass production (167-1)(167-17)
    
    #'''''''''''''''''''''''''''''''''''
    #night time specific parameters
    #'''''''''''''''''''''''''''''''''''
    Dglymaxmax=lsmaxmax*(1+E_nh4)*1.5       #(molC m-3 s-1) maximum glycogen decomposition rate (tunable) (177-28)

    O2CH=O2CH           #(molO2 / molC CH) oxygen consumption to CH consumption ratio (177-50)(value: 72-04)
  
    #-----------------------------------
    #photosynthesis related parameters
    #-----------------------------------
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll (from a 704)
    Pmax0=7                     #(g C /(g Chl h) Maximum production rate per chlorophyll (around 6 by Cullen 1990) (from a704-16 (in Healey 1985 model))
    Pmax=Pmax0*Mchl/12/3600/55       #(mol C s-1 mol chl-1) carbon fixing rate (156-10) (156-15) for unit conversion) (from a704-16)
    O=0.01  #(m2 / umol photon) absorption cross section (from a704_16)
    Tau=5    #(s) handling time (from a704_16) (Irina: half life of RNA is less than 1 minute (30 (s)) or so)
    PI=Pmax*(1-exp(-O*I*Tau))     #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)
    
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

    Chlmaxg=0.048 #(gchl gC-1)Chlorophyll max of the cell(186-40)
    Chlmaxmol=Chlmaxg*12*55/868  #(molCchl molC-1) Chlorophyll max of the cell in mol (186-40)
    Rnitrogenasedegradation=1/3000*10 #(dimensionless) Rate constant for nitrogenase degradation (Rspeed4 in 186-36)

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
    Nleak=copy(o)   #(molN s-1 cell-1) leaking nitrogen (186-28)
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
    cy_v[0]=18333/5/2/2*0.8     #(molN m-3) original cyanophysin concentration (per volume)
    if Datanumber==0:
        cy_v[0]=18333/5/2/2*1/4
    if Datanumber==1:
        cy_v[0]=18333/5/2/2*1/4
    if Datanumber==2:
        cy_v[0]=18333/5/2/2*1/4     #(molN m-3) original cyanophysin concentration (per volume)
    if Datanumber==3:
        cy_v[0]=18333/5/2/2*1/4     #(molN m-3) original cyanophysin concentration (per volume)
    chl_v[0]=Chlmin/V    #(molCchl m-3) original chlorophyll concentration (per volume)

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #----------------------------
    #each parameter per volume
    #----------------------------
    gly_v[0]=0    #(molC m-3) original glycogen concentration    (per volume)
    if Datanumber==4:
        gly_v[0]=18333/5/2/2*0.8*3
    
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
        O2_bioC_ratio=RM_n2.O2_bioC_ratio
        O2_bio_ratio_nh4_balanced=RM_nh4.O2_bioC_ratio
        #print(E_n2,E_nh4)
        
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
            lightfactor=(I[i]/(I[i]+KIn))**lsin 
            Rchlproduction[i]=((daytime[i]+dT/2)/(4*60*60))*(1/1500)*20/800    #(186-42) learned from Rnitrogenaseproduction
            Rchlproduction[i]=Rchlproductionmax
            E[i]=E_nh4
            lsmax[i]=lsmaxmax*cy[i]/(cy[i]+Kcy/Ycnbio[i]/Qnc)     #(177-26)(167-1)(167-2)

            Pglymax[i]=(glymax-gly[i])*Rpgly   #(177-26)
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
            Rnitrogenaseproduction[i]=((nighttime[i])/(4*60*60))**3*(1/1500)*3*((0.5-1)/(0.186-0.046)*(O2c_v[i]-0.046)+1)*(gly[i]/(gly[i]+Kgly*5))**3    #(186-42)(187-7)
            #=====================
            
            resmax[i]=res[i-1]+nighttime[i]*4.5e-18/(60*60*4)/10/10        #(186-49)
            resmax[i]=nighttime[i]**2*2e-26
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
            if i==0:
                O2c_v0=O2c_v[0]
            else:
                O2c_v0=O2c_v[i-1]
            if O2c_v0<0:
                O2c_v0=0
            
            if O2c_v0>Ocri:        #updated O2 influence
                N2fixpotential[i]=0
                which2[i]=1
            else:
                N2fixpotential[i]=Fenitrogenase[i]*N2fixperFe*(Ocri-O2c_v0)/Ocri    #(187-4)(186-29)(updated O2 influence)
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
               #     if O2c_v[i]>Ocri: 
               #         N2fixideal0[i]=0
               #     else:
                    N2fixideal0[i]=((cymax-cy[i])/cymax)**N2fixidealpowerfactor*N2fixmax*gly[i]/(gly[i]+Kgly)#*(Ocri-O2c_v[i])/Ocri
                    #=====play==========
                    # N2fixideal0[i]=1*N2fixmax*gly[i]/(gly[i]+Kgly)
                    #===================
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
            #print(An2fix,8/4.28)
     #       if res[i]<Resmin[i]:
            Resoriginal[i]=copy(res[i])
            res[i]=copy(Resmin[i])
            Vo2[i]=res[i]
            Dgly[i]=Dgly[i]+res[i]-Resoriginal[i]
            O2c_v[i]=((-1)*Vo2[i])/(4*pi*R*Do2)+O2_v
            
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
    fmolCmin=Qc*ones(size(Th))*1e15      #(molC cell-1)
    fmolN=TrueQn*1e15   #(fmolN cell-1)
    fmolNmin=Qc*Qnc*ones(size(Th))*1e15  #(fmolN cell-1)
    
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
    print("Ychl-Fe_photo",1/RFeChl)
    print("lmax",lsmaxmax)
    print("Nstore_initial",cy[0])
    print("Knstore", Kcy)
    print("Cmax_store",glymax)
    print("Rcstore",Rpgly)
    print("PmaxI",Pmax)
    print("Apho",O*Tau)
    print("DmaxCtore",Dglymaxmax)
    print("Kdecom_Cstore",Kgly)
    print("P3",2)
    print("Cpotential_O2",2e-26)
    print("CN2fix_Fenitroge",N2fixperFe)
    print("CPhoto-Fe_Buffer-Fe",(4*60*60)*(1/1500)*20/800)
    print("Kfe",KFepho)
    print("RBuffer-Fe_photo-Fe",Rchldecomposition)
    print("RBuffer-Fe_Nitrogen_Fe",1/3000*10)
    print("CNitroge-Fe_Buffer-Fe",(1/(4*60*60))**3*(1/1500)*3*((0.5-1)/(0.186-0.046))*(-1))
    print("P4",3)
    print("O2Cell_nitroge_cri",(-0.046+1/((0.5-1)/(0.186-0.046)))*(-1))
    print("Knitroge_Cstore",Kgly*5)
    print("Nmax_store",cymax)
    print("N2max_fix",N2fixmax)
    print("KN2fix_Ctore",Kgly)
    
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Plot
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

    nfigure=1
    #================================
    #Mohr 2010 comparison
    #================================
    Xlabel = '$t$ (h)'
    nfigure=nfigure+1
    if Datanumber==0:
        figure(nfigure)
        plot(Th,PC2H4)
        plot(data.Mohr2010Fig1Tethylenproduction,data.Mohr2010Fig1ethylenproduction,'o',color='black')
        xlabel('T (h)')
        ylabel('Ethylene production (fmolC$_{2}$H$_{4}$ cell$^{-1}$ h$^{-1}$)')
        title('Ethylen production')
        xlim([0,24*Days])
        xticks(arange(0,24*Days+1,2*Days))

    nfigure=nfigure+1
    if Datanumber==0:
        figure(nfigure)
        plot(Th,MolarCN)
        plot(data.Mohr2010Fig1TCNratio,data.Mohr2010Fig1CNratio,'o',color='black')
        xlabel('T (h)')
        ylabel('Molar C:N ratio (molC molN-1)')
        title('Molar C:N ratio')
        xlim([0,24*Days])
        xticks(arange(0,24*Days+1,2*Days))

    #================================
    #GroBkopf 2012 comparison
    #================================
    
    color05='black'
    color20='red'
    
    color051='gray'
    color201='pink'
    savefolder='Croco\\NoResPro02'
    dpi=600
    Xlabel = '$t$ (h)'
    xedge = 0.2
    
    def dp051(filename):
        da=genfromtxt('data02\\'+filename+'.txt',)
        errorbar(da[:5,0],da[:5,1],(da[:5,5],da[:5,4]),fmt='o',color=color051,elinewidth=1.5,ecolor=color051,capthick=1.5,capsize=5)

    def dp201(filename):
        da=genfromtxt('data02\\'+filename+'.txt',)
        errorbar(da[:5,0],da[:5,1],(da[:5,5],da[:5,4]),fmt='o',color=color201,elinewidth=1,ecolor=color201,capthick=1,capsize=5)
    
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

    if (Datanumber==1) or (Datanumber==2):
        figure(1)
        if Datanumber==1:
            plot(Th,PC2H4/C2H4Nratio,color=color05)
            dp05('N2fix05all')
        if Datanumber==2:
            plot(Th,PC2H4/C2H4Nratio,color=color20)
            dp20('N2fix20all')
        xlabel(Xlabel)
        ylabel('N$_{2}$fix (fmolN cell$^{-1}$ h$^{-1}$)')
        xlim([0,24*Days+xedge])
        xticks(arange(0,24*Days+1,2*Days))
        if Datanumber==2:
            Savefig2(savefolder,1,dpi)
    
    if (Datanumber==1) or (Datanumber==2):
        figure(2)

        if Datanumber==1:
            plot(Th,CPo2,color=color05)
            plot(data.GroBkopf2012Fig3Tphotosynthesis_05O2,data.GroBkopf2012Fig3photosynthesis_05O2,'o',color=color05)
        if Datanumber==2:
            plot(Th,CPo2,color=color20)
            plot(data.GroBkopf2012Fig3Tphotosynthesis_20O2,data.GroBkopf2012Fig3photosynthesis_20O2,'o',color=color20)
        xlabel(Xlabel)
        ylabel('O$_2$ production (fmolO$_2$ cell$^{-1}$ h$^{-1}$)')
        xlim([0,24*Days+xedge])
        xticks(arange(0,24*Days+1,2*Days))
        if Datanumber==2:
            Savefig2(savefolder,2,dpi)
    
    if (Datanumber==1) or (Datanumber==2):
        figure(3)
        
        if Datanumber==1:
            plot(Th,CPo2,color=color051)
            dp051('O2fluxPho05all')
            plot(Th,Cres,color=color05)
            dp05('O2fluxRes05all')
        if Datanumber==2:
            plot(Th,CPo2,color=color201)
            dp201('O2fluxPho20all')
            plot(Th,Cres,color=color20)
            dp20('O2fluxRes20all')
        xlabel(Xlabel)
        ylabel('O$_2$ flux (fmol O$_2$ cell$^{-1}$ h$^{-1}$)')
      #  title('Respiration')
        xlim([0,24*Days+xedge])
        xticks(arange(0,24*Days+1,2*Days))
        if Datanumber==2:
            Savefig2(savefolder,3,dpi)

    if (Datanumber==1) or (Datanumber==2):
        figure(4)
       
        if Datanumber==1:
            plot(Th,fmolC,color=color05)
            dp05('Ccell05all')
        if Datanumber==2:
            plot(Th,fmolC,color=color20)
            dp20('Ccell20all')
        xlabel(Xlabel)
        ylabel('C per cell (fmol C cell$^{-1}$)')
       # title('C per cell')
        xlim([0,24*Days+xedge])
        xticks(arange(0,24*Days+1,2*Days))
        if Datanumber==2:
            Savefig2(savefolder,4,dpi)

    if (Datanumber==1) or (Datanumber==2):
        figure(5)
        
        if Datanumber==1:
            
            plot(Th,fmolN,color=color05)
            dp05('Ncell05all')
        if Datanumber==2:  
            plot(Th,fmolN,color=color20)  
            dp20('Ncell20all')
        xlabel(Xlabel)
        ylabel('N per cell (fmol N cell$^{-1}$)')
      #  title('N per cell')
        xlim([0,24*Days+xedge])
        xticks(arange(0,24*Days+1,2*Days))
        if Datanumber==2:
            Savefig2(savefolder,5,dpi)
    
    if (Datanumber==1) or (Datanumber==2):
        figure(6)
        
        if Datanumber==1:
            plot(Th,MolarCN,color=color05)
            dp05('CNcell05all')
        if Datanumber==2:  
            plot(Th,MolarCN,color=color20)
            dp20('CNcell20all')
        xlabel(Xlabel)
        ylabel('Molar C:N ratio (molC molN$^{-1}$)')
        xlim([0,24*Days+xedge])
        xticks(arange(0,24*Days+1,2*Days))
        if Datanumber==2:
            Savefig2(savefolder,6,dpi)
        
        figure(7)
        if Datanumber==1:
            plot(Th,O2c_v,color=color05)
        if Datanumber==2:
            plot(Th,O2c_v,color=color20)
        xlabel(Xlabel)
        ylabel('O$_2$ cell (molO$_2$ m$^{-3}$)')
        xlim([0,24*Days+xedge])
        xticks(arange(0,24*Days+1,2*Days))
        ylim([0,1])
        yticks(arange(0,1.01,0.2))   
        if Datanumber==2:
            Savefig2(savefolder,7,dpi)
    
    #================================
    #Dron 2012 comparison
    #================================
    

    nfigure=nfigure+1
    if Datanumber==3:
        figure(nfigure)
        plot(Th,fmolC)
    #    plot(Th,fmolCmin)
        plot(data.Dron2012Fig3TCpercell_black,data.Dron2012Fig3Cpercell_black,'o',color='black')
        plot(data.Dron2012Fig3TCpercell_gray,data.Dron2012Fig3Cpercell_gray,'o',color='gray')
        xlabel('T (h)')
        ylabel('C per cell (fmolC cell$^{-1}$)')
        title('C per cell')
        xlim([0,24*Days])
        xticks(arange(0,24*Days+1,2*Days))
    
    nfigure=nfigure+1
    if Datanumber==3:
        figure(nfigure)
        plot(Th,fmolN)
    #    plot(Th,fmolNmin)
        plot(data.Dron2012Fig3TNpercell_black,data.Dron2012Fig3Npercell_black,'o',color='black')
        plot(data.Dron2012Fig3TNpercell_gray,data.Dron2012Fig3Npercell_gray,'o',color='gray')
        xlabel('T (h)')
        ylabel('N per cell (fmolN cell$^{-1}$)')
        title('N per cell')
        xlim([0,24*Days])
        xticks(arange(0,24*Days+1,2*Days))
        
    nfigure=nfigure+1
    if Datanumber==3:
        figure(nfigure)
        plot(Th,MolarCN)
        plot(data.Dron2012Fig2TCNratioblack,data.Dron2012Fig2CNratioblack,'o',color='black')
        plot(data.Dron2012Fig2TCNratiogray,data.Dron2012Fig2CNratiogray,'o',color='gray')
        xlabel('T (h)')
        ylabel('Molar C:N ratio (molC molN$^{-1}$)')
        title('Molar C:N ratio')
        xlim([0,24*Days])
        xticks(arange(0,24*Days+1,2*Days))
    
    #====================================
    #Saito 2011 comparison
    #====================================
    
    nfigure=nfigure+1
    if Datanumber==4:
        StackPlotColors=('#00FF00','#FF0000','#000000')
        figure(nfigure)
        stackplot(Th,NFebuffer,NFepho,NFenitrogenase,colors=StackPlotColors)
        plot(data.Saito2011FigS4TNitrogenasePSIpercent,data.Saito2011FigS4NitrogenasePSIpercent*4.87e4,'o',color='#FFFFFF')
        plot(data.Saito2011FigS4TPSIbufferpercent,data.Saito2011FigS4PSIbufferpercent*4.87e4,'o',color='yellow')
        ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        xlabel('T (h)')
        ylabel('Fe (molFe cell$^{-1}$)')
        title('Cell Fe')
        xlim([0,24*Days])
        xticks(arange(0,24*Days+1,2*Days))
    
    return
#done for crocosphaera t00function
#############################################################################

a=array([25,50,75,100,150,200])
a=array([100])
#Dron average-------------
a=array([130*2/pi])
#--------------------------

Labdata=labdata()

for aa in Labdata:
    croco(aa)

show()
