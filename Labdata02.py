'''
Created on Jul 23, 2015
Here we separate date from each paper
@author: Keisuke
'''

from pylab import *

class Mohr2010:
    def __init__(self,Tag,Irradiance,Days,Daytime,Nighttime,Temperature,O2,\
                 Mohr2010Fig1Tethylenproduction,Mohr2010Fig1ethylenproduction,\
                 Mohr2010Fig1TCNratio,Mohr2010Fig1CNratio):
        #=====================================
        #common variables
        #=====================================
        self.Tag=Tag
        self.Irradiance=Irradiance
        self.Days=Days
        self.Daytime=Daytime
        self.Nighttime=Nighttime
        self.Temperature=Temperature
        self.O2=O2
        #=====================================
        
        self.Mohr2010Fig1Tethylenproduction=Mohr2010Fig1Tethylenproduction
        self.Mohr2010Fig1ethylenproduction=Mohr2010Fig1ethylenproduction
        self.Mohr2010Fig1TCNratio=Mohr2010Fig1TCNratio
        self.Mohr2010Fig1CNratio=Mohr2010Fig1CNratio

class GroBkopf2012_05O2:
    def __init__(self,Tag,Irradiance,Days,Daytime,Nighttime,Temperature,O2,\
                 GroBkopf2012Fig1TC2H2_05O2,GroBkopf2012Fig1C2H2_05O2,\
                 GroBkopf2012Fig3Tphotosynthesis_05O2,GroBkopf2012Fig3photosynthesis_05O2,\
                 GroBkopf2012Fig3Trespiration_05O2,GroBkopf2012Fig3respiration_05O2,\
                 GroBkopf2012Fig4TCpercell_05O2,GroBkopf2012Fig4Cpercell_05O2,\
                 GroBkopf2012Fig4TNpercell_05O2,GroBkopf2012Fig4Npercell_05O2,\
                 GroBkopf2012Fig4TCNratio_05O2,GroBkopf2012Fig4CNratio_05O2):

        #=====================================
        #common variables
        #=====================================
        self.Tag=Tag
        self.Irradiance=Irradiance
        self.Days=Days
        self.Daytime=Daytime
        self.Nighttime=Nighttime
        self.Temperature=Temperature
        self.O2=O2
        #=====================================

        self.GroBkopf2012Fig1TC2H2_05O2=GroBkopf2012Fig1TC2H2_05O2
        self.GroBkopf2012Fig1C2H2_05O2=GroBkopf2012Fig1C2H2_05O2
        self.GroBkopf2012Fig3Tphotosynthesis_05O2=GroBkopf2012Fig3Tphotosynthesis_05O2
        self.GroBkopf2012Fig3photosynthesis_05O2=GroBkopf2012Fig3photosynthesis_05O2
        self.GroBkopf2012Fig3Trespiration_05O2=GroBkopf2012Fig3Trespiration_05O2
        self.GroBkopf2012Fig3respiration_05O2=GroBkopf2012Fig3respiration_05O2
        self.GroBkopf2012Fig4TCpercell_05O2=GroBkopf2012Fig4TCpercell_05O2
        self.GroBkopf2012Fig4Cpercell_05O2=GroBkopf2012Fig4Cpercell_05O2
        self.GroBkopf2012Fig4TNpercell_05O2=GroBkopf2012Fig4TNpercell_05O2
        self.GroBkopf2012Fig4Npercell_05O2=GroBkopf2012Fig4Npercell_05O2
        self.GroBkopf2012Fig4TCNratio_05O2=GroBkopf2012Fig4TCNratio_05O2
        self.GroBkopf2012Fig4CNratio_05O2=GroBkopf2012Fig4CNratio_05O2

class GroBkopf2012_20O2:
    def __init__(self,Tag,Irradiance,Days,Daytime,Nighttime,Temperature,O2,\
                 GroBkopf2012Fig1TC2H2_20O2,GroBkopf2012Fig1C2H2_20O2,\
                 GroBkopf2012Fig3Tphotosynthesis_20O2,GroBkopf2012Fig3photosynthesis_20O2,\
                 GroBkopf2012Fig3Trespiration_20O2,GroBkopf2012Fig3respiration_20O2,\
                 GroBkopf2012Fig4TCpercell_20O2,GroBkopf2012Fig4Cpercell_20O2,\
                 GroBkopf2012Fig4TNpercell_20O2,GroBkopf2012Fig4Npercell_20O2,\
                 GroBkopf2012Fig4TCNratio_20O2,GroBkopf2012Fig4CNratio_20O2):

        #=====================================
        #common variables
        #=====================================
        self.Tag=Tag
        self.Irradiance=Irradiance
        self.Days=Days
        self.Daytime=Daytime
        self.Nighttime=Nighttime
        self.Temperature=Temperature
        self.O2=O2
        #=====================================

        self.GroBkopf2012Fig1TC2H2_20O2=GroBkopf2012Fig1TC2H2_20O2
        self.GroBkopf2012Fig1C2H2_20O2=GroBkopf2012Fig1C2H2_20O2
        self.GroBkopf2012Fig3Tphotosynthesis_20O2=GroBkopf2012Fig3Tphotosynthesis_20O2
        self.GroBkopf2012Fig3photosynthesis_20O2=GroBkopf2012Fig3photosynthesis_20O2
        self.GroBkopf2012Fig3Trespiration_20O2=GroBkopf2012Fig3Trespiration_20O2
        self.GroBkopf2012Fig3respiration_20O2=GroBkopf2012Fig3respiration_20O2
        self.GroBkopf2012Fig4TCpercell_20O2=GroBkopf2012Fig4TCpercell_20O2
        self.GroBkopf2012Fig4Cpercell_20O2=GroBkopf2012Fig4Cpercell_20O2
        self.GroBkopf2012Fig4TNpercell_20O2=GroBkopf2012Fig4TNpercell_20O2
        self.GroBkopf2012Fig4Npercell_20O2=GroBkopf2012Fig4Npercell_20O2
        self.GroBkopf2012Fig4TCNratio_20O2=GroBkopf2012Fig4TCNratio_20O2
        self.GroBkopf2012Fig4CNratio_20O2=GroBkopf2012Fig4CNratio_20O2

class Dron2012:
    def __init__(self,Tag,Irradiance,Days,Daytime,Nighttime,Temperature,O2,\
                 Dron2012Fig2TCNratioblack,Dron2012Fig2CNratioblack,\
                 Dron2012Fig2TCNratiogray,Dron2012Fig2CNratiogray,\
                 Dron2012Fig3TCpercell_black,Dron2012Fig3Cpercell_black,\
                 Dron2012Fig3TCpercell_gray,Dron2012Fig3Cpercell_gray,\
                 Dron2012Fig3TNpercell_black,Dron2012Fig3Npercell_black,\
                 Dron2012Fig3TNpercell_gray,Dron2012Fig3Npercell_gray):

        #=====================================
        #common variables
        #=====================================
        self.Tag=Tag
        self.Irradiance=Irradiance
        self.Days=Days
        self.Daytime=Daytime
        self.Nighttime=Nighttime
        self.Temperature=Temperature
        self.O2=O2
        #=====================================

        self.Dron2012Fig2TCNratioblack=Dron2012Fig2TCNratioblack
        self.Dron2012Fig2CNratioblack=Dron2012Fig2CNratioblack
        self.Dron2012Fig2TCNratiogray=Dron2012Fig2TCNratiogray
        self.Dron2012Fig2CNratiogray=Dron2012Fig2CNratiogray
        self.Dron2012Fig3TCpercell_black=Dron2012Fig3TCpercell_black
        self.Dron2012Fig3Cpercell_black=Dron2012Fig3Cpercell_black
        self.Dron2012Fig3TCpercell_gray=Dron2012Fig3TCpercell_gray
        self.Dron2012Fig3Cpercell_gray=Dron2012Fig3Cpercell_gray
        self.Dron2012Fig3TNpercell_black=Dron2012Fig3TNpercell_black
        self.Dron2012Fig3Npercell_black=Dron2012Fig3Npercell_black
        self.Dron2012Fig3TNpercell_gray=Dron2012Fig3TNpercell_gray
        self.Dron2012Fig3Npercell_gray=Dron2012Fig3Npercell_gray
        
class Saito2011:
    def __init__(self,Tag,Irradiance,Days,Daytime,Nighttime,Temperature,O2,\
                 Saito2011FigS4TNitrogenasePSIpercent,Saito2011FigS4NitrogenasePSIpercent,\
                 Saito2011FigS4TPSIbufferpercent,Saito2011FigS4PSIbufferpercent):

        #=====================================
        #common variables
        #=====================================
        self.Tag=Tag
        self.Irradiance=Irradiance
        self.Days=Days
        self.Daytime=Daytime
        self.Nighttime=Nighttime
        self.Temperature=Temperature
        self.O2=O2
        #=====================================

        self.Saito2011FigS4TNitrogenasePSIpercent=Saito2011FigS4TNitrogenasePSIpercent
        self.Saito2011FigS4NitrogenasePSIpercent=Saito2011FigS4NitrogenasePSIpercent
        self.Saito2011FigS4TPSIbufferpercent=Saito2011FigS4TPSIbufferpercent
        self.Saito2011FigS4PSIbufferpercent=Saito2011FigS4PSIbufferpercent

def labdata():
    
    a=genfromtxt('Data_from_multiple_sources.csv',delimiter=',')
    Tag=[0,1,2,3,4]
    Irradiance=[50,150,150,130*2/pi,150]
    Days=[1,1,1,4,2]
    Daytime=[43200,43200,43200,43200,50400]
    Nighttime=[43200,43200,43200,43200,36000]
    Kelvin=273.15
    Temperature=[28+Kelvin,28+Kelvin,28+Kelvin,28+Kelvin,27+Kelvin]
    O2=[0.01,0.046,0.186,0.186,0.189]
    
    Mohr2010c=Mohr2010(Tag[0],Irradiance[0],Days[0],Daytime[0],Nighttime[0],Temperature[0],O2[0],a[:,1],a[:,2],a[:,4],a[:,5])
    GroBkopf2012_05O2c=GroBkopf2012_05O2(Tag[1],Irradiance[1],Days[1],Daytime[1],Nighttime[1],Temperature[1],O2[1],a[:,7],a[:,8],a[:,13],a[:,14],a[:,19],a[:,20],a[:,25],a[:,26],a[:,31],a[:,32],a[:,37],a[:,38])
    GroBkopf2012_20O2c=GroBkopf2012_20O2(Tag[2],Irradiance[2],Days[2],Daytime[2],Nighttime[2],Temperature[0],O2[2],a[:,10],a[:,11],a[:,16],a[:,17],a[:,22],a[:,23],a[:,28],a[:,29],a[:,34],a[:,35],a[:,40],a[:,41])
    Dron2012c=Dron2012(Tag[3],Irradiance[3],Days[3],Daytime[3],Nighttime[3],Temperature[3],O2[3],a[:,43],a[:,44],a[:,46],a[:,47],a[:,49],a[:,50],a[:,52],a[:,53],a[:,55],a[:,56],a[:,58],a[:,59])
    Saito2011c=Saito2011(Tag[4],Irradiance[4],Days[4],Daytime[4],Nighttime[4],Temperature[4],O2[4],a[:,61],a[:,62],a[:,64],a[:,65])
    
    return(Mohr2010c,GroBkopf2012_05O2c,GroBkopf2012_20O2c,Dron2012c,Saito2011c)
    
    
    
    
    