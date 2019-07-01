'''
Created on Oct 11, 2018
Here connetct temperature to dissolved O2 saturation
in (mmol m-3) from Benson 1984
@author: Keisuke
'''


from pylab import *

def DissolvedO2Saturation(T):
    O2Array = transpose(genfromtxt('..\\Functions\\O2saturation.csv',delimiter=','))
    
    for T0 in arange(0,41,1):
     T1 = T0 + 1
     if T>=T0 and T<T1:
         return  O2Array[1,T0] + (O2Array[1,T1]-O2Array[1,T0])*(T-T0)
     
    return 0.776*T^2 - 7.6641*T + 353.09
 