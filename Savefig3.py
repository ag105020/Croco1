'''
Created on May 18, 2014
This one reads dpi
@author: Keisuke
'''
from pylab import * 

def Savefig3(savefolder,fignumber,Dpi):
    First_part="C:\\Users\\keiin\\OneDrive\\Desktop\\figures\\01\\"
    Second_part=savefolder+"\\"
    Figure_number=str(fignumber)
    Last_part=".png"
    savefig(First_part+Second_part+Figure_number+Last_part,dpi=Dpi)
