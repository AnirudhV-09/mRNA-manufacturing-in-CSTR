# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 12:43:42 2022

@author: Veliy
"""


import numpy as np
import scipy.integrate as scint
import xlwings
from scipy.optimize import fsolve

wb = xlwings.Book("mRNA transcripts.xlsx")
sheet = wb.sheets("ODE")
def read_parameters(boolcalc = False):
    
    RNA_0 = sheet.range("rna").value
    mg_0 = sheet.range("mgleft").value
    ntp = sheet.range("ntptotal").value
    tmax = sheet.range("tmax").value
    t7_0=sheet.range("polymerase").value
    h_0=sheet.range("hplus").value
    oh_0=sheet.range("hydroxyl").value
    mgntp_0=sheet.range("mgntp").value
    HEPES_0=sheet.range("HEPES").value
    hntp=sheet.range("hntp").value
    mg2ntp=sheet.range("mg2ntp").value
    mghntp=sheet.range("mghntp").value
    Nall=sheet.range("nall").value
   
    return RNA_0,mg_0,ntp,t7_0,h_0,oh_0,mgntp_0,HEPES_0,Nall,tmax,hntp,mg2ntp,mghntp
    
k1=sheet.range("firstk").value
k2=sheet.range("secondk").value
kapp=sheet.range("kapp").value
kac=sheet.range("kac").value
print (k1)
print (k2)
print (kapp)
print (kac)

def model(A,t,t7,oh,HEPES,nall,hntp,mg2ntp,mghntp,mgntp,mg):
    [rna,ntptot,htot]=A
    h=(10**-7.5)
    ntp=0.0001
    S0=[h,ntp,hntp,mg,mgntp,mg2ntp,mghntp]
    S0=fsolve(massbalance,S0,args=(ntptot,htot))
    vtr=kapp*(t7)*(mg*mgntp/(1+(k1*mg)+(k2*mgntp)))
    vdeg=kac*h*rna
    RNAd=vtr-vdeg
    ntpd=-1*nall*vtr
    hd=(nall-1)*vtr
    print ("model executed")
    return [RNAd,ntpd,hd]

def solve():
     rna1,mg1,ntp1,t7,h1,oh1,mgntp1,HEPES1,nall1,tm,hntp,mg2ntp,mghntp= read_parameters(boolcalc=True)
     npts = len(sheet.range("tstart:tend"))
     tend = sheet.range("tmax").value
     t = np.linspace(0, tend, npts)
     arr0=[rna1,ntp1,h1]
     solution = scint.odeint(model,arr0,t,args=(t7,oh1,HEPES1,nall1,hntp,mg2ntp,mghntp,mgntp1,mg1))
    
     rnat = solution[:,0]
     ntpt = solution[:,1]
     ht=solution[:,2]
     
    
     sheet.range("tstart").value = t.reshape((npts,1))
     sheet.range("rnastart").value = rnat.reshape((npts,1))
     sheet.range("ntpstart").value = ntpt.reshape((npts,1))
     sheet.range("hstart").value = ht.reshape((npts,1))
     print ("solve executed")
     
def massbalance(B,ntptot,htot):
    [hplus,ntpleft,hntp,mg2ntp,mghntp,mgntp,mgleft]=B
    z1= hplus*ntpleft-((10**-6.9)*hntp)
    z2=(mgleft*ntpleft)-(mgntp*10**-4.42)
    z3=(mgleft*mgntp)-(mg2ntp*10**-1.69)
    z4=(hntp*mgleft)-(mghntp*10**-1.49)
    z5=0.085-mgleft-mgntp-mghntp-(2*mg2ntp)
    z6=ntptot-ntpleft-hntp-mgntp-mg2ntp-mghntp
    z7=htot-hplus-hntp-mghntp
    return z1,z2,z3,z4,z5,z6,z7
                       
     

if __name__ == "__main__":
    parameters = np.array(read_parameters())
    boolgo = True
    while boolgo:
        val = sheet.range("A1").value
        if val == 1:
            sheet.range("A1").value = ''
            boolgo = False
        else:
           parametersnow = np.array(read_parameters())
           if np.linalg.norm(parametersnow - parameters) > 1e-10:
            parameters = parametersnow
            solve()
