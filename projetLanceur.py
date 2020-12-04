# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 19:31:40 2020

@author: 33643
"""
import numpy as np

RE = 6378.137
mu = 3.986005*(10**5)
OE = 6.300387486749/(24*60*60)
g = 9.80665

def losses(z):
    return 2.452*(10**(-3))*(z**2) + 1.051*z + 1387.5

def Oi(K):
    O = []
    for k in K:
        O.append(k/(k+1))
    return(O)

def DV(B,ISP):
    g = 9.80665
    dv = 0
    for i in range(0,len(B)):
        dv += g*ISP[i]*np.ln(B[i])
    return dv

def Bjs(bn,n):
    B = []
    B.insert(0,bn)
    for i in np.arange(0,n,1):
        bj = 
        

LaunchP = {'Baikonur':45.6, 'Vandenberg':34.7, 'Kourou':5.2,'Cap Canaveral':28.5}

LaunchPchoice = 'Vandenberg'
i = 90
za = 340
zp = 340
LPlat = LaunchP[LaunchPchoice]

azi = np.rad2deg(np.arcsin(np.cos(np.deg2rad(i))/np.cos(np.deg2rad(LPlat))))
if azi < 10**(-10):
    azi = 0

print("azimut : ",azi)

a = za + zp + 2*RE

vp = np.sqrt(mu*((2/(zp+RE))-(1/a)))
print("injection speed : ",round(vp*1000,2)," m/s")

los = losses(zp)
print("losses : ",round(los,2)," m/s")

vi = OE*RE*np.cos(LPlat)*np.sin(azi)
print("vitesse initiale : ",round(vi,2)," m/s")

Dvreq = vp - vi + los
print("Dv required : ", round(Dvreq,2), " m/s")

