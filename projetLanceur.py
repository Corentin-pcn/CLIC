import numpy as np

RE = 6378.137
mu = 3.986005*(10**5)
OE = 6.300387486749/(24*60*60)
g = 9.80665

#format matrices caracteristique : [1,2,3]

engines ={ "LOX-RP1": {"stage": [1,2,3], "ISP": 330, "mean_ISP": 287, "index": 0.15},
           "LOX-LK2": {"stage": [2,3], "ISP": 440, "mean_ISP": np.nan, "index": 0.22},
           "Solid": {"stage": [1], "ISP": 300, "mean_ISP": 260, "index": 0.10} }
      
configs = [
    ["Solid", "LOX-RP1"],
#    ["Solid", "LOX-LK2"],
#    ["LOX-RP1", "LOX-RP1"],
#    ["LOX-RP1", "LOX-LK2"],
#    ["Solid", "LOX-RP1", "LOX-RP1"],
#    ["Solid", "LOX-RP1", "LOX-LK2"],
#    ["Solid", "LOX-LK2", "LOX-RP1"],
#    ["Solid", "LOX-LK2", "LOX-LK2"],
#    ["LOX-RP1", "LOX-RP1", "LOX-RP1"],
#    ["LOX-RP1", "LOX-RP1", "LOX-LK2"],
#    ["LOX-RP1", "LOX-LK2", "LOX-RP1"],
#    ["LOX-RP1", "LOX-LK2", "LOX-LK2"]
]

def losses(z):
    return 2.452*(10**(-3))*(z**2) + 1.051*z + 1387.5
  
def azimut(i,phi):
  i = np.deg2rad(i)
  phi = np.deg2rad(phi)
  return np.arcsin(np.cos(i)/np.cos(phi))

def Oi(K):
    return([k/(k+1) for k in K])

def DV(B,ISP):
    g = 9.80665
    dv = 0
    for j in range(0,len(B)):
      dv += g*ISP[j]*np.log(B[j])
    return dv

def Bjs(bn,O,ISP):
    if len(O) == 2:
      B = [0,bn]
    elif len(O) == 3:
      B = [0,0,bn]
    for j in np.arange(len(O)-2,-1,-1):
        B[j] = 1/O[j]*(1 - (ISP[j+1]/ISP[j])*(1-O[j+1]*B[j+1]))
    return(B)        

def Ajs(B,K):
  return([((1 + K[j])/B[j])-K[j] for j in np.arange(0,len(K))])

def Mi(A,K,mu,O): 
  if len(O) == 2:
      Mi = [0,mu/A[-1]]
  elif len(O) == 3:
      Mi = [0,0,mu/A[-1]]
  for j in np.arange(len(O)-2,-1,-1):
    Mi[j] = Mi[j+1]/A[j]
  return(Mi)

def Me(A,Mi,K,Oi):
  return([((1-A[j])/(1 + K[j]))*Mi[j] for j in np.arange(0,len(Oi))])

def Ms(K,Me):
  return([K[j]*Me[j] for j in np.arange(0,len(K))])




def lagrange(b0,mu,Dvreq):
  for config in configs :
    print("\nConfiguration moteurs : ", config)
    isp = []
    kj = []
    nb_stage = len(config)
    for j in range (nb_stage):
      k = engines[config[j]]["index"]
      kj.append(k)
      if j == 0:
        isp.append(engines[config[j]]["mean_ISP"])
      else:
        isp.append(engines[config[j]]["ISP"])
      
    omega = Oi(kj)
    b = Bjs(b0,omega,isp)
    dv = DV(b,isp)
    
    n_iter = 0
    lim_sup = 10
    lim_inf = 1
    
    while True : 
        n_iter += 1
        b = Bjs(b0,omega,isp)
        dv = DV(b,isp)
        err = np.abs(Dvreq - dv)
        
        if err < 0.1 and dv != np.nan : 
            a = Ajs(b,kj)
            mi = Mi(a,kj,mu,omega)
            me = Me(a,mi,kj,omega)
            ms = Ms(kj,me)
            break
        
        else:
            if Dvreq > dv:
                tmp = (b0 + lim_sup)/2
                lim_inf = b0
                b0 = tmp
            elif Dvreq < dv:
                tmp = (b0 + lim_inf) / 2
                lim_sup = b0
                b0 = tmp

        if n_iter > 1000:
            print("Not possible")
            break
    
  print(f"\n Omega : {omega} \n Isp : {isp} \n Kj : {kj} \n Dv : {dv} m/s")
  print(f"\n Aj : {a}\n Bj : {b}  \n Mi : {b} \n Mej : {me}\n Msj : {ms}")
  return(mi,me,ms,dv)



LaunchP = {'Baikonur':45.6, 'Vandenberg':34.7, 'Kourou':5.2,'Cap Canaveral':28.5}

LaunchPchoice = 'Vandenberg'
i = 90
za = 340
zp = 340
mu = 290
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

vi = OE*RE*np.cos(np.deg2rad(LPlat)*np.sin(np.deg2rad(azi)))
print("vitesse initiale : ",round(vi,2)," m/s")

Dvreq = vp - vi + los
print("Dv required : ", round(Dvreq,2), " m/s")

#Optimal Staging

b0 = 2
L = lagrange(b0,mu,Dvreq)

"""
a = 1.5
b = 10.0
c = 0.5
Dv = 10000
while (Dv-Dvreq >= 100 and c > 0.0001):
  for i in np.arange(a,b,c):
    M, D, I = [], [], []
    mi,Dv = lagrange(i,mu)
    if Dv > Dvreq:
      M.append(mi)
      D.append(Dv)
      I.append(i)
      print("mass:", M)
  mini = min(M)
  j = I[M.index(mini)]
  if j > ((b+a)/2):
    a, b, c = (b+a)/2, b, c/2
  else:
    a, b, c = a, (b+a)/2, c/2

print(mi)
print(Dv)"""









