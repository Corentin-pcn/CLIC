import numpy as np

g0 = 9.80665
RE = 6378.137
OE = 6.300387486749 / (86400)
mu = 3.986005e5

engines ={ "LOX-RP1": {"stage": [1,2,3], "ISP": 330, "mean_ISP": 287, "index": 0.15},
           "LOX-LK2": {"stage": [2,3], "ISP": 440, "mean_ISP": np.nan, "index": 0.22},
           "Solid": {"stage": [1], "ISP": 300, "mean_ISP": 260, "index": 0.10} }
      
configs = [
     ["Solid", "LOX-RP1"],
     ["Solid", "LOX-LK2"],
     ["LOX-RP1", "LOX-RP1"],
     ["LOX-RP1", "LOX-LK2"],
     ["Solid", "LOX-RP1", "LOX-RP1"],
     ["Solid", "LOX-RP1", "LOX-LK2"],
     ["Solid", "LOX-LK2", "LOX-RP1"],
     ["Solid", "LOX-LK2", "LOX-LK2"],
     ["LOX-RP1", "LOX-RP1", "LOX-RP1"],
     ["LOX-RP1", "LOX-RP1", "LOX-LK2"],
     ["LOX-RP1", "LOX-LK2", "LOX-RP1"],
     ["LOX-RP1", "LOX-LK2", "LOX-LK2"]
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


def Ms(K,Me):
  return([K[j]*Me[j] for j in np.arange(0,len(K))])

def lagrange(b0,mu,Dvreq):
  Aconf, bconf = [], []
  for config in configs :
    #print("\nConfiguration moteurs : ", config)
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
            mi = [0 for i in range (nb_stage)]
            me = [0 for i in range (nb_stage)]

            for k in reversed(range (nb_stage)) : 
              if k == nb_stage-1:
                mi[k] = mass_u/a[k]
              else:
                mi[k] = mi[k+1] / a[k]
           
            for k in range(nb_stage):
              me[k] = (1-a[k])/(1+kj[k])*mi[k]

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
            #print("Not possible")
            break
    Aconf.append(sum(mi))
    bconf.append((dv,b,mi,me,ms,n_iter,err))
  Aconf = np.array(Aconf)
  min = np.where(Aconf == np.amin(Aconf))
  dv,b,mi,me,ms,n_iter,err = bconf[int(min[0][0])]

  print(f"\n Omega : {omega} \n Isp : {isp} \n Kj : {kj} \n Dv : {dv} m/s  (found in {n_iter} iterations, with error = {err}). ")
  print(f"\n Aj : {a}\n Bj : {b}  \n Mi : {mi} \n Mej : {me}\n Msj : {ms}")
  return(mi,me,ms,dv)


missions ={ "1": {"zp": 410, "za": 410, "i": 51.6, "mu": 32000, "LPlat": 45.6},
           "2": {"zp": 340, "za": 340, "i": 90, "mu": 290, "LPlat": 34.7},
           "3": {"zp": 200, "za": 35786, "i": 5.2, "mu": 3800, "LPlat": 5.2},
           "4": {"zp": 567, "za": 567, "i": 97.6, "mu": 1150, "LPlat": 28.5}}

choice = ""
choice = input("Enter the mission number of your choice (1,2,3,4) : ")

zp = missions[choice]["zp"]
za = missions[choice]["za"]
i = missions[choice]["i"]
mass_u = missions[choice]["mu"]
LPlat = missions[choice]["LPlat"]


loss = losses(zp)
print("loss : ", loss)
azi = azimut(i, LPlat)
print("azimut : ", azi)
 
if zp==za : # circular orbit 
  vp = np.sqrt(mu / (RE + zp)) *1000  # in m/s
else : # not circular orbit 
  a = RE + (zp + za)/2
  vp = np.sqrt(mu * (2 / (RE + zp) - 1 / a))

vi = OE * RE *  np.cos(np.deg2rad(LPlat)) * np.sin(azi) * 1000
Dvreq = vp - vi + loss

print("speeds : ")
print(f"\t injection speed : {vp} m/s")
print(f"\t initial velocity : {vi} m/s")
print(f"\t Delta_V required : {Dvreq} m/s")

#Optimal staging

b0 = 2
L = lagrange(b0,mu,Dvreq)



"""a = 1.5
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