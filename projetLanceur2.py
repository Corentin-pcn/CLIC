import numpy as np

g0 = 9.80665
RE = 6378.137
OE = 6.300387486749 / (86400)
muE = 3.986005e5

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
  return np.rad2deg(np.arcsin(np.cos(i)/np.cos(phi)))

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

def Mi(A,mu): 
  if len(A) == 2:
      Mi = [0,mu/A[-1]]
  elif len(A) == 3:
      Mi = [0,0,mu/A[-1]]
  for j in np.arange(len(A)-2,-1,-1):
    Mi[j] = Mi[j+1]/A[j]
  return(Mi)

def Me(A,Mi,K):
  return([((1-A[j])/(1 + K[j]))*Mi[j] for j in np.arange(0,len(A))])

def Ms(K,Me):
  return([K[j]*Me[j] for j in np.arange(0,len(K))])

def lagrange(b0,mu,Dvreq,configs,engines):
    Aconf, bconf = [], []
    for config in configs :
      #print("\nConfiguration moteurs : ", config)
      isp = []
      kj = []
      nb_stage = len(config)
      #print(nb_stage)
      for j in range (nb_stage):
          k = engines[config[j]]["index"]
          kj.append(k)
          if j == 0:
              isp.append(engines[config[j]]["mean_ISP"])
          else:
              isp.append(engines[config[j]]["ISP"])
      omega = Oi(kj)
      linf = 1
      lsup = 10
      c = 0.5
      dv = Dvreq + 1000
      while c > 0.000001 or dv-Dvreq > 100:
          B, dV, MI = [], [], []
          for i in np.arange(linf,lsup+c,c):
              b = Bjs(i,omega,isp)
              dv = DV(b,isp)
              if (dv > Dvreq):
                  a = Ajs(b,kj)
                  mi = Mi(a,mu)
                  me = Me(a,mi,kj)
                  ms = Ms(kj,me)
                  if (min(me) > 0 and min(ms) > 0 and mi[1] <= mi[0]-mi[1] and mi[-1] <= mi[-2]-mi[-1] and mi[0] < 1500000):
                      B.append(b)
                      dV.append(dv)
                      MI.append(mi)
          ind = np.where(dV == np.amin(dV))
          ind = ind[0][0]
          if (lsup+linf)/2 < B[ind][-1]:
              linf = (linf+lsup)/2
              lsup = lsup
              c /= 2
          elif (lsup+linf)/2 > B[ind][-1]:
              lsup = (lsup+linf)/2
              linf = linf
              c /= 2
          else:
              c /= 2
      Aconf.append(sum(mi))
      bconf.append((dv,b,mi,me,ms,omega,kj,isp))
    Aconf = np.array(Aconf)
    mini = np.where(Aconf == np.amin(Aconf))
    dv,b,mi,me,ms,omega,kj,isp = bconf[int(mini[0][0])]
    print(f"\n Omega : {omega} \n Isp : {isp} \n Kj : {kj} \n Dv : {dv} m/s")
    print(f"\n Aj : {a}\n Bj : {b}  \n Mi : {mi} \n Mej : {me}\n Msj : {ms}")
    return(1)


LaunchP = {'Baikonur':45.6, 'Vandenberg':34.7, 'Kourou':5.2,'Cap Canaveral':28.5}
LaunchPchoice = 'Vandenberg'

missions ={ "1": {"zp": 410, "za": 410, "i": 51.6, "mu": 32000, "LPlat": 45.6},
           "2": {"zp": 340, "za": 340, "i": 90, "mu": 290, "LPlat": 34.7},
           "3": {"zp": 200, "za": 35786, "i": 5.2, "mu": 3800, "LPlat": 5.2},
           "4": {"zp": 567, "za": 567, "i": 97.6, "mu": 1150, "LPlat": 28.5}}

choice = ""
choice = input("Enter the mission number of your choice (1,2,3,4) : ")

zp = missions[choice]["zp"]
za = missions[choice]["za"]
i = missions[choice]["i"]
mu = missions[choice]["mu"]
LPlat = missions[choice]["LPlat"]


loss = losses(zp)
print("loss : ", loss)
azi = azimut(i, LPlat)
print("azimut : ", azi)
 
if zp==za : # circular orbit 
  vp = np.sqrt(muE / (RE + zp)) *1000  # in m/s
else : # not circular orbit 
  a = RE + (zp + za)/2
  vp = np.sqrt(muE * (2 / (RE + zp) - 1 / a))*1000

vi = OE * RE *  np.cos(np.deg2rad(LPlat)) * np.sin(np.deg2rad(azi)) * 1000
Dvreq = vp - vi + loss

print("speeds : ")
print(f"\t injection speed : {vp} m/s")
print(f"\t initial velocity : {vi} m/s")
print(f"\t Delta_V required : {Dvreq} m/s")

#Optimal staging

b0 = 2

L = lagrange(b0,mu,Dvreq,configs,engines)


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