# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from numpy.linalg import inv
from tkinter import *

# =============================================================================

window = Tk()
window.title("Near Field Parameters")
window.configure(background = "#FFFFAA")



# ========================number of unknowns===================================
 
#label_1 = Label(window, text = "Number of unknowns in the system (31): ", bg = "#FFFFAA", font = "Ariel")
#nun = StringVar()
#entry_1 = Entry(window, textvariable = nun)

#label_1.grid(row = 0, sticky = E)
#entry_1.grid(row = 0, column = 1)

def execute():
#    print(nun.get())
    onbekendes = nun.get()
    
    
    return(onbekendes)



# ==========================Feed Voltage=======================================

label_2 = Label(window, text = "Antenna Feed voltage: ", bg = "#FFFFAA", font = "Ariel")
fv = StringVar()
entry_2 = Entry(window, textvariable = fv)

label_2.grid(row = 1, sticky = E)
entry_2.grid(row = 1, column = 1)

def execute1():
#    print(nun.get())
    onbekendes = fv.get()
    
    
    return(onbekendes)

# ===============================Frequency=====================================

label_3 = Label(window, text = "Frequency in GHz: ", bg = "#FFFFAA", font = "Ariel")
freQ = StringVar()
entry_3 = Entry(window, textvariable = freQ)

label_3.grid(row = 2, sticky = E)
entry_3.grid(row = 2, column = 1)

def execute2():
#    print(nun.get())
    onbekendes = freQ.get()
    
    
    return(onbekendes)

# ==============================Length half dipole=============================

#label_4 = Label(window, text = " Half of dipole length (0.2242 Lamda): ", bg = "#FFFFAA", font = "Ariel")
#lengt = StringVar()
#entry_4 = Entry(window, textvariable = lengt)

#label_4.grid(row = 3, sticky = E)
#entry_4.grid(row = 3, column = 1)

def execute3():
#    print(nun.get())
    onbekendes = lengt.get()
    
    
    return(onbekendes) 


# =============================Radius of dipole================================
 
#label_5 = Label(window, text = "Radius of Dipole (0.005 lambda): ", bg = "#FFFFAA", font = "Ariel")
#raD = StringVar()
#entry_5 = Entry(window, textvariable = raD)

#label_5.grid(row = 4, sticky = E)
#entry_5.grid(row = 4, column = 1)

def execute4():
#    print(nun.get())
    onbekendes = raD.get()
    
    
    return(onbekendes)


# =============================Wave Impedance==================================
 
#label_6 = Label(window, text = "Wave Impedance in free space(ohms): ", bg = "#FFFFAA", font = "Ariel")
#wi = StringVar()
#entry_6 = Entry(window, textvariable = wi)

#label_6.grid(row = 5, sticky = E)
#entry_6.grid(row = 5, column = 1)

def execute5():
#    print(nun.get())
    onbekendes = wi.get()
    
    
    return(onbekendes)


# =============================free space permitivity===========================

#label_7 = Label(window, text = "Free space permittivity (F/m): ")
#Eps = StringVar()
#entry_7 = Entry(window, textvariable = Eps)
#
#label_7.grid(row = 6, sticky = E)
#entry_7.grid(row = 6, column = 1)
#
#def execute6():
##    print(nun.get())
#    onbekendes = Eps.get()
#    
#    
#    return(onbekendes)
#

 
# =============================feedpoint=======================================
 
#label_8 = Label(window, text = "Feed point of the antenna (zg=0 for Center-fed): ", bg = "#FFFFAA", font = "Ariel")
#zG = StringVar()
#entry_8 = Entry(window, textvariable = zG)

#label_8.grid(row = 7, sticky = E)
#entry_8.grid(row = 7, column = 1)

def execute7():
#    print(nun.get())
    onbekendes = zG.get()
    
    
    return(onbekendes)

# =============================================================================
# =============================================================================

button_1 = Button(window, text = "Run", command = window.destroy )
button_1.grid(row = 8, column = 1)

    
    
# =============================================================================
window.mainloop()

#1
onb = execute()
gnunkns = 31
print('Unknowns: ',gnunkns)


#2
fv1 = execute1()
FV = float(fv1)
print('Feed Voltage: ',FV)

#3
freq1 = execute2()
Freq = (float(freq1) * (10**9))
print('Frequency: ',Freq)

#4
length1 = execute3()
Len = 0.2242
print('Length of dipole: ',Len)

#5
rad1 = execute4()
Rad = 0.005
print('Radius of dipole: ',Rad)

#6
wi1 = execute5()
WI = 377
print('Wave Impedance: ',WI)

##7
#eps1 = execute6()
#EPS = float(eps1)
#print('Free space permittivit: ',EPS)

#8
Zg1 = execute7()
ZG = 0
print('Feed point of the antenna (zg=0 for Center-fed): ',ZG)

Freq1 = (float(Freq/10**9))


# =============================================================================







# ==================### RMHVECTOR ###==========================================
def rmhvector(rx,ry,rz,p):
    
    rmhx = (rx[p+1]+rx[p])/2
    rmhy = (ry[p+1]+ry[p])/2
    rmhz = (rz[p+1]+rz[p])/2
    
    return(rmhx,rmhy,rmhz)
# =============================================================================
# ==================### SUNIT ###==============================================
def sunit(rx,ry,rz,p):
    
    rx1 = rx[p+1]
    ry1 = ry[p+1]
    rz1 = rz[p+1]
    rx2 = rx[p]
    ry2 = ry[p]
    rz2 = rz[p]
    rmag = np.sqrt((rx1 - rx2)**2 + (ry1 - ry2)**2 + (rz1 - rz2)**2)
    snx = (rx1 - rx2)/rmag
    sny = (ry1 - ry2)/rmag
    snz = (rz1 - rz2)/rmag

    return(snx, sny, snz, rmag)
 
# =============================================================================
# =========================##  SCALARFUN  ##===================================
def scalarfun(rx, ry, rz, wk, rad, m, n, del1):  #Needs to be finished
    
    rmx = rx[m]
    rmy = ry[m]
    rmz = rz[m]
    
    if (del1 == 0.5):
        rnx = rx[n]
        rny = ry[n]
        rnz = rz[n]
        (sx1, sy1, sz1, rmag) = sunit(rx, ry, rz, n)
        
    elif (del1 == -0.5):
        (rnx, rny, rnz) = rmhvector(rx,ry,rz,n-1)
        (sx1, sy1, sz1, rmag) = sunit(rx, ry, rz, n-1)
        
    delta2 = rmag/2.0
    
    if (m==n):
        def F(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
            return(np.exp(-1j*wk*np.sqrt((rmx-rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2))/np.sqrt((rmx - rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2))
        psi = integrate.quad(F, 0, delta2,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
      
#       F = np.exp(-1j*wk*np.sqrt((rmx-rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2))/np.sqrt((rmx - rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2)

    elif (m!=n):
        def F(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
            return(np.exp(-1j*wk*np.sqrt((rmx-rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2))/np.sqrt((rmx - rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2))
        psi = integrate.quad(F, 0, delta2,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))

    return(psi[0])
    
   
# =============================================================================
# 
# =======================## SCALARPOT ##=======================================
def scalarpot(rx, ry, rz, wk, rad, m, del1, n, q):
    
    rnx = rx[n]
    rny = ry[n]
    rnz = rz[n]
    
    rnx2 = rx[q]
    rny2 = ry[q]
    rnz2 = rz[q]
    
    (sx1, sy1, sz1, rmag) = sunit(rx, ry, rz, n)
    if (del1 == 0.5):
        (rmx, rmy, rmz) = rmhvector(rx, ry, rz, m)
        
    elif (del1 == -0.5):
        (rmx, rmy, rmz) = rmhvector(rx, ry, rz, m-1)
        
    delta = rmag
#    delta2=delta/2
    
    
    dist1 = np.sqrt((rmx - rnx)**2 + (rmy - rny)**2 + (rmz - rnz)**2)
    dist2 = np.sqrt((rmx - rnx2)**2 + (rmy - rny2)**2 + (rmz - rnz2)**2)
    
    if ((dist1 < delta) & (dist2 < delta)):       #possible problem with comparisons
#        def F(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
#            return(np.exp(-1j*wk*np.sqrt((rmx-rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2))/np.sqrt((rmx - rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2))
#        
#        psi = integrate.quad(F, 0, delta,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
        def Fi(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
            return (np.imag(np.exp(-1j*wk*np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2))/np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2)))
        def Fr(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
            return (np.real(np.exp(-1j*wk*np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2))/np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2)))
        
        psir = integrate.quad(Fr, 0.0, delta,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
        psii = integrate.quad(Fi, 0.0, delta,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
        psi=psir[0]+1j*psii[0]
      
    else:
#        def F(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
#            return(np.exp(-1j*wk*np.sqrt((rmx-rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2))/np.sqrt((rmx - rnx - s*sx1)**2 + (rmy - rny - s*sy1)**2 + (rmz - rnz - s*sz1)**2 + rad**2))
#        
#        psi = integrate.quad(F, 0, delta,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
        def Fi(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
            return (np.imag(np.exp(-1j*wk*np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2))/np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2)))
        def Fr(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
            return (np.real(np.exp(-1j*wk*np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2))/np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2)))
        
        psir = integrate.quad(Fr, 0.0, delta,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
        psii = integrate.quad(Fi, 0.0, delta,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
        psi=psir[0]+1j*psii[0]

    return(psi)
        
# =============================================================================
# ==========================## VECPOT ##=======================================

def vecpot(rx, ry, rz, wk, rad, m, n, del1):
    
    rmx = rx[m]
    rmy = ry[m]
    rmz = rz[m]
    
    if (del1 == 0.5): 
        rnx = rx[n]
        rny = ry[n]
        rnz = rz[n]
        (sx1, sy1, sz1, rmag) = sunit(rx, ry, rz, n)
    
    elif (del1 == -0.5):
        (rnx, rny, rnz) = rmhvector(rx,ry,rz,n-1)
        (sx1, sy1, sz1, rmag) = sunit(rx, ry, rz, n-1)
        
    delta2 = rmag/2.0
    
    if (m==n):
       def Fi(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
           return (np.imag(np.exp(-1j*wk*np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2))/np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2)))
       def Fr(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
           return (np.real(np.exp(-1j*wk*np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2))/np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2)))
        
       psir = integrate.quad(Fr, 0.0, delta2,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
       psii = integrate.quad(Fi, 0.0, delta2,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
       psi=psir[0]+1j*psii[0]
       
    elif (m!=n):
        def Fi(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
           return(np.imag(np.exp(-1j*wk*np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2))/np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2)))
        def Fr(s,rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad):
           return(np.real(np.exp(-1j*wk*np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2))/np.sqrt((rmx-rnx-s*sx1)**2+(rmy-rny-s*sy1)**2+(rmz-rnz-s*sz1)**2+rad**2)))
        psir = integrate.quad(Fr, 0.0, delta2,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
        psii = integrate.quad(Fi, 0.0, delta2,(rmx,rnx,wk,sx1,rmy,rny,sy1,rmz,rnz,sz1,rad))
        psi=psir[0]+1j*psii[0]
       
    return(psi)

# =========================### EFIELD ###======================================
def efieldnear(rx, ry, rz, xdd, ydd, zdd, solvector, wk, rad, factor, delta, wavelength, ix, iy, iz, nunkns):
#    rx = np.arange(nunkns2,dtype=np.float)
    
    nunkns2 = nunkns+2 
    mp1 = nunkns2 +1
    rxx= np.arange(mp1+3,dtype=np.float)
    ryy= np.arange(mp1+3,dtype=np.float)
    rzz= np.arange(mp1+3,dtype=np.float)
    for o in range (0,33):
        rxx[o]=rx[o]
        ryy[o]=ry[o]
        rzz[o]=rz[o]
    if (ix == 1):
        
        rxx[mp1] = xdd - 0.0005*wavelength
        ryy[mp1] = ydd
        rzz[mp1] = zdd
        
        rxx[mp1+1] = xdd 
        ryy[mp1+1] = ydd
        rzz[mp1+1] = zdd
        
        rxx[mp1+2] = xdd + 0.0005*wavelength
        ryy[mp1+2] = ydd
        rzz[mp1+2] = zdd
        
    elif (iy == 1):
        
        rxx[mp1] = xdd 
        ryy[mp1] = ydd - 0.0005*wavelength
        rzz[mp1] = zdd
        
        rxx[mp1+1] = xdd 
        ryy[mp1+1] = ydd
        rzz[mp1+1] = zdd
        
        rxx[mp1+2] = xdd 
        ryy[mp1+2] = ydd + 0.0005*wavelength
        rzz[mp1+2] = zdd
    
    else:
        
        rxx[mp1] = xdd 
        ryy[mp1] = ydd 
        rzz[mp1] = zdd - 0.0005*wavelength
        
        rxx[mp1+1] = xdd 
        ryy[mp1+1] = ydd
        rzz[mp1+1] = zdd
        
        rxx[mp1+2] = xdd 
        ryy[mp1+2] = ydd 
        rzz[mp1+2] = zdd + 0.0005*wavelength
        
    
    (rx1, ry1, rz1) = rmhvector(rxx,ryy,rzz,mp1+1)
    (rx2, ry2, rz2) = rmhvector(rxx,ryy,rzz,mp1)
    
    diffx = rx1 - rx2
    diffy = ry1 - ry2
    diffz = rz1 - rz2
    
    esum = 0.0
    
    for n in range(0,nunkns):
        np1 = n + 1
        
        psi1 = vecpot(rxx, ryy, rzz, wk, rad, mp1+1, np1, 0.5)
        psi2 = vecpot(rxx, ryy, rzz, wk, rad, mp1+1, np1, -0.5)
        
        psi3 = scalarpot(rxx, ryy, rzz, wk, rad, mp1+1, +0.5, np1, np1+1)
        psi4 = scalarpot(rxx, ryy, rzz, wk, rad, mp1+1, +0.5, np1-1, np1)
        psi5 = scalarpot(rxx, ryy, rzz, wk, rad, mp1+1, -0.5, np1, np1+1)
        psi6 = scalarpot(rxx, ryy, rzz, wk, rad, mp1+1, -0.5, np1-1, np1)
        
        
        # s unit vectors
        (sx1, sy1, sz1, rmag) = sunit(rx, ry, rz, np1)
        (sx2, sy2, sz2, rmag) = sunit(rx, ry, rz, n)
        
        #dotties 
        dot1 = psi1*(diffx*sx1+diffy*sy1+diffz*sz1)
        dot2 = psi2*(diffx*sx2+diffy*sy2+diffz*sz2)
        dotprod = (wk**2)*(dot1 + dot2)
        
        matrix = factor*(dotprod - psi3/delta + psi4/delta + psi5/delta - psi6/delta)
        esum = esum + matrix*solvector[n]                                               #change vector operation
        
    return( -esum/0.001/np.sqrt(2))  #x components of electric field
        




# z-comp of Electric Field (Near Field Calculation for a Wire
#Antenna) for near field calculation August, 2007

#******************** INPUT *******************************                 % Variables with asterik (*) on the comment on the right side are INPUT
nunkns = gnunkns 
feedVoltage = FV                                        # ******** Number of Unknowns (for current) on the Antenna
nunkns2 = nunkns + 2
freq = Freq                               # ******************* Frequency (Hz)
vel = 3 * 10 ** 8                                        # velocity of light in free space
omega = 2 * np.pi * freq                             # angular frequency (rad/sec)
wk = omega / vel                                     # propagation constant (2*pi/lambda)
wavelength = vel / freq                              # wavelength (m)
length = Len * wavelength                        # ****************** Half dipole length (in lambda)
dlength = 2 * length
rad = Rad * wavelength                             # ******************** Radius of Dipole (in lambda)
waveimp = WI                                       # ********** Wave Impedance in free space(ohms)
eps = 1 / (36 * np.pi * 10 ** 9 )                   # ********** Free space permittivity (F/m)
zg = ZG                                              # ********* Feed point of the antenna (zg=0 for Center-fed)
delta = 2 * length / (nunkns + 1 )                  # Antenna segment length (2*L/(N+1))


# Grid points (x,y,z) where near field needs to be computed
dx = 0
dy = 0
dz = 0.001
ndx = 1
ndy = 1
ndz = 300
xstart = 0.0 * wavelength

ystart_array = np.array([0.01, 0.03, 0.05])
zstart = 0.001 * wavelength


# ****** For Antenna case, the Dipole is Center-fed with 1 Volt ******
# Forcing Function "vmvector" is Unity at feedpoint
vmvector= np.arange(nunkns,dtype=np.float).reshape(nunkns,1)

for i in range(1,nunkns):
    matchpnt = ((nunkns - 1) / 2) + 0
    
    if (i == matchpnt):
         
        vmvector[i] = feedVoltage                 # ******* Antenna feed voltage(change the value inside the brackets)
    else:
        vmvector[i] = 0
    


# ************************** End of Input Data********************
# The array rx,ry,rz gives the x,y,z components of wire segment end points
        
rx = np.arange(nunkns2,dtype=np.float)
ry = np.arange(nunkns2,dtype=np.float)
rz = np.arange(nunkns2,dtype=np.float)

for n in range(0, nunkns2):
    rx[n] = 0.0
    ry[n] = 0.0
    rz[n] = -(length) + delta*(n-0)
    
    
     
# Calculation of the Impedance Matrix "zmatrix"
#factor = np.dtype(np.complex128)
    
factor =  -1 / (1j * 4 * np.pi * omega * eps)
zmatrix = np.zeros(shape=(nunkns,nunkns),dtype=np.complex)
for m in range(0,nunkns):
    mp1 = m + 1
    (rx1,ry1,rz1) = rmhvector(rx,ry,rz,mp1)
    (rx2,ry2,rz2) = rmhvector(rx,ry,rz,m)
    
    diffx = rx1 - rx2
    diffy = ry1 - ry2
    diffz = rz1 - rz2
    
    for l in range(0,nunkns):
        np1 = l + 1
    

        # Contribution due to vector potential
#

        psi1 = vecpot(rx, ry, rz, wk, rad, mp1, np1, 0.5)
        psi2 = vecpot(rx, ry, rz, wk, rad, mp1, np1, -0.5)

        #Contribution due to scalar potential

        psi3 = scalarpot(rx, ry, rz, wk, rad, mp1, 0.5, np1, np1 + 1)
        psi4 = scalarpot(rx, ry, rz, wk, rad, mp1, 0.5, np1 - 1, np1)
        psi5 = scalarpot(rx, ry, rz, wk, rad, mp1, -0.5, np1, np1 + 1)
        psi6 = scalarpot(rx, ry, rz, wk, rad, mp1, -0.5, np1 - 1, np1)

        # S unit vectors
        (sx1, sy1, sz1, rmag) = sunit(rx, ry, rz, np1)
        (sx2, sy2, sz2, rmag) = sunit(rx, ry, rz, l)

#        # Dot products
        dot1 = psi1 * (diffx * sx1 + diffy * sy1 + diffz * sz1)
        dot2 = psi2 * (diffx * sx2 + diffy * sy2 + diffz * sz2)
        dotprod = wk ** 2 * (dot1 + dot2)
#
        ytt=factor * (dotprod - psi3 / delta + psi4 / delta + psi5 / delta - psi6 / delta)
        zmatrix[m][l] = factor * (dotprod - psi3 / delta + psi4 / delta + psi5 / delta - psi6 / delta)
#
#    

#
## SOLUTION BY MATRIX INVERSION
#
solstep = inv(zmatrix) * vmvector                 # solvector is the solution vector
solvector=solstep[15]
##with current on the antenna
## Plot the Current Distribution on the Wire Antenna
#
rsolvec = np.real(solvector)
isolvec = np.imag(solvector)
rrealpart = np.arange(nunkns,dtype=np.float)
realpart = np.arange(nunkns2,dtype=np.float)
imagpart = np.arange(nunkns2,dtype=np.float)
iimagpart = np.arange(nunkns,dtype=np.float)
zpart = np.arange(nunkns,dtype=np.float)
for ip in range (0,nunkns):
    rrealpart[ip] = rsolvec[ip]
    iimagpart[ip] = isolvec[ip]
    zpart[ip] = rz[ip]
#end
#
for ipp in range(0,nunkns2-1):
    if ipp == 0:
        realpart[ipp] = 0
        imagpart[ipp] = 0
    elif ipp == nunkns2:
        realpart[ipp] = 0
        imagpart[ipp] = 0
    else:
        realpart[ipp] = rrealpart[ipp - 1]
        imagpart[ipp] = iimagpart[ipp - 1]
realpart[nunkns2-1] = 0
imagpart[nunkns2-1] = 0    

#
#
## Near Field Computations













mp1 = nunkns2 + 1
nearfield = np.zeros(shape=(3,ndz),dtype=np.float)
xd = np.arange(ndx,dtype=np.float)
yd = np.arange(ndy,dtype=np.float)
zd = np.arange(ndz,dtype=np.float)
for kk in range(0,3):
    nfield_points = 0
    ystart = ystart_array[kk]
    for i in range(0,ndx):
        xd[i] = xstart + (i - 1) * dx
        xdi = xd[i]
        for j in range(0,ndy):
            yd[j]= ystart + (j - 1) * dy
            ydj = yd[j]
            for k in range(0,ndz):
                nfield_points = nfield_points + 1
                zd[k] = zstart + (k - 0) * dz
                zdk = zd[k]
                # efieldx,efieldy and efieldz are the x,y,z components of the
                #Electric Field

                # **** NOTE ***** In this example only the z-comp of E-field is
                #being  
                # therefore, efieldx and efieldy are commented out

                #efieldx=efieldnear(rx,ry,rz,xdi,ydj,zdk,solvector,wk,rad,factor,delta,
                #wavelength,1,0,0,nunkns);

                #efieldy=efieldnear(rx,ry,rz,xdi,ydj,zdk,solvector,wk,rad,factor,delta,
                #wavelength,0,1,0,nunkns);

                efieldz = efieldnear(rx, ry, rz, xdi, ydj, zdk, solvector, wk, rad, factor, delta, wavelength, 0, 0, 1, nunkns)

                # nearfield(kk,nfield_points)=abs(efieldx); % if x-comp of near
                #field is needed
                # nearfield(kk,nfield_points)=abs(efieldy); % if y-comp of near
                #field is needed
                nearfield[kk][nfield_points-1] = abs(efieldz)            # if z-comp of near
                newnearfield=nearfield.astype(str)
                
                #field is needed
#-0.145305793270531 - 0.069104594965050i
#            end        # loop over k (z values)
#        end    # loop for j (y values)
#    end# loop over i (x values)
#end# loop over kk
#
#
#
#If following three lines are commented out The Antenna Current is
#not plotted
plt.plot(rz,realpart,'k-',rz,imagpart,'o')
plt.grid()
plt.title('Current Distribution over a lamda/2 dipole @ %f GHz' %(Freq1))
plt.xlabel('Length (m)')
plt.ylabel('Current Distribution (A)')

#
#
## Plot of Electric Field vs z/lambda for three values of rho/lambda
## % % 
zd_to_meter=zd*wavelength

plt.figure(0)
plt.grid()
plt.plot(zd_to_meter, nearfield[0],'k-')
plt.plot(zd_to_meter, nearfield[1],'b-')
plt.plot(zd_to_meter, nearfield[2],'r-')
plt.title('Z-component of Electric Field for lamda/2 dipole at: %f GHz' %(Freq1))
plt.xlabel('Z (m)')
plt.ylabel('Electric Field in (V/m)')

print(zd)
print("hmmmmmmmmmmmmmmm")
print(newnearfield)
#plt.plot(zd(mslice[:]), nearfield(2, mslice[:]), mstring('b-'))
#hold(mstring('on'))
#grid(mstring('on'))

#
#plot(zd(mslice[:]), nearfield(3, mslice[:]), mstring('r-'))

#xlabel(mstring('Z/lambda'))
#ylabel(mstring('Electric Field in (V/m)'))
#gtext(mstring('y (lambda) =0.01'))
#gtext(mstring('0.03'))
#gtext(mstring('0.05'))
#hold(mstring('off'))


