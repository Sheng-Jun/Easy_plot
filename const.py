#! /usr/bin/env python
#Written on 1/12/2013 by Sheng-Jun Lin
from numpy import *

#the constants and the functions list
const = 'Qe   [C]\nqe   [esu=statC=g**0.5 cm**1.5 s**-1=erg**0.5 cm**0.5]\n\
Me   [kg]\nme   [g]\nMp   [kg]\nmp   [g]\n\
Msun [kg]\nRsun [m]\n\
msun [g]\nrsun [cm]\neV   [J]\nev   [erg]\n\
eps  [C**2 s**2 kg**-1 m**-3]\nmu   [kg m C**-2]\n\
G    [m**3 s**-2 kg**-1]\nGG   [dyn cm**2 g**-2]\n\
k    [J K**-1]\nkk   [erg K**-1]\n\
h    [J s=kg m**2 s**-1]\nh_bar[J s]\nhh   [erg s]\nhe   [eV s]\n\
c    [m/s]\ncc   [cm/s]\nAU   [m]\nau   [cm]\nLY   [m]\nly   [cm]\n\
PC   [m]\npc   [cm]\nyr_s [s]\nmole\nH\ncs10 [cm/s]'

Qe = 1.602176565*10**-19  #electron charge[C]
qe = 4.80320451*10**-10   #electron charge[esu=statC=g**0.5 cm**1.5 s**-1=erg**0.5 cm**0.5]
Me = 9.10095*10**-31    #electron mass[kg]
me = 9.10095*10**-28    #electron mass[g]
Mp = 1.6726*10**-27     #proton mass[kg]
mp = 1.6726*10**-24     #proton mass[g]
Msun = 1.9891*10**30    #solar mass[kg]
msun = 1.9891*10**33    #solar mass[g]
Rsun = 1.392*10**9      #solar radius[m]
rsun = 1.392*10**11      #solar radius[cm]
eV = 1.602*10**-19      #Joule
ev = 1.602*10**-12      #Erg
eps = 8.8542*10**-12    #permittivity[C**2 s**2 kg**-1 m**-3]
mu = 4*pi*10**-7        #permeability of free space[kg m C**-2]
G = 6.67384*10**-11       #gravitational constant[m**3 s**-2 kg**-1]
GG = 6.67384*10**-8       #gravitational constant[dyn cm**2 g**-2]
k = 1.3807*10**-23      #Boltzmann constant[J K**-1]
kk = 1.3807*10**-16     #Boltzmann constant[erg K**-1]
h = 6.62606957*10**-34   #Planck constant[J s=kg m**2 s**-1]
h_bar = 1.054571726*10**-34  #h/(2pi)[J s]
hh = 6.62606957*10**-27  #Planck constant[erg s]
he = 4.135667516*10**-15 #Planck constant[eV s]
c = 299792458.           #speed of light[m/s]
cc = 29979245800.        #speed of light[cm/s]
LY = 9.4607*10**15      #meter
ly = 9.4607*10**17      #cm
AU = 1.496*10**11       #meter
au = 1.496*10**13       #cm
PC = 3.0857*10**16      #meter
pc = 3.0857*10**18      #cm
yr_s = 60.*60.*24.*365.242 #second
mole = 6.02214129*10**23
H = 1.00794
cs10 = sqrt(412424231.1507937) #cm/s

