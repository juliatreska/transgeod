#!/usr/bin/python
# -*- coding: UTF-8 -*-

import numpy as np
from transformacje import *

geo = Transformacje(model = "wgs84")
#import sys
#sys.path.append('C:\\Users\\kinga\\projekt')



plik = "wsp_inp.txt"
# odczyt z pliku: https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.genfromtxt.html
tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)

flh = np.zeros((12,3))
xygk = np.zeros((12,3))
xy2000 = np.zeros((12,2))
xy1992 = np.zeros((12,2))
neu = np.zeros((12,3))
elaz = np.zeros((12,2))

for indeks, wiersz in enumerate(tablica):
    f, l, h = geo.xyz2flh(wiersz[0], wiersz[1], wiersz[2])
    flh[indeks, 0] = f
    flh[indeks, 1] = l
    flh[indeks, 2] = h
    
    L0 = geo.lambda0(l)
    
    xgk, ygk = geo.fl2xygk(f, l, L0)
    
    xygk[indeks, 0] = xgk
    xygk[indeks, 1] = ygk
    xygk[indeks, 2] = L0
    
    x2000, y2000 = geo.u2000(xgk, ygk, L0)
    
    xy2000[indeks, 0] = x2000
    xy2000[indeks, 1] = y2000
    
    xgk92, ygk92 = geo.fl2xygk(f, l, L0 = 19)
    
    x1992, y1992 = geo.u1992(xgk92, ygk92)
    
    xy1992[indeks, 0] = x1992
    xy1992[indeks, 1] = y1992
    
    #układ neu wokoł punktu: 
    X = 3664940.500
    Y = 1409153.590
    Z = 5009571.170
    
    n, e, u = geo.neu(f, l, h, X, Y, Z)
    
    neu[indeks, 0] = n
    neu[indeks, 1] = e
    neu[indeks, 2] = u
    
    el, az = geo.neu2elaz(n, e, u)
    
    elaz[indeks, 0] = el
    elaz[indeks, 1] = az
    
dane_raport = np.hstack((tablica, flh, xy2000, xy1992, neu, elaz))
    
# zapis: https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.savetxt.html
out = np.savetxt("wsp_out.txt", dane_raport, delimiter='--- ', fmt = ['%10.3f', '%10.3f', '%10.3f', '%10.5f', '%10.5f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.8f', '%10.8f', '%10.8f', '%10.5f', '%10.5f'], header = 'konwersja współrzednych geodezyjnych \\ julia treska \n X [m] --------- Y [m] ------- Z [m] ------- fi [st] ------ lam [st] ---- h [m] ------ x2000 [m] ----- y2000 [m] ----- x92 [m] ------ y92 [m] ------ n --------- e ----------- u ------- elewacja [st] --- azymut [st]')
