# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:36:36 2022

@author: User
"""
import math as m
import numpy as np

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        if model == "wgs84":
            self.a = 6378137.0
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a #spłaszczenie
        self.ecc2 = (2 *self.flattening -self.flattening ** 2) #mimosród
        
    def xyz2flh(self, X, Y, Z):
        '''
        Przeliczenie wspolrzednych geocentrycznych X, Y, Z na wspolrzedne geodezyjne f, l, h.
        
        INPUT:
            X [float] - wspolrzedna X w metrach
            Y [float] - wspolrzedna Y w metrach
            Z [float] - wspolrzedna Z w metrach
        
        OUTPUT:
            fi [float] - wspolrzedna fi w stopniach
            lam [float] - wspolrzedna lam w stopniach
            h [float] - wspolrzedna h w metrach
        '''
        r = m.sqrt(X**2 + Y**2) # promień
        fi = m.atan(Z / (r * (1 - self.ecc2))) # pierwsze przyblizenie
        j = 0
        fi_next = fi
        while True:
            j+=1
            fi_prev = fi_next
            N = self.a / m.sqrt(1 - self.ecc2 * (m.sin(fi_prev) ** 2))
            h = (r / m.cos(fi_prev)) - N
            fi_next = m.atan(Z / (r *(1 - (self.ecc2 * (N / (N + h))))))
            fi_next = fi_next
            if abs(fi_next - fi_prev) < (0.0000001/206265): # rho'' =206265
                break
        fi = fi_prev
        lam   = m.atan(Y / X)
        N = self.a / m.sqrt(1 - self.ecc2 * (m.sin(fi) ** 2))
        h = (r / m.cos(fi)) - N
        return(m.degrees(fi), m.degrees(lam), h)
    
    def flh2xyz(self, fi, lam, h):
        '''
        Funkcja służąca przeliczeniu wspolrzednych geodezyjnych f, l, h na wspolrzedne geocentryczne XYZ.
        
        INPUT:
            fi [float] - szerokosc geodezyjna w stopniach dziesietnych
            lam [float] - dlugosc geodezyjna w stopniach dziesietnych
            h [float] - wsyokosc elipsoidalna w metrach
        OUTPUT:
            X [float] - wspolrzedna geocentryczna X w metrach
            Y [float] - wspolrzedna geocentryczna Y w metrach
            Z [float] - wspolrzedna geocentryczna Z w metrach
        '''
        fi = m.radians(fi)
        lam = m.radians(lam)
        
        N = self.a / m.sqrt(1 - self.ecc2 * (m.sin(fi) ** 2))
        X = (N + h) * m.cos(fi) * m.cos(lam)
        Y = (N + h) * m.cos(fi) * m.sin(lam)
        Z = (N * (1 - self.ecc2) + h) * m.sin(fi)
        
        return(X, Y, Z)
    
    def neu(self, fi, lam, h, X, Y, Z):
        '''
        Funkcja służąca przeliczeniu wspolrzednych geodezyjnych f, l do układu wspolrzednych topocentrycznych NEU.
        
        INPUT:
        	fi [float] - szerokosc geodezyjna poczatku ukladu w stopniach dziesietnych
        	lam [float] - dlugosc geodezyjna poczatku ukladu w stopniach dziesietnych
        	h [float] - wsyokosc elipsoidalna poczatku ukladu w metrach
        	X [float] - wspolrzedna geocentryczna X w metrach
        	Y [float] - wspolrzedna geocentryczna Y w metrach
        	Z [float] - wspolrzedna geocentryczna Z w metrach
        OUTPUT:
        	n [float] - wspolrzedna topocentryczna n w metrach
        	e [float] - wspolrzedna topocentryczna e w metrach
        	u [float] - wspolrzedna topocentryczna u w metrach
        '''
        fi = m.radians(fi)
        lam = m.radians(lam)
        
        N = self.a / m.sqrt(1 - self.ecc2 * (m.sin(fi) ** 2))
        
        Xp = (N + h) * m.cos(fi) * m.cos(lam)
        Yp = (N + h) * m.cos(fi) * m.sin(lam)
        Zp = (N * (1 - self.ecc2) + h) * m.sin(fi)
        
        XYZp = np.array([Xp, Yp, Zp])
        XYZs = np.array([X, Y, Z])
        
        XYZ = XYZs - XYZp
        XYZ = np.array(XYZ)
        	
        Rneu = np.array([[-m.sin(fi)*m.cos(lam), -m.sin(lam), m.cos(fi)*m.cos(lam)],
                         [-m.sin(fi)*m.sin(lam), m.cos(fi), m.cos(fi)*m.sin(lam)],
                         [m.cos(fi), 0, m.sin(fi)]])
        
        n, e, u = Rneu.T @ XYZ
        return(n, e, u)
    
    def neu2elaz(self, n, e, u):
        '''
        Funkcja wyznaczajaca kąt elewacji oraz kat azymutu.
        INPUT:
            n [float] - wspolrzedna topocentryczna n w metrach
            e [float] - wspolrzedna topocentryczna e w metrach
            u [float] - wspolrzedna topocentryczna u w metrach
        OUTPUT:
            el [float] - kat elewacji w stopniach dziesietnych
            az [float] - kat azymutu w stopniach dziesietnych
        '''
        #neu = np.array([n, e, u])
        
        az = m.atan2(e, n)
        el = m.asin((u) / m.sqrt(n**2 + e**2 + u**2))
        return(m.degrees(el), m.degrees(az))
    
    def lambda0(self, lam):
        L0 = (m.floor( (lam + 1.5) / 3)) *3
        return(L0)
    
    def fl2xygk(self, fi, lam, L0):
        '''
        Funkcja przeliczajaca wspolrzedne geodezyjne f, l na plaszczyzne Gaussa-Krugera.
        INPUT:
            fi [float] - szerokosc geodezyjna poczatku ukladu w stopniach dziesietnych
            lam [float] - dlugosc geodezyjna poczatku ukladu w stopniach dziesietnych
        OUTPUT:
            xgk [float] - wspolrzedna x na plaszczyznie GK w metrach
            ygk [float] - wspolrzedna y na plaszczyznie GK w metrach
            L0 [float] - poludnik osiowy w stopniach dziesietnych
        '''
        
        fi = m.radians(fi)
        lam = m.radians(lam)
        
        b2 = (self.a ** 2) * (1 - self.ecc2)
        ep2 = ((self.a ** 2) - (b2)) / (b2)
        t = m.tan(fi)
        n2 = ep2 * (m.cos(fi) ** 2)
        N = self.a / m.sqrt(1 - self.ecc2 * (m.sin(fi) ** 2))
        #sigma:
        A0 = 1 - (self.ecc2/4) - (3/64) * (self.ecc2 ** 2) - (5/256) * (self.ecc2 ** 3)
        A2 = (3/8) * (self.ecc2 + (self.ecc2 ** 2)/4 + (15/128) * (self.ecc2 ** 3))
        A4 = (15/256) * (self.ecc2 **2 + (3/4) * (self.ecc2 ** 3))
        A6 = (35/3072) * self.ecc2 **3
        si = self.a * (A0 * fi - A2 * m.sin(2 * fi) + A4 * m.sin(4 * fi) - A6 * m.sin(6 * fi))
        dL = lam - m.radians(L0)
        
        xgk = si + ((dL**2)/2)*N*m.sin(fi)*m.cos(fi)*(1+((dL**2)/12)*((m.cos(fi))**2)*(5-t**2+9*n2+4*(n2**2))+((dL**4)/360)*((m.cos(fi))**4)*(61-58*(t**2)+t**4+270*n2-330*n2*(t**2)))
        ygk = dL * N * m.cos(fi) * (1+((dL**2)/6)*((m.cos(fi))**2)*(1-t**2+n2)+((dL**4)/120)*((m.cos(fi))**4)*(5-18*(t**2)+t**4+14*n2-58*n2*(t**2)))

        return(xgk, ygk)
    
    def u1992(self, xgk, ygk):
    	'''
    	Funkcja przeliczajaca wspolrzedne x, y z plaszczyzny GK do ukladu PL1992.
    	INPUT:
    		xgk [float] - wspolrzedna x na plaszczyznie GK w metrach
    		ygk [float] - wspolrzedna y na plaszczyznie GK w metrach
    	OUTPUT:
    		x1992 [float] - wspolrzedna x w ukladzie PL1992 w metrach
    		y1992 [float] - wspolrzedna y w ukladzie PL1992 w metrach
    	'''
    	x1992 = xgk * 0.9993 - 5300000
    	y1992 = ygk * 0.9993 + 500000
    	return(x1992, y1992)
    
    def u2000(self, xgk, ygk, L0):
    	'''
    	Funkcja przeliczajaca wspolrzedne x, y z plaszczyzny GK do ukladu PL2000.
    	INPUT:
    		xgk [float] - wspolrzedna x na plaszczyznie GK w metrach
    		ygk [float] - wspolrzedna y na plaszczyznie GK w metrach
    		L0 [float] - poludnik osiowy w stopniach dziesietnych
    	OUTPUT:
    		x2000 [float] - wspolrzedna x w ukladzie PL2000 w metrach
    		y2000 [float] - wspolrzedna y w ukladzie PL2000 w metrach
    	'''
    	x2000 = xgk * 0.999923
    	y2000 = ygk * 0.999923 + (L0/3) * 1000000 + 500000
    	return(x2000, y2000)
    
    def u1992togk(self, x1992, y1992):
    	'''
    	Funkcja przeliczajaca wspolrzedne x, y z ukladu PL1992 na plaszczyzne GK.
    	INPUT:
    		x1992 [float] - wspolrzedna x w ukladzie PL1992 w metrach
    		y1992 [float] - wspolrzedna y w ukladzie PL1992 w metrach
    	OUTPUT:
    		xgk [float] - wspolrzedna x na plaszczyznie GK w metrach
    		ygk [float] - wspolrzedna y na plaszczyznie GK w metrach
    	'''
    	xgk = (x1992 + 5300000) / 0.9993
    	ygk = (y1992 - 500000) / 0.9993
    	return(xgk, ygk)
    
    def u2000togk(self, x2000, y2000, L0):
    	'''
    	Funkcja przeliczajaca wspolrzedne x, y z ukladu PL2000 na plaszczyzne GK.
    	INPUT:
    		x2000 [float] - wspolrzedna x w ukladzie PL2000 w metrach
    		y2000 [float] - wspolrzedna y w ukladzie PL2000 w metrach
    		L0 [float] - poludnik osiowy w stopniach dziesietnycH
    	OUTPUT:
    		xgk [float] - wspolrzedna x na plaszczyznie GK w metrach
    		ygk [float] - wspolrzedna y na plaszczyznie GK w metrach
    	'''
    	xgk = x2000 / 0.999923
    	ygk = (y2000 - (L0*180/m.pi/3) * 1000000 - 500000) /0.999923
    	return(xgk, ygk)
    
    def xygk2fl(self, xgk, ygk, L0):
        '''
        Funkcja przeliczajaca wspolrzedne x, y z plaszczyzny GK na wspolrzedne geodezyjne f, l.
        INPUT:
            xgk [float] - wspolrzedna x na plaszczyznie GK w metrach
            ygk [float] - wspolrzedna y na plaszczyznie GK w metrach
            L0 [float] - poludnik osiowy w stopniach dziesietnycH
        OUTPUT:
            fi [float] - szerokosc geodezyjna poczatku ukladu w stopniach dziesietnych
            lam [float] - dlugosc geodezyjna poczatku ukladu w stopniach dziesietnych
        '''
        Azero = 1 - (self.ecc2 / 4) - (3 * (self.ecc2 ** 2) / 64) - (5 * (self.ecc2 **3)/256)
        
        fi = xgk / (self.a * Azero);
        
        while 1:
            A0 = 1 - (self.ecc2/4) - (3/64) * (self.ecc2 ** 2) - (5/256) * (self.ecc2 ** 3)
            A2 = (3/8) * (self.ecc2 + (self.ecc2 ** 2)/4 + (15/128) * (self.ecc2 ** 3))
            A4 = (15/256) * (self.ecc2 **2 + (3/4) * (self.ecc2 ** 3))
            A6 = (35/3072) * self.ecc2 **3
            si = self.a * (A0 * fi - A2 * m.sin(2 * fi) + A4 * m.sin(4 * fi) - A6 * m.sin(6 * fi))
            
            fs = fi
            fi = fi + (xgk-si)/(self.a*A0)
            if abs(fi - fs) < (0.000001/206265):
                break
            
            N = self.a / m.sqrt(1 - self.ecc2 * (m.sin(fi) ** 2))
            M = (self.a * (1 - self.ecc2)) / m.sqrt((1 - self.ecc2 * (m.sin(fi)**2)) ** 3)
            fp = fi
            t1 = m.tan(fp)
            b2 = (self.a ** 2) * (1 - self.ecc2)
            ep2 = ((self.a ** 2) - (b2)) / (b2)
            n12 = ep2 * ((m.cos(fp))**2)
            
            f = fp - ((ygk**2)*t1/(2*M*N))*(1 - ((ygk**2)/(12*(N**2)))*(5+3*(t1**2)+n12-9*n12*(t1**2)-4*(n12**2))+((ygk**4)/(360*(N**4)))*(61+90*(t1**2)+45*(t1**4)))
            lam = m.radians(L0) + (ygk/(N*m.cos(fp)))*(1 - ((ygk**2)/(6*(N**2)))*(1+2*(t1**2)+n12)+((ygk**4)/(120*(N**4)))*(5+28*(t1**2)+24*(t1**4)+6*n12+8*n12*(t1**2)))
            return(m.degrees(f), m.degrees(lam))
    
    def odleglosc2D(self, Xp, Yp, Xk, Yk):
    	'''
    	Funkcja sluzaca obliczeniu odleglosci 2D.
    	INPUT:
    		Xp [float] - wspolrzedna X punktu poczatkowego w metrach
    		Yp [float] - wspolrzedna Y punktu poczatkowego w metrach
    		Xk [float] - wspolrzedna X punktu koncowego w metrach
    		Yk [float] - wspolrzedna Y punktu koncowego w metrach
    	OUTPUT:
    		odleglosc [float] - wartosc odleglosci w metrach
    	'''
    	odleglosc = m.sqrt((Xk - Xp)**2+(Yk - Yp)**2)
    	return(odleglosc)
    
    def ogledlosc3D(self, Xp, Yp, Zp, Xk, Yk, Zk):
    	'''
    	Funkcja sluzaca obliczeniu odleglosci 3D.
    	INPUT:
    		Xp [float] - wspolrzedna X punktu poczatkowego w metrach
    		Yp [float] - wspolrzedna Y punktu poczatkowego w metrach
    		Zp [float] - wspolrzedna Z punktu poczatkowego w metrach
    		Xk [float] - wspolrzedna X punktu koncowego w metrach
    		Yk [float] - wspolrzedna Y punktu koncowego w metrach
    		Zk [float] - wspolrzedna Z punktu koncowego w metrach
    	OUTPUT:
    		odleglosc [float] - wartosc odleglosci w metrach
    	'''
    	odleglosc = m.sqrt((Xp - Xk)**2 + (Yp - Yk)**2 + (Zp - Zk)**2)
    	return(odleglosc)

    def dms(self, radiany):
        '''
        Funkcja przeliczajaca wartosci z radianow na stopnie, minuty, sekundy.
        INPUT:
            radiany [float] - wartosc w radianach
        OUTPUT:
            wynik [string] - ta sama wartosc w stopniach, minutach, sekundach
        '''
        i = 1
        z = ''
        
        if radiany < 0:
            i = -1
            z = '-'
        
        radiany = abs(radiany)
        
        radiany = radiany*180/m.pi
        d = i * m.floor(radiany)
        minuty = m.floor(60*(radiany - d))
        s = (radiany - d - minuty / 60) * 3600
        wynik = f"{z} {d: .0f} {minuty: .0f} {s: 0.5f}"
        return(wynik)


if __name__ == '__main__':
    #utworzenie obiektu
    geo = Transformacje(model = "grs80")
    
    #dane XYZ geocentryczne
    X= 3664940.500
    Y= 1409153.590
    Z= 5009571.170
    
    #przeliczenie na flh
    f, l, h = geo.xyz2flh(X, Y, Z)
    print(f, l, h)
    