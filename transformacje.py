# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:36:36 2022

@author: User
"""
import math as m

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
    
#WSZYSTKIE FUNKCJE NA KOLEJNE ZAJĘCIAAAA
#uruchom git bash, inicjalizacja, ls -all, git status, git add transformacje.py, git commit -m"dodano modul transformacje", repozytorium na githubie, dodajemy repozytorium do folderu i push