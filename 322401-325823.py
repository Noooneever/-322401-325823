import numpy as np
from math import sin, cos, tan, pi, degrees
import sys
class Transformacje:
    
    def __init__(self,  model = 'wgs84'):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        """
        if  model == 'wgs84':
            self.a = 6378137.0
            self.b = 6356752.31424518
        elif model == 'grs80':
            self.a = 6378137
            self.b = 6356752.31414036
        elif model == 'Krasowski':
            self.a = 6378245.000
            self.b = 6356863.019
        else: 
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.e = np.sqrt(2 * self.flat - self.flat ** 2) 
        self.e2 = (2 * self.flat - self.flat ** 2)
        # PODAWANIE JEDNOSTKI
    unit = sys.argv[2]
    if unit == "dec_degree" or unit == "dms" or unit == "no_unit":
            pass
    else:
            raise NotImplementedError(f"Nieprawidłowa jednostka")
            
    def s2r(self, stopnie, minuty, sekundy):
                """
                Funkcja przelicza wartość podana w stopniach, minutach, sekundach na wartość w radianach.
                
                INPUT:
                    stopnie  [int] : liczba stopni
                    minuty   [int] : liczba minut
                    sekundy  [float] : liczba sekund
                
                OUTPUT:
                    radiany :[float] : wartość podana w radianach
                """
                kat_stopnie = stopnie + minuty/60 + sekundy/3600
                radiany = kat_stopnie * pi/180
                return(radiany)
        
    def xyz2flh(self, x, y, z):
            """
            Funkcja służy do transformacji współrzędnych ortokartezjańskich (prostokątnych) x, y, z 
            na współrzędne geodezyjne f, l, h.
            INPUT:
                self : rodzaj elipsoidy
                x [float] : współrzędna geocentryczna x (metry)
                y [float] : współrzędna geocentryczna y (metry)
                z [float] : współrzędna geocentryczna z (merty)
                    
            OUTPUT:
                f [float] : szerokość geodezyjna (radiany)
                l [float] : długość geodezyjna (radiany)
                h [float] : wysokość elipsoidalna (metry)
            """
            l = np.arctan2(y, x)
            p = np.sqrt(x**2 + y**2)
            f = np.arctan(z / (p * (1 - self.e2)))
            while True:
                N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
                h = p/np.cos(f) - N
                fs = f
                f = np.arctan(z / (p * (1 - (self.e2 * (N / (N + h))))))
                if np.abs(fs - f) < (0.000001 / 206265):
                    break
            return f, l, h
        
    def flh2xyz(self, f, l, h):
            """
            Funkcja służy do transformacji współrzędnych współrzędnych geodezyjnych f, l  
            na współrzędne w układzie ortokartozjańskim (prostokątnym) x, y, z
            INPUT:
                self : rodzaj elipsoidy
                f [float] : szerokość geodezyjna (radiany)
                l [float] : długość geodezyjna (radiany)
                h [float] : wysokość elipsoidalna (metry)
                    
            OUTPUT:
                x [float] : współrzędna geocentryczna x (metry)
                y [float] : współrzędna geocentryczna y (metry)
                z [float] : współrzędna geocentryczna z (metry)
            """
            N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
            x = (N + h)*np.cos(f)*np.cos(l)
            y = (N + h)*np.cos(f)*np.sin(l)
            z = ((N * (1 - self.e2)) + h) * np.sin(f)
            return x, y, z
        
    def XYZ2neu(self, f, l, h, X2, Y2, Z2):
            
            """
                    Funkcja służy do transformacji współrzędnych ortokartezjańskich (prostokątnych) x, y, z 
                    na współrzędne geodezyjne w układzie 
                    INPUT:
                    self : rodzaj elipsoidy
                    x [float] : współrzędna geocentryczna x (metry)
                    y [float] : współrzędna geocentryczna y (metry)
                    z [float] : współrzędna geocentryczna z (merty)
                    
                OUTPUT:
                    f [float] : szerokość geodezyjna (radiany)
                    l [float] : długość geodezyjna (radiany)
                    h [float] : wysokość elipsoidalna (metry)
                
                """
            R = np.array([[-np.sin(f)*np.cos(l), -np.sin(f), np.cos(f)*np.cos(l)],
                          [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                          [np.cos(f), 0, np.sin(f)]])
            X1, Y1, Z1 = self.flh2xyz(f, l, h)
            dx = np.array([X1, Y1, Z1]) - np.array([X2, Y2, Z2])
            dX = R @ dx
            n = dX[0]
            e = dX[1]
            u = dX[2]
            return n, e, u
        
    def fl_2000(self, f, l):
        """
        Funkcja służy do transformacji współrzędnych geodezyjnych f, l 
        do Państwowy Układ Współrzędnych Geodezyjnych PL-2000 x_2000, y_2000
        INPUT:
            self : rodzaj elipsoidy
            f  [float] : szerokość geodezyjna (radiany)
            l  [float] : długość geodezyjna (radiany)
                
        OUTPUT:
            x_2000 [float] : współrzędna x układu współrzędnych 2000 (metry)
            y_2000 [float] : współrzędna y układu współrzędnych 2000 (metry)
        """
        b2 = (self.a ** 2) * (1 - self.e2)
        ep2 = (self.a ** 2 - b2) / b2
        l_0 = 0
        n = 0
        if l > self.s2r(13, 30, 0) and l < self.s2r(16, 30, 0):
            l_0 += self.s2r(15, 0, 0)
            n += 5
        elif l > self.s2r(16, 30, 0) and l < self.s2r(19, 30, 0):
            l_0 += self.s2r(18, 0, 0)
            n += 6
        elif l > self.s2r(19, 30, 0) and l < self.s2r(22, 30, 0):
            l_0 += self.s2r(21, 0, 0)
            n += 7
        elif l > self.s2r(22, 30, 0) and l < self.s2r(25, 30, 0):
            l_0 += self.s2r(24, 0, 0)
            n += 8
        dlambda = l - l_0
        t = tan(f)
        eta2 = ep2 * (cos(f) ** 2)
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        A0 = 1 - (self.e2 / 4) - ((3 * self.e2 ** 2) / 64) - ((5 * self.e2 ** 3) / 256)
        A2 = (3 / 8) * (self.e2 + (self.e2 ** 2) / 4 + (15 * self.e2 ** 3) / 128)
        A4 = (15 / 256) * (self.e2 ** 2 + (3 * self.e2 ** 3) / 4)
        A6 = (35 * self.e2 ** 3) / 3072
        s = self.a * ((A0 * f) - (A2 * sin(2 * f)) +
                         (A4 * sin(4 * f)) - (A6 * sin(6 * f)))
        x_gk = s + ((dlambda ** 2 / 2) * N * sin(f) * cos(f) * (1 + (((dlambda ** 2)/12) * (cos(f) ** 2) * (5 - t ** 2 + 9 * eta2 + 4 *
                        eta2 ** 2)) + (((dlambda ** 4) / 360) * (cos(f) ** 4) * (61 - 58 * (t ** 2) + t ** 4 + 270 * eta2 - 330 * eta2 * (t ** 2)))))
        y_gk = dlambda * N * cos(f) * (1 + (((dlambda ** 2)/6) * (cos(f) ** 2) * (1 - t ** 2 + eta2)) + (
                ((dlambda ** 4) / 120) * (cos(f) ** 4) * (5 - 18 * t ** 2 + t ** 4 + 14 * eta2 - 58 * eta2 * t ** 2)))
        m = 0.999923
        x_2000 = x_gk * m
        y_2000 = y_gk * m + n * 1000000 + 500000
        return(x_2000, y_2000)
    def fl_1992(self, f, l):
        """
        Funkcja służy do transformacji współrzędnych geodezyjnych f, l 
        do Państwowy Układ Współrzędnych Geodezyjnych PL-1992 x_1992, y_1992
        INPUT:
            self : rodzaj elipsoidy
            f  [float] : szerokość geodezyjna (radiany)
            l  [float] : długość geodezyjna (radiany)
                    
        OUTPUT:
            x_1992 [float] : współrzędna x układu współrzędnych 1992 (metry)
            y_1992 [float] : współrzędna y układu współrzędnych 1992 (metry)
        """
        b2 = (self.a ** 2) * (1 - self.e2)
        ep2 = (self.a ** 2 - b2) / b2
        l_0 = self.s2r(19, 0, 0)
        dlambda = l - l_0
        t = tan(f)
        eta2 = ep2 * (cos(f) ** 2)
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        A0 = 1 - (self.e2 / 4) - ((3 * self.e2 ** 2) / 64) - ((5 * self.e2 ** 3) / 256)
        A2 = (3 / 8) * (self.e2 + (self.e2 ** 2) / 4 + (15 * self.e2 ** 3) / 128)
        A4 = (15 / 256) * (self.e2 ** 2 + (3 * self.e2 ** 3) / 4)
        A6 = (35 * self.e2 ** 3) / 3072
        s = self.a * ((A0 * f) - (A2 * sin(2 * f)) +
                         (A4 * sin(4 * f)) - (A6 * sin(6 * f)))
        x_gk = s + ((dlambda ** 2 / 2) * N * sin(f) * cos(f) * (1 + (((dlambda ** 2)/12) * (cos(f) ** 2) * (5 - t ** 2 + 9 * eta2 + 4 *
                        eta2 ** 2)) + (((dlambda ** 4) / 360) * (cos(f) ** 4) * (61 - 58 * (t ** 2) + t ** 4 + 270 * eta2 - 330 * eta2 * (t ** 2)))))
        y_gk = dlambda * N * cos(f) * (1 + (((dlambda ** 2)/6) * (cos(f) ** 2) * (1 - t ** 2 + eta2)) + (
                ((dlambda ** 4) / 120) * (cos(f) ** 4) * (5 - 18 * t ** 2 + t ** 4 + 14 * eta2 - 58 * eta2 * t ** 2)))
        m = 0.9993
        x_1992 = x_gk * m - 5300000
        y_1992 = y_gk * m + 500000
        return(x_1992, y_1992)

        
if __name__ == "__main__":

    geo = Transformacje(model = sys.argv[3])

    plik = sys.argv[4]
    metoda = sys.argv[1]
    

    with open(plik, 'r') as t:
        text = t.readlines()
        dane = text[4:]
        wsp = []
        if metoda == "xyz2flh" or metoda == "XYZ2neu":
            for linia in dane:
                elem = linia.replace(",", " ")
                x, y, z = elem.split()
                wsp.append([float(x), float(y), float(z)])
            wsp = np.array(wsp)
            X = wsp[:,0]
            Y = wsp[:,1]
            Z = wsp[:,2]
            wyniki_1 = []
            for i in range(0, len(wsp)):
                f, l, h = geo.xyz2flh(X[i], Y[i], Z[i])
                wyniki_1.append([f, l, h])
            wyniki_1 = np.array(wyniki_1)
            
        elif metoda == "flh2xyz" or metoda == "fl_2000" or metoda == "fl_1992":
            wsp_flh = []
            for linia in dane:
                elem = linia.replace(",", " ")
                f, l, h = elem.split()
                wsp_flh.append([float(f), float(l), float(h)])
            wsp_flh = np.array(wsp_flh)
            f = wsp_flh[:,0]
            l = wsp_flh[:,1]
            h = wsp_flh[:,2]
            if metoda == "flh2xyz":
                wyniki2 = []
                for i in range(0, len(wsp_flh)):
                    X1, Y1, Z1 = geo.flh2xyz(f[i]*pi/180, l[i]*pi/180, h[i])
                    wyniki2.append([X1, Y1, Z1])
                wyniki2 = np.array(wyniki2)
            else: 
                pass
        
    if metoda == "xyz2flh":
        with open('raport_xyz2flh.txt', 'w') as p:
            unit = sys.argv[2]
            if unit == "dec_degree":
                p.write('      f        |       l        |  h [m]     \n')
                for f, l, h in wyniki_1:
                    f = degrees(f)
                    l = degrees(l)
                    p.write(f'{f:.12f}  {l:.12f}  {h:.3f} \n')
            elif unit == "dms":
                p.write('        f         |       l         |  h [m]     \n')
                for f, l, h in wyniki_1:
                    f = geo.rad2dms(f)
                    l = geo.rad2dms(l)
                    p.write(f'{f} {l}   {h:.3f} \n')
                
    elif metoda == "flh2xyz":
        with open('raport_flh2xyz.txt', 'w') as p:
            p.write('    X [m]   |    Y [m]   |    Z [m]    \n')
            for X, Y, Z in wyniki2:
                p.write(f'{X:.3f} {Y:.3f} {Z:.3f} \n')
                          
    elif metoda == "fl_2000":
        wyniki_3 = []
        for f, l, h in wsp_flh:
            x_2000, y_2000 = geo.fl_2000(f*pi/180, l*pi/180)
            wyniki_3.append([x_2000, y_2000])
        wyniki_3 = np.array(wyniki_3)
        with open('raport_fl_2000.txt', 'w') as r:
            r.write('       X [m]     |       Y[m]    \n')
            for x, y in wyniki_3:
                r.write(f'{x:.9f} {y:.9f} \n')
            
    elif metoda == "fl_1992":
        wyniki_4 = []
        for f, l, h in wsp_flh:
            x_1992, y_1992 = geo.fl_1992(f*pi/180, l*pi/180)
            wyniki_4.append([x_1992, y_1992])
        wyniki_4 = np.array(wyniki_4)
        with open('raport_fl_1992.txt', 'w') as r:
            r.write('      X [m]     |      Y [m]    \n')
            for x, y in wyniki_4:
                r.write(f'{x:.2f} {y:.2f} \n')
                            
    elif metoda == "XYZ2neu":
        wyniki_5 = []
        X2 = float(input("Podaj współrzędną X punktu początkowego wektora przestrzennego: "))
        Y2 = float(input("Podaj współrzędne Y punktu początkowego wektora przestrzennego: "))
        Z2 = float(input("Podaj współrzędne Z punktu początkowego wektora przestrzennego: "))
        for f, l, h in wyniki_1:
            dx = geo.XYZ2neu(f, l, h, X2, Y2, Z2)
            wyniki_5.append(dx)
                    
        wyniki_5 = np.array(wyniki_5)
        with open('raport_XYZ2neu.txt', 'w') as f:
            f.write('  N [m]  |   E [m]  |   U [m] \n')
            for dx in wyniki_5:
                n = dx[0]
                e = dx[1]
                u = dx[2]
                f.write(f'{n:^10.4f} {e:^10.4f} {u:^10.4f} \n')
              
print("Program został wykonany poprawnie")