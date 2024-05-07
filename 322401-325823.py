import numpy as np
from math import sin, cos, tan, pi
class Transformacje:
    
    def __init__(self,  model = 'wgs84'):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
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
        
        
    def s2r(stopnie, minuty, sekundy):
        """
        Funkcja przelicza wartosc podaną w stopniach na warto
        """
        kat_stopnie = stopnie + minuty/60 + sekundy/3600
        radiany = kat_stopnie * pi/180
        return(radiany)
        
    def xyz2flh(self, x, y, z):
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
            N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
            x = (N + h)*np.cos(f)*np.cos(l)
            y = (N + h)*np.cos(f)*np.sin(l)
            z = ((N * (1 - self.e2)) + h) * np.sin(f)
            return x, y, z
    def XYZ2neu(dX, f, l):
            R = np.array([[-sin(f)*cos(l), -sin(l), cos(f)*cos(l)],
                          [-sin(f)*sin(l), cos(l), cos(f)*cos(l)],
                          [cos(f), 0, sin(f)]])
            dx = R.T @ dX
            n = dx[0]
            e = dx[1]
            u = dx[2]
            return n, e, u 
    def fl_2000(self, f, l):
        b2 = (self.a ** 2) * (1 - self.e2)
        ep2 = (self.a ** 2 - b2) / b2
        l_0 = 0
        n = 0
        if l > s2r(13, 30, 0) and l < s2r(16, 30, 0):
            l_0 += s2r(15, 0, 0)
            n += 5
        elif l > s2r(16, 30, 0) and l < s2r(19, 30, 0):
            l_0 += s2r(18, 0, 0)
            n += 6
        elif l > s2r(19, 30, 0) and l < s2r(22, 30, 0):
            l_0 += s2r(21, 0, 0)
            n += 7
        elif l > s2r(22, 30, 0) and l < s2r(25, 30, 0):
            l_0 += s2r(24, 0, 0)
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
        b2 = (self.a ** 2) * (1 - self.e2)
        ep2 = (self.a ** 2 - b2) / b2
        l_0 = 0
        n = 0
        if l > s2r(13, 30, 0) and l < s2r(16, 30, 0):
                l_0 += s2r(15, 0, 0)
                n += 5
        elif l > s2r(16, 30, 0) and l < s2r(19, 30, 0):
                l_0 += s2r(18, 0, 0)
                n += 6
        elif l > s2r(19, 30, 0) and l < s2r(22, 30, 0):
                l_0 += s2r(21, 0, 0)
                n += 7
        elif l > s2r(22, 30, 0) and l < s2r(25, 30, 0):
                l_0 += s2r(24, 0, 0)
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
        m = 0.9993
        x_1992 = x_gk * m - 5300000
        y_1992 = y_gk * m + 500000
        return(x_1992, y_1992)