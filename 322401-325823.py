import numpy as np
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
        else: 
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.e = np.sqrt(2 * self.flat - self.flat ** 2) 
        self.e2 = (2 * self.flat - self.flat ** 2)
        
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