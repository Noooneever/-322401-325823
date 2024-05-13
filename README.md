
INSTRUKCJA KORZYSTANIA Z PLIKU TRANSFORMACJE.PY

1. DO CZEGO SŁUŻY PROGRAM?
Program powstałw celu transformacji współrzedych podabych w poliku tekstowym dla wybranej z podanych elipsoid (WGS84, GRS80, elisoida Krasowskiego).
Możliwe transformacje: XYZ -> FLH ; FLH -> XYZ ; XYZ -> NEU ; FL -> PL2000 ; FL -> PL1992.

2. JAKIE SĄ WYMAGANIA, BY PROGRAM DZIAŁAŁ NA KOMPUTERZE?
Komputer powinien mieć pobraną aplikacje Python 3.12 (ewentualnie inne wersje, nie było sprawdzane, ale powinno działać) oraz dobrze skonfigurowane 
zmienne środowiskowe.

3. DLA JAKIEGO SYSTEMU OPERACYJNEGO ZOSTAŁ NAPISANY PROGRAM?
Program został napisany i testowany tylko dla systemu Windows. 
	
4. JAK UŻYWAĆ PROGRAMU?
Aby użyć programu trzeba spełnić wymagania podane w drugim punkcie. Po spełnieniu ich należy uruchomić program w wierszu poleceń w folderze instalacyjnym programu.
Po uruchomieniu programu należy wpisać: python transformacje.py [transformacja] [jednostka] [model_elipsoidy] [plik_ze_wspolrzednymi]
Gdzie:
transformacja to: xyz2flh, plh2xyz, XYZ2neu, fl2xygk2000, fl2xygk1992
jednostka to: w przypadku zmiany z fl: dec_degree, dms   ;   w przypadku zaminy z XYZ: no_unit (domyślnie są metry)
Model_elipsoidy: wgs84, grs80, Krasowski
plik_ze_wspolrzednymi: nazwa pliku z jednostakami

5. PRZYKŁADOWE WYWOŁANIA PROGRAMU:
python transformacje.py xyz2flh dms grs84 wsp_xyz.txt
python transformacje.py flh2xyz dec_degree Krasowski wsp_flh.txt
python transformacje.py XYZ2neu no_unit grs80 wsp_xyz.txt
python transformacje.py fl2xygk2000 dec_degree Krasowski4 wsp_flh.txt
python transformacje.py fl2xygk1992 dec_degree grs80 wsp_flh.txt
	
Wprzypadku gdy komenda zostanie poprawnie napisana pojawi sie komunitak potwierdzający: 
"Program został wykonany poprawnie".

6. PRZYKŁADOWY PLIK ZE WSPÓŁRZĘDNYMI
Program został napisany dla przykładowego plików ze współrzędnymi.
Program pozbywa się wstępu, dzieli poniższe linijki ze współrzędnymi na części oraz przypisuje do odpowiednich zmiennych.
Linie ze współrzędnymi są w następującej formie: X,Y,Z.

7. BŁĘDY I NIETYPOWE ZACHOWANIA PROGRAMU:
Podczas pisania programu nie napotkaliśmy żadnych problemów.
MIŁEGO KORZYSTANIA Z PROGRAMU! :)

