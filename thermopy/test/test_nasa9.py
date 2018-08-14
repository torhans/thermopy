# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 17:34:44 2015

@author: monteiro
"""
from thermopy import nasa9polynomials as nasa9
from thermopy.iapws import Water
from numpy import array, dot , linspace , stack , ones, zeros


def nist_enthalpy_enginering(Tlow,Thigh,Nt,coefs):
    temp= linspace(Tlow,Thigh,Nt)
    T2=temp.reshape(-1,1)/1000
    Ta= stack((T2,T2**2/2,T2**3/3,T2**4/4,-1/T2, ones(T2.shape), zeros(T2.shape),zeros(T2.shape)),axis=-1)
    return( dot(Ta,coefs))

def test_enthalpy_tests():
    """Test for various elements enthalpies checked against a literature
    source. Relative error <= 1%.\n
    RE = abs(delta_cp_burcat - delta_cp_lit) / delta_cp_burcat
    """
    # Relative Error
    RE = 1/100
    # Initialization
    database = nasa9.Database()

    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    
    # 298 K to 1000 K
    h2_nist_cc = [33.066178,-11.363417,11.432816,-2.772874,-0.158558,-9.980797,172.707974,0.0]
    
    hydrogen = database.set_compound('H2')
    delta_h_nasa9 = (hydrogen.enthalpy(700)
                      - hydrogen.enthalpy(600))
    h_nist=nist_enthalpy_enginering(600,700,2,h2_nist_cc)
    
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 
    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE

    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    
    # 100 K to 700 K
    o2_nist_cc = [31.32234,-20.23531,57.86644,-36.50624,-0.007374,-8.903471,246.7945,0.0]
    oxygen = database.set_compound('O2')
    delta_h_nasa9 = oxygen.enthalpy(400) - oxygen.enthalpy(350)

    h_nist=nist_enthalpy_enginering(350,400,2,o2_nist_cc)
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 

    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE

    
    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    h2o_nist_cc= [-203.6060,1523.290,-3196.413,2474.455,3.855326,-256.5478,-488.7163,-285.8304]
    
    water = database.set_compound('H2O(L)')
    delta_h_nasa9 = water.enthalpy(400) - water.enthalpy(350)

    h_nist=nist_enthalpy_enginering(350,400,2,h2o_nist_cc)
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 

    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE

    
    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    # 298 K to 1200 K
    no_nist_cc=[23.83491,12.58878,-1.139011,-1.497459,0.214194,83.35783,237.1219,90.29114]
    
    no = database.set_compound('no')
    delta_h_nasa9 = no.enthalpy(1450) - no.enthalpy(400)
    
    h_nist=nist_enthalpy_enginering(400,1450,2,no_nist_cc)
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 

    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE

    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    # 298 K to 6000 K
    ar_nist_cc=[20.78600,2.825911e-7,-1.464191e-7,1.092131e-8,-3.661371e-8,-6.197350,179.9990,0.000000]
    
    ar = database.set_compound('Ar')
    delta_h_nasa9 = ar.enthalpy(1100) - ar.enthalpy(300)

    h_nist=nist_enthalpy_enginering(300,1100,2,ar_nist_cc)
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 

    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE

    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    # 298 K to 1000 K
    h2_nist_cc=[33.066178,-11.363417,11.432816,-2.772874,-0.158558,-9.980797,172.707974,0.0]

    h2 = database.set_compound('H2')
    delta_h_nasa9 = h2.enthalpy(250) - h2.enthalpy(200)

    h_nist=nist_enthalpy_enginering(200,250,2,h2_nist_cc)
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 

    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE

    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    # 298 K to 2327 K
    al2o3_nist_cc=[102.4290,38.74980,-15.91090,2.628181,-3.007551,-1717.930,146.9970,-1675.690]

    al2o3 = database.set_compound('Al2O3(a)')
    delta_h_nasa9 = al2o3.enthalpy(2200) - al2o3.enthalpy(300)
    
    h_nist=nist_enthalpy_enginering(300,2200,2,al2o3_nist_cc)
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 

    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE
    
    
    # CRN(S) CHROMIUM NITRIDE CONDENSED
    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    # 298 K to 600 K
    crn_nist_cc=[-2256.160,10531.10,-17686.50,10370.30,41.42940,362.1770,-4904.890,-117.1520]
    
    crn = database.set_compound('CrN(cr)')
    delta_h_nasa9 = crn.enthalpy(600) - crn.enthalpy(300)

    h_nist=nist_enthalpy_enginering(300,600,2,crn_nist_cc)
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 

    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE
    
    # FeCL2(L)
    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    # 950 K to 2000 K
    fecl2l_nist_cc=[102.1733,-1.078819e-8,8.424401e-9,-2.096737e-9,-1.846362e-10,-341.7989,263.5225,-311.3365]
    
    fecl2l = database.set_compound('FeCL2(L)')
    delta_h_nasa9 = fecl2l.enthalpy(2000) - fecl2l.enthalpy(950)
    
    h_nist=nist_enthalpy_enginering(950,2000,2,fecl2l_nist_cc)
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 

    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE
    
    
    # FeS(L)
    # D.R. Burgess, "Thermochemical Data" in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J.  Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, doi:10.18434/T4D303, (retrieved August 12, 2018)
    # 1463 K to 3800 K
    fesl_nist_cc=[62.55080,0.000002,-6.720620e-7,6.411921e-8,-4.303011e-7,-84.51170,166.2660,-68.81140]

    fesl = database.set_compound('FeS(L)')
    delta_h_nasa9 = fesl.enthalpy(3800) - fesl.enthalpy(1465)
    
    h_nist=nist_enthalpy_enginering(1465,3800,2,fesl_nist_cc)
    delta_h_literature = (h_nist[1]-h_nist[0]) * 1e3 

    re = abs(delta_h_nasa9 - delta_h_literature) / delta_h_nasa9
    assert re < RE


def match_elements(compound,elements):
    ret=True
    for celem,telem in zip(compound.elements,elements):
        ret &= str(celem[0]).upper()==str(telem[0]).upper()
        ret &= celem[1]==telem[1]
    if not ret:
        print(compound.elements,elements)
    return ret
    
def test_elements_tests():

    database = nasa9.Database()
    
    fesl = database.set_compound('FeS(L)')
    assert match_elements(fesl,(('Fe',1),('S',1)))
    
    al2so3 = database.set_compound('AL2O3(a)')
    assert match_elements(al2so3,(('Al',2),('O',3)))

    cap = database.set_compound('Ca+')
    assert match_elements(cap,(('Ca',1),('E',-1)))

    dm = database.set_compound('D-')
    assert match_elements(dm,(('D',1),('E',1)))


def test_reactions_test():
    
    database = nasa9.Database()
    
    N2O=database.set_compound("N2O")
    N2=database.set_compound("N2")
    O2=database.set_compound("O2")
    
    prods=(O2,N2)
    reacts=(N2O,)
    reacts_coefs=(2,)
    prods_coefs=(1,2)
    
    reaction = nasa9.Reaction(1000, reacts, prods, reacts_coefs,prods_coefs)
    print(reaction.equilibrium_constant)
   
    assert abs(reaction.equilibrium_constant() - 1.9675138130075416e+16)/1.9675138130075416e+16 < 0.01
    
