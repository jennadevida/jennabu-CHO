# -*- coding: utf-8 -*-
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt

Hf_29815_H2O = -241.83  # this is Hf (kJ/(mol*K))
s_29815_H2O = 188.84  # (J/(mol*K))

Hf_29815_CH4 = -74.87  # this is Hf (kJ/(mol*K))
s_29815_CH4 = 186.25  # (J/(mol*K)) P = 1bar

Hf_29815_CO = -110.53  # this is Hf (kJ/(mol*K))
s_29815_CO = 197.66  # (J/(mol*K)) P = 1bar

Hf_29815_CO2 = -393.52  # this is Hf (kJ/(mol*K))
s_29815_CO2 = 213.79  # (J/(mol*K)) P = 1bar

Hf_29815_C2H2 = 226.73  # this is Hf (kJ/(mol*K))
s_29815_C2H2 = 200.93  # (J/(mol*K)) P = 1bar

Hf_29815_H2 = 0.0  # this is Hf (kJ/(mol*K))
s_29815_H2 = 130.68  # (J/(mol*K)) P = 1bar

R = 8.314*10**(-3) # kJ/(mol K)
P_0 = 1.0 #bar
P_1 = 1.0 #bar

'''
# SOLAR ELEMENTARY ABUNDANCES: (Lodders, 2003)
n_H = 1.0
n_He = 0.0968277856261
n_Ti = 9.54992586021*10**(-8)
n_V = 1.09647819614*10**(-8)
n_O = 0.000602559586074
n_C = 0.000275422870334
n_N = 8.12830516164*10**(-5)
n_S = 1.62181009736*10**(-5)
n_Na = 2.23872113857*10**(-6)
n_K = 1.44543977075*10**(-7)
n_Fe = 3.2359365693*10**(-5)
'''

# Using Looders Data: ratio_CO_solar = 0.46

# If you want to generate abundances with C / O = X, for example, simply increase the carbon by a factor of (X / [C / O_solar]), so the new ratio of carbon to oxygen will be X

# Kevin et al. used n_0 = 5*10**-4 and rate_CO_solar = 0.5 (hence n_C = 2.5 *10**-4)

n_O = 5.0 * 10**(-4)
n_C = 2.5 * 10**(-4)

rate_CO_solar = n_C/n_O  # This number is C/O = 0.46 for abundances of Lodders and C/O = 0.5 for Kevin et al.

def dH(A,B,C,D,E,F,H,t): # kJ/(mol*K)
    '''
    :param A, B, C, D, E , F, H:
    Shomate equation parameters for thermochemical functions. Data taken from the NIST database.
    :param t: Time. It's define as t=T/1000 (type = array)
    :return: d Enthalpy
    '''
    dH_MOLECULE = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H
    return dH_MOLECULE

def s(A,B,C,D,E,G,t): # J/(mol*K)
    '''
    :param A, B, C, D, E , F, H:
    Shomate equation parameters for thermochemical functions. Data taken from the NIST database.
    :param t: Time. It's define as t=T/1000 (type = array)
    :return: Entropy
    '''
    s_MOLECULE = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    return s_MOLECULE

def dG(dH,t,dS):
    '''
    :param dH: Enthalpy
    :param t: Time. It's define as t=T/1000 (type = array)
    :param dS: Entropy
    :return: Gibbs Fre Energy
    '''
    dG = dH - t*dS
    return dG
    
# enthalpy and free energy of reaction at 298.15 K for one reaction (e.g. CO2 + H2 <-> CO + H2O)
    
def Hrxn_29815(Hf_29815_reactant1,Hf_29815_reactant2,Hf_29815_product1,Hf_29815_product2):
    Hrxn_29815 = Hf_29815_product1 + Hf_29815_product2 - Hf_29815_reactant1 - Hf_29815_reactant2
    return Hrxn_29815
    
def Srxn_29815(S_29815_reactant1,S_29815_reactant2,S_29815_product1,S_29815_product2):
    Srxn_29815 = S_29815_product1 + S_29815_product2 - S_29815_reactant1 - S_29815_reactant2
    return Srxn_29815
    
def Grxn_29815(Hrxn_29815,Srxn_29815):
    Grxn_29815 = Hrxn_29815 - 298.15*(Srxn_29815)/1000
    return Grxn_29815

def molecule_H_S(Hf_29815_molecule, s_29815_molecule,A,B,C,D,E,F,G,H,t):
    dH_molecule = dH(A, B, C, D, E, F, H,t)
    h_molecule = dH_molecule + Hf_29815_molecule
    s_molecule = s(A, B, C, D, E, G,t)
    ds_molecule = s_molecule - s_29815_molecule
    return dH_molecule, h_molecule, s_molecule, ds_molecule

def insert_one_temp(temperature, P_wanted):

    file_T_H_S_molecules = open('T_H_S_molecules.dat', 'w')

    file_T_H_S_molecules.write("\n        Molecule  |  Temperature |       H       |       S       |")
    file_T_H_S_molecules.write("\n     ----------------------------------------------------------------------------------")

    t = temperature/1000

    #       H2O

    if temperature >= 500.0 and temperature <= 1700.0:

        # 500-1700 K valid temperature range
        A = 30.09200
        B = 6.832514
        C = 6.793435
        D = -2.534480
        E = 0.082139
        F = -250.8810
        G = 223.3967
        H = -241.8264

    elif temperature > 1700.0 and temperature <= 6000.0:

        # 1700-6000 K valid temperature range
        A = 41.96426
        B = 8.622053
        C = -1.499780
        D = 0.098119
        E = -11.15764
        F = -272.1797
        G = 219.7809
        H = -241.8264

    else:
        print("ERROR: Enter a valid temperature between 500K and 6000K!")

    dH_H2O, h_H2O, s_H2O, ds_H2O = molecule_H_S(Hf_29815_H2O, s_29815_H2O, A, B, C, D, E, F, G, H, t)

    file_T_H_S_molecules.write("\n \t" + "H2O" + str(temperature) + "  |" + str(h_H2O) + "|" + str(s_H2O) + "|")

    #       CH4

    if temperature >= 298.0 and temperature <= 1300.0:
        # 298-1300 K valid temperature range
        A = -0.73029
        B = 108.4773
        C = -42.52157
        D = 5.862788
        E = 0.678565
        F = -76.84376
        G = 158.7163
        H = -74.87310

    elif temperature > 1300.0 and temperature <= 6000.0:

        # 1300-6000 K valid temperature range
        A = 85.81217
        B = 11.26467
        C = -2.114146
        D = 0.138190
        E = -26.42221
        F = -153.5327
        G = 224.4143
        H = -74.87310

    else:
        print("ERROR: Enter a valid temperature between 298K and 6000K!")

    dH_CH4, h_CH4, s_CH4, ds_CH4 = molecule_H_S(Hf_29815_CH4, s_29815_CH4, A, B, C, D, E, F, G, H, t)

    file_T_H_S_molecules.write("\n \t" + "CH4" + str(temperature) + "  |" + str(h_CH4) + "|" + str(s_CH4) + "|")

    #       CO

    if temperature >= 298.0 and temperature <= 1300.0:
        # 298-1300 K valid temperature range
        A = 25.56759
        B = 6.096130
        C = 4.054656
        D = -2.671301
        E = 0.131021
        F = -118.0089
        G = 227.3665
        H = -110.5271

    elif temperature > 1300.0 and temperature <= 6000.0:

        # 1300-6000 K valid temperature range
        A = 35.15070
        B = 1.300095
        C = -0.205921
        D = 0.013550
        E = -3.282780
        F = -127.8375
        G = 231.7120
        H = -110.5271

    else:
        print("ERROR: Enter a valid temperature between 298K and 6000K!")

    dH_CO, h_CO, s_CO, ds_CO = molecule_H_S(Hf_29815_CO, s_29815_CO, A, B, C, D, E, F, G, H, t)

    file_T_H_S_molecules.write("\n \t" + "CO" + str(temperature) + "  |" + str(h_CO) + "|" + str(s_CO) + "|")

    #       CO2

    if temperature >= 298.0 and temperature <= 1200.0:
        # 298-1200 K valid temperature range
        A = 24.99735
        B = 55.18696
        C = -33.69137
        D = 7.948387
        E = -0.136638
        F = -403.6075
        G = 228.2431
        H = -393.5224

    elif temperature > 1200.0 and temperature <= 6000.0:

        # 1200-6000 K valid temperature range
        A = 58.16639
        B = 2.720074
        C = -0.492289
        D = 0.038844
        E = -6.447293
        F = -425.9186
        G = 263.6125
        H = -393.5224

    else:
        print("ERROR: Enter a valid temperature between 298K and 6000K!")

    dH_CO2, h_CO2, s_CO2, ds_CO2 = molecule_H_S(Hf_29815_CO2, s_29815_CO2, A, B, C, D, E, F, G, H, t)

    file_T_H_S_molecules.write("\n \t" + "CO2" + str(temperature) + "  |" + str(h_CO2) + "|" + str(s_CO2) + "|")


    #       C2H2

    if temperature >= 298.0 and temperature <= 1100.0:

        # 298-1100 K valid temperature range
        A = 40.68697
        B = 40.73279
        C = -16.17840
        D = 3.669741
        E = -0.658411
        F = 210.7067
        G = 235.0052
        H = 226.7314

    elif temperature > 1100.0 and temperature <= 6000.0:

        # 1100-6000 K valid temperature range
        A = 67.47244
        B = 11.75110
        C = -2.021470
        D = 0.136195
        E = -9.806418
        F = 185.4550
        G = 253.5337
        H = 226.7314

    else:
        print("ERROR: Enter a valid temperature between 298K and 6000K!")

    dH_C2H2, h_C2H2, s_C2H2, ds_C2H2 = molecule_H_S(Hf_29815_C2H2, s_29815_C2H2, A, B, C, D, E, F, G, H, t)

    file_T_H_S_molecules.write("\n \t" + "C2H2" + str(temperature) + "  |" + str(h_C2H2) + "|" + str(s_C2H2) + "|")

    #       H2

    if temperature >= 298.0 and temperature <= 1000.0:
        # 298-1000 K valid temperature range
        A = 33.066178
        B = -11.363417
        C = 11.432816
        D = -2.772874
        E = -0.158558
        F = -9.980797
        G = 172.707974
        H = 0.0

    elif temperature > 1000.0 and temperature <= 2500.0:
        # 1000-2500 K valid temperature range
        A = 18.563083
        B = 12.257357
        C = -2.859786
        D = 0.268238
        E = 1.977990
        F = -1.147438
        G = 156.288133
        H = 0.0

    elif temperature > 2500.0 and temperature <= 6000.0:

        # 2500-6000 K valid temperature range
        A = 43.413560
        B = -4.293079
        C = 1.272428
        D = -0.096876
        E = -20.533862
        F = -38.515158
        G = 162.081354
        H = 0.0

    else:
        print("ERROR: Enter a valid temperature between 298K and 6000K!")

    dH_H2, h_H2, s_H2, ds_H2 = molecule_H_S(Hf_29815_H2, s_29815_H2, A, B, C, D, E, F, G, H, t)

    file_T_H_S_molecules.write("\n \t" + "H2" + str(temperature) + "  |" + str(h_H2) + "|" + str(s_H2) + "|")

    file_T_H_S_molecules.close()

    # reaction 1: CH4 + H2O <-> CO + 3 H2

    print("\n# Reaction 1: CH4 + H2O <-> CO + 3 H2")

    Hrxn1_29815 = Hrxn_29815(Hf_29815_CH4, Hf_29815_H2O, Hf_29815_CO, 3 * Hf_29815_H2)
    Srxn1_29815 = Srxn_29815(s_29815_CH4, s_29815_H2O, s_29815_CO, 3 * s_29815_H2)
    Grxn1_29815 = Grxn_29815(Hrxn1_29815, Srxn1_29815)

    Hrxn1 = Hrxn1_29815 + dH_CO + 3 * dH_H2 - dH_CH4 - dH_H2O
    Grxn1 = Hrxn1 - temperature * (s_CO + 3 * s_H2 - s_CH4 - s_H2O) / 1000

    # Equilibrium constant calculation (K')

    K1_P_wanted = (P_0 / P_wanted) ** 2 * np.exp(-Grxn1 / (R * temperature))

    print("Ready :D! \n")

    # reaction 2: CO2 + H2 <-> CO + H2O

    print("#reaction 2: CO2 + H2 <-> CO + H2O")

    Hrxn2_29815 = Hrxn_29815(Hf_29815_CO2, Hf_29815_H2, Hf_29815_CO, Hf_29815_H2O)
    Srxn2_29815 = Srxn_29815(s_29815_CO2, s_29815_H2, s_29815_CO, s_29815_H2O)
    Grxn2_29815 = Grxn_29815(Hrxn2_29815, Srxn2_29815)

    Hrxn2 = Hrxn2_29815 + dH_CO + dH_H2O - dH_CO2 - dH_H2
    Grxn2 = Hrxn2 - temperature * (s_CO + s_H2O - s_CO2 - s_H2) / 1000


    # Equilibrium constant calculation (K')

    K2_sin_presion = np.exp(-Grxn2 / (R * temperature))

    print("Ready :D! \n")

    # reaction 3: 2 CH4 <-> C2H2 + 3 H2

    print("#reaction 3: 2 CH4 <-> C2H2 + 3 H2")

    Hrxn3_29815 = Hrxn_29815(2 * Hf_29815_CH4, 0.0, Hf_29815_C2H2, 3 * Hf_29815_H2)
    Srxn3_29815 = Srxn_29815(2 * s_29815_CH4, 0.0, s_29815_C2H2, 3 * s_29815_H2)
    Grxn3_29815 = Grxn_29815(Hrxn3_29815, Srxn3_29815)

    Hrxn3 = Hrxn3_29815 + dH_C2H2 + 3 * dH_H2 - 2 * dH_CH4
    Grxn3 = Hrxn3 - temperature * (s_C2H2 + 3 * s_H2 - 2 * s_CH4) / 1000

    # Equilibrium constant calculation (K')

    K3_P_wanted = (P_0 / P_wanted) ** 2 * np.exp(-Grxn3 / (R * temperature))

    print("Ready :D! \n")

    return K1_P_wanted, K2_sin_presion, K3_P_wanted

def classic_calculation():

    T = np.linspace(500, 1700, 13)  # degrees K
    t = T / 1000

    #       H2O

    # 500-1700 K valid temperature range
    A = 30.09200
    B = 6.832514
    C = 6.793435
    D = -2.534480
    E = 0.082139
    F = -250.8810
    G = 223.3967
    H = -241.8264

    dH_H2O_1 = dH(A, B, C, D, E, F, H,t)
    h_H2O_1 = dH_H2O_1 + Hf_29815_H2O
    s_H2O_1 = s(A, B, C, D, E, G,t)
    ds_H2O_1 = s_H2O_1 - s_29815_H2O

    T = np.linspace(1800, 3000, 13)  # degrees K
    t = T / 1000

    # 1700-6000 K valid temperature range
    A = 41.96426
    B = 8.622053
    C = -1.499780
    D = 0.098119
    E = -11.15764
    F = -272.1797
    G = 219.7809
    H = -241.8264

    dH_H2O_2 = dH(A, B, C, D, E, F, H,t)
    h_H2O_2 = dH_H2O_2 + Hf_29815_H2O
    s_H2O_2 = s(A, B, C, D, E, G,t)
    ds_H2O_2 = s_H2O_2 - s_29815_H2O

    dH_H2O = np.array(list(dH_H2O_1) + list(dH_H2O_2))
    h_H2O = np.array(list(h_H2O_1) + list(h_H2O_2))
    s_H2O = np.array(list(s_H2O_1) + list(s_H2O_2))
    ds_H2O = np.array(list(ds_H2O_1) + list(ds_H2O_2))

    file_T_H_S_H2O = open('T_H_S_H2O.dat', 'w')

    file_T_H_S_H2O.write("\n                                            WATER (H2O)")
    file_T_H_S_H2O.write("\n      Temperature |      dH       |       H       |       dS      |       S        ")
    file_T_H_S_H2O.write(
        "\n     ---------------------------------------------------------------------------------------------------------------")
    for i in np.arange(26):
        file_T_H_S_H2O.write(
            "\n \t" + str(np.linspace(500, 3000, 26)[i]) + "  |" + str(dH_H2O[i]) + "|" + str(h_H2O[i]) + "|" + str(
                ds_H2O[i]) + "|" + str(s_H2O[i]) + "|")

    file_T_H_S_H2O.close()

    #       CH4

    T = np.linspace(500, 1300, 9)  # degrees K
    t = T / 1000

    # 298-1300 K valid temperature range
    A = -0.73029
    B = 108.4773
    C = -42.52157
    D = 5.862788
    E = 0.678565
    F = -76.84376
    G = 158.7163
    H = -74.87310

    dH_CH4_1 = dH(A, B, C, D, E, F, H,t)
    h_CH4_1 = dH_CH4_1 + Hf_29815_CH4
    s_CH4_1 = s(A, B, C, D, E, G,t)
    ds_CH4_1 = s_CH4_1 - s_29815_CH4

    T = np.linspace(1400, 3000, 17)  # degrees K
    t = T / 1000

    # 1300-6000 K valid temperature range
    A = 85.81217
    B = 11.26467
    C = -2.114146
    D = 0.138190
    E = -26.42221
    F = -153.5327
    G = 224.4143
    H = -74.87310

    dH_CH4_2 = dH(A, B, C, D, E, F, H,t)
    h_CH4_2 = dH_CH4_2 + Hf_29815_CH4
    s_CH4_2 = s(A, B, C, D, E, G,t)
    ds_CH4_2 = s_CH4_2 - s_29815_CH4

    dH_CH4 = np.array(list(dH_CH4_1) + list(dH_CH4_2))
    h_CH4 = np.array(list(h_CH4_1) + list(h_CH4_2))
    s_CH4 = np.array(list(s_CH4_1) + list(s_CH4_2))
    ds_CH4 = np.array(list(ds_CH4_1) + list(ds_CH4_2))

    file_T_H_S_CH4 = open('T_H_S_CH4.dat', 'w')

    file_T_H_S_CH4.write("\n                                            METHANE (CH4)")
    file_T_H_S_CH4.write("\n      Temperature |      dH       |       H       |       dS      |       S        ")
    file_T_H_S_CH4.write("\n     ----------------------------------------------------------------------------------")

    for i in np.arange(26):
        file_T_H_S_CH4.write(
            "\n \t" + str(np.linspace(500, 3000, 26)[i]) + "  |" + str(dH_CH4[i]) + "|" + str(h_CH4[i]) + "|" + str(
                ds_CH4[i]) + "|" + str(s_CH4[i]) + "|")

    file_T_H_S_CH4.close()

    #       CO

    T = np.linspace(500, 1300, 9)  # degrees K
    t = T / 1000

    # 298-1300 K valid temperature range
    A = 25.56759
    B = 6.096130
    C = 4.054656
    D = -2.671301
    E = 0.131021
    F = -118.0089
    G = 227.3665
    H = -110.5271


    dH_CO_1 = dH(A, B, C, D, E, F, H,t)
    h_CO_1 = dH_CO_1 + Hf_29815_CO
    s_CO_1 = s(A, B, C, D, E, G,t)
    ds_CO_1 = s_CO_1 - s_29815_CO

    T = np.linspace(1400, 3000, 17)  # degrees K
    t = T / 1000

    # 1300-6000 K valid temperature range
    A = 35.15070
    B = 1.300095
    C = -0.205921
    D = 0.013550
    E = -3.282780
    F = -127.8375
    G = 231.7120
    H = -110.5271

    dH_CO_2 = dH(A, B, C, D, E, F, H, t)
    h_CO_2 = dH_CO_2 + Hf_29815_CO
    s_CO_2 = s(A, B, C, D, E, G, t)
    ds_CO_2 = s_CO_2 - s_29815_CO

    dH_CO = np.array(list(dH_CO_1) + list(dH_CO_2))
    h_CO = np.array(list(h_CO_1) + list(h_CO_2))
    s_CO = np.array(list(s_CO_1) + list(s_CO_2))
    ds_CO = np.array(list(ds_CO_1) + list(ds_CO_2))

    file_T_H_S_CO = open('T_H_S_CO.dat', 'w')

    file_T_H_S_CO.write("\n                                           CARBON MONOXIDE (CO)")
    file_T_H_S_CO.write("\n      Temperature |      dH       |       H       |       dS      |       S        ")
    file_T_H_S_CO.write("\n     ---------------------------------------------------------------------------------")

    for i in np.arange(26):
        file_T_H_S_CO.write(
            "\n \t" + str(np.linspace(500, 3000, 26)[i]) + "  |" + str(dH_CO[i]) + "|" + str(h_CO[i]) + "|" + str(
                ds_CO[i]) + "|" + str(s_CO[i]) + "|")

    file_T_H_S_CO.close()

    #       H2


    T = np.linspace(500, 1000, 6)  # degrees K
    t = T / 1000

    # 298-1000 K valid temperature range
    A = 33.066178
    B = -11.363417
    C = 11.432816
    D = -2.772874
    E = -0.158558
    F = -9.980797
    G = 172.707974
    H = 0.0


    dH_H2_1 = dH(A, B, C, D, E, F, H, t)
    h_H2_1 = dH_H2_1 + Hf_29815_H2
    s_H2_1 = s(A, B, C, D, E, G, t)
    ds_H2_1 = s_H2_1 - s_29815_H2

    T = np.linspace(1100, 2500, 15)  # degrees K
    t = T / 1000

    # 1000-2500 K valid temperature range
    A = 18.563083
    B = 12.257357
    C = -2.859786
    D = 0.268238
    E = 1.977990
    F = -1.147438
    G = 156.288133
    H = 0.0

    dH_H2_2 = dH(A, B, C, D, E, F, H, t)
    h_H2_2 = dH_H2_2 + Hf_29815_H2
    s_H2_2 = s(A, B, C, D, E, G, t)
    ds_H2_2 = s_H2_2 - s_29815_H2

    T = np.linspace(2600, 3000, 5)  # degrees K
    t = T / 1000

    # 2500-6000 K valid temperature range
    A = 43.413560
    B = -4.293079
    C = 1.272428
    D = -0.096876
    E = -20.533862
    F = -38.515158
    G = 162.081354
    H = 0.0

    dH_H2_3 = dH(A, B, C, D, E, F, H, t)
    h_H2_3 = dH_H2_3 + Hf_29815_H2
    s_H2_3 = s(A, B, C, D, E, G, t)
    ds_H2_3 = s_H2_3 - s_29815_H2

    dH_H2 = np.array(list(dH_H2_1) + list(dH_H2_2) + list(dH_H2_3))
    h_H2 = np.array(list(h_H2_1) + list(h_H2_2) + list(h_H2_3))
    s_H2 = np.array(list(s_H2_1) + list(s_H2_2) + list(s_H2_3))
    ds_H2 = np.array(list(ds_H2_1) + list(ds_H2_2) + list(ds_H2_3))

    file_T_H_S_H2 = open('T_H_S_H2.dat', 'w')

    file_T_H_S_H2.write("\n                                            HYDROGEN (H2)")
    file_T_H_S_H2.write("\n      Temperature |      dH       |       H       |       dS      |       S        ")
    file_T_H_S_H2.write("\n     -------------------------------------------------------------------------------------")
    for i in np.arange(26):
        file_T_H_S_H2.write(
            "\n \t" + str(np.linspace(500, 3000, 26)[i]) + "  |" + str(dH_H2[i]) + "|" + str(h_H2[i]) + "|" + str(
                ds_H2[i]) + "|" + str(s_H2[i]) + "|")

    file_T_H_S_H2.close()

    #       CO2

    T = np.linspace(500, 1200, 8)  # degrees K
    t = T / 1000

    # 298-1200 K valid temperature range
    A = 24.99735
    B = 55.18696
    C = -33.69137
    D = 7.948387
    E = -0.136638
    F = -403.6075
    G = 228.2431
    H = -393.5224


    dH_CO2_1 = dH(A, B, C, D, E, F, H,t)
    h_CO2_1 = dH_CO2_1 + Hf_29815_CO2
    s_CO2_1 = s(A, B, C, D, E, G,t)
    ds_CO2_1 = s_CO2_1 - s_29815_CO2

    T = np.linspace(1300, 3000, 18)  # degrees K
    t = T / 1000

    # 1200-6000 K valid temperature range
    A = 58.16639
    B = 2.720074
    C = -0.492289
    D = 0.038844
    E = -6.447293
    F = -425.9186
    G = 263.6125
    H = -393.5224

    dH_CO2_2 = dH(A, B, C, D, E, F, H,t)
    h_CO2_2 = dH_CO2_2 + Hf_29815_CO2
    s_CO2_2 = s(A, B, C, D, E, G,t)
    ds_CO2_2 = s_CO2_2 - s_29815_CO2

    dH_CO2 = np.array(list(dH_CO2_1) + list(dH_CO2_2))
    h_CO2 = np.array(list(h_CO2_1) + list(h_CO2_2))
    s_CO2 = np.array(list(s_CO2_1) + list(s_CO2_2))
    ds_CO2 = np.array(list(ds_CO2_1) + list(ds_CO2_2))

    file_T_H_S_CO2 = open('T_H_S_CO2.dat', 'w')

    file_T_H_S_CO2.write("\n                                           CARBON DIOXIDE (CO2)")
    file_T_H_S_CO2.write("\n      Temperature |      dH       |       H       |       dS      |       S        ")
    file_T_H_S_CO2.write("\n     --------------------------------------------------------------------------------")

    for i in np.arange(26):
        file_T_H_S_CO2.write(
            "\n \t" + str(np.linspace(500, 3000, 26)[i]) + "  |" + str(dH_CO2[i]) + "|" + str(h_CO2[i]) + "|" + str(
                ds_CO2[i]) + "|" + str(s_CO2[i]) + "|")

    file_T_H_S_CO2.close()

    #       C2H2

    T = np.linspace(500, 1100, 7)  # degrees K
    t = T / 1000

    # 298-1100 K valid temperature range
    A = 40.68697
    B = 40.73279
    C = -16.17840
    D = 3.669741
    E = -0.658411
    F = 210.7067
    G = 235.0052
    H = 226.7314


    dH_C2H2_1 = dH(A, B, C, D, E, F, H,t)
    h_C2H2_1 = dH_C2H2_1 + Hf_29815_C2H2
    s_C2H2_1 = s(A, B, C, D, E, G,t)
    ds_C2H2_1 = s_C2H2_1 - s_29815_C2H2

    T = np.linspace(1200, 3000, 19)  # degrees K
    t = T / 1000

    # 1100-6000 K valid temperature range
    A = 67.47244
    B = 11.75110
    C = -2.021470
    D = 0.136195
    E = -9.806418
    F = 185.4550
    G = 253.5337
    H = 226.7314

    dH_C2H2_2 = dH(A, B, C, D, E, F, H,t)
    h_C2H2_2 = dH_C2H2_2 + Hf_29815_C2H2
    s_C2H2_2 = s(A, B, C, D, E, G,t)
    ds_C2H2_2 = s_C2H2_2 - s_29815_C2H2

    dH_C2H2 = np.array(list(dH_C2H2_1) + list(dH_C2H2_2))
    h_C2H2 = np.array(list(h_C2H2_1) + list(h_C2H2_2))
    s_C2H2 = np.array(list(s_C2H2_1) + list(s_C2H2_2))
    ds_C2H2 = np.array(list(ds_C2H2_1) + list(ds_C2H2_2))

    file_T_H_S_C2H2 = open('T_H_S_C2H2.dat', 'w')

    file_T_H_S_C2H2.write("\n                                           ACETYLENE (C2H2)")
    file_T_H_S_C2H2.write("\n      Temperature |      dH       |       H       |       dS      |       S        ")
    file_T_H_S_C2H2.write("\n     ----------------------------------------------------------------------------------")

    for i in np.arange(26):
        file_T_H_S_C2H2.write(
            "\n \t" + str(np.linspace(500, 3000, 26)[i]) + "  |" + str(dH_C2H2[i]) + "|" + str(h_C2H2[i]) + "|" + str(
                ds_C2H2[i]) + "|" + str(s_C2H2[i]) + "|")

    file_T_H_S_C2H2.close()

    T = np.linspace(500, 3000, 26)

    # reaction 1: CH4 + H2O <-> CO + 3 H2

    print("# Reaction 1: CH4 + H2O <-> CO + 3 H2")

    Hrxn1_29815 = Hrxn_29815(Hf_29815_CH4, Hf_29815_H2O, Hf_29815_CO, 3 * Hf_29815_H2)
    Srxn1_29815 = Srxn_29815(s_29815_CH4, s_29815_H2O, s_29815_CO, 3 * s_29815_H2)
    Grxn1_29815 = Grxn_29815(Hrxn1_29815, Srxn1_29815)

    Hrxn1 = Hrxn1_29815 + dH_CO + 3 * dH_H2 - dH_CH4 - dH_H2O
    Grxn1 = Hrxn1 - T * (s_CO + 3 * s_H2 - s_CH4 - s_H2O) / 1000

    '''

    #print(Grxn1)

    plt.figure()
    plt.plot(T,Grxn1, label='$\Delta G_{rxn1}$')
    plt.plot(T,Hrxn1, label='$\Delta H_{rxn1}$')
    plt.xlabel('Temperature (K)')
    plt.ylabel('(kJ/mol)')
    plt.legend( loc='best')
    plt.savefig("abundances-nist-1.png")
    #plt.show()

    '''

    # Equilibrium constant calculation (K')

    K1_1bar = (P_0 / P_1) ** 2 * np.exp(-Grxn1 / (R * T))
    K1_P_wanted = (P_0 / P_wanted) ** 2 * np.exp(-Grxn1 / (R * T))

    print("Ready :D! \n")

    '''

    plt.figure()
    plt.plot(T,K1_1bar,label='$K1_{1bar}$')
    plt.plot(T,K1_P_wanted,label='$K1_{P_wanted}$')
    plt.xlim([500, 3000])
    plt.yscale('log')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Equilibrium constant rxn1')
    plt.legend( loc='best')
    plt.savefig('abundances-nist-1-K.png')
    #plt.show()

    '''

    # reaction 2: CO2 + H2 <-> CO + H2O

    print("#reaction 2: CO2 + H2 <-> CO + H2O")

    Hrxn2_29815 = Hrxn_29815(Hf_29815_CO2, Hf_29815_H2, Hf_29815_CO, Hf_29815_H2O)
    Srxn2_29815 = Srxn_29815(s_29815_CO2, s_29815_H2, s_29815_CO, s_29815_H2O)
    Grxn2_29815 = Grxn_29815(Hrxn2_29815, Srxn2_29815)

    Hrxn2 = Hrxn2_29815 + dH_CO + dH_H2O - dH_CO2 - dH_H2
    Grxn2 = Hrxn2 - T * (s_CO + s_H2O - s_CO2 - s_H2) / 1000

    '''

    #print(Grxn2)

    plt.figure()
    plt.plot(T,Grxn2, label='$\Delta G_{rxn2}$')
    plt.plot(T,Hrxn2, label='$\Delta H_{rxn2}$')
    plt.xlabel('Temperature (K)')
    plt.ylabel('(kJ/mol)')
    plt.legend( loc='best')
    plt.savefig("abundances-nist-2.png")
    #plt.show()

    '''

    # Equilibrium constant calculation (K')

    K2_sin_presion = np.exp(-Grxn2 / (R * T))

    print("Ready :D! \n")

    '''

    plt.figure()
    plt.plot(T,K2_1bar,label='$K´_{2}(1bar)$')
    plt.plot(T,K2_P_wanted,label='$K´_{2}(P_wanted)$',ls='--')
    plt.plot(T,K2_sin_presion,label='$K´_{2}$',ls='--')
    plt.xlim([500, 3000])
    plt.yscale('log')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Equilibrium constant rxn2')
    plt.legend( loc='best')
    plt.savefig('abundances-nist-2-K.png')
    #plt.show()

    plt.figure()
    plt.plot(T,1/K2_sin_presion,label='$1/K´_{2}$',color='r')
    plt.xlim([500, 3000])
    plt.yscale('log')
    plt.xlabel('Temperature (K)')
    plt.ylabel('$1/K´_{2}$')
    plt.legend( loc='best')
    plt.savefig('abundances-nist-fig1.png')
    #plt.show()

    '''

    # reaction 3: 2 CH4 <-> C2H2 + 3 H2

    print("#reaction 3: 2 CH4 <-> C2H2 + 3 H2")

    Hrxn3_29815 = Hrxn_29815(2 * Hf_29815_CH4, 0.0, Hf_29815_C2H2, 3 * Hf_29815_H2)
    Srxn3_29815 = Srxn_29815(2 * s_29815_CH4, 0.0, s_29815_C2H2, 3 * s_29815_H2)
    Grxn3_29815 = Grxn_29815(Hrxn3_29815, Srxn3_29815)

    Hrxn3 = Hrxn3_29815 + dH_C2H2 + 3 * dH_H2 - 2 * dH_CH4
    Grxn3 = Hrxn3 - T * (s_C2H2 + 3 * s_H2 - 2 * s_CH4) / 1000

    '''

    #print(Grxn3)

    plt.figure()
    plt.plot(T,Grxn3, label='$\Delta G_{rxn3}$')
    plt.plot(T,Hrxn3, label='$\Delta H_{rxn3}$')
    plt.xlabel('Temperature (K)')
    plt.ylabel('(kJ/mol)')
    plt.legend( loc='best')
    plt.savefig("abundances-nist-3.png")
    #plt.show()

    '''

    # Equilibrium constant calculation (K')

    K3_1bar = (P_0 / P_1) ** 2 * np.exp(-Grxn3 / (R * T))
    K3_P_wanted = (P_0 / P_wanted) ** 2 * np.exp(-Grxn3 / (R * T))

    print("Ready :D! \n")

    '''

    plt.figure()
    plt.plot(T,K3_1bar,label='$K3_{1bar}$')
    plt.plot(T,K3_P_wanted,label='$K3_{P_wanted}$')
    plt.xlim([500, 3000])
    plt.yscale('log')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Equilibrium constant rxn3')
    plt.legend( loc='best')
    plt.savefig('abundances-nist-3-K.png')
    #plt.show()

    '''

    '''

    # This generates a graphic of the constants (K1, K2 and K3) with different pressures

    plt.figure()
    plt.plot(T,K1_1bar,label='$K´_{1}(1bar)$',color='orange',ls='-',linewidth=1)
    plt.plot(T,K1_P_wanted,label='$K´_{1}('+P_wanted_str+'bar)$',color='orange',ls='-',linewidth=3)
    plt.plot(T,K2_sin_presion,label='$K´_{2}$',color='r',ls='-.')
    plt.plot(T,K3_1bar,label='$K´_{3}(1bar)$',color='k',ls='--',linewidth=1)
    plt.plot(T,K3_P_wanted,label='$K´_{3}('+P_wanted_str+'bar)$',color='k',ls='--',linewidth=3)
    plt.xlim([500, 3000])
    plt.ylim([10**-22, 10**14])
    plt.yscale('log')
    plt.xticks(np.linspace(500,3000,6,endpoint=True))
    plt.yticks(10**(np.linspace(-22,14,13)))
    plt.tick_params(direction='in',bottom='on', top='on', left='on', right='on')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Normalised Equilibrium Constants')
    plt.legend( loc='best')
    plt.savefig('abundances-nist_P_'+P_wanted_str+'.png')
    plt.show()
    '''
    return K1_1bar, K1_P_wanted, K2_sin_presion, K3_1bar, K3_P_wanted


def abundances_norm_H(CO_wanted, P_wanted_str, selection):
    n_C_wanted = (CO_wanted / rate_CO_solar) * n_C

    if selection == '1':
        # METHANE

        T = np.linspace(500, 3000, 26)
        count_1bar = 0
        count_P_wanted = 0

        file_roots_calculation = open('roots_calculation.dat', 'w')

        for j in np.arange(2):
            A = []
            if j == 0:
                print("P = 1 bar")
                K1 = K1_1bar
                K2 = K2_sin_presion
                K3 = K3_1bar
                n_CH4_1bar = []
            if j == 1:
                print('P=' + P_wanted_str)
                K1 = K1_P_wanted
                K2 = K2_sin_presion
                K3 = K3_P_wanted
                n_CH4_P_wanted = []
            A.append(-2.0 * n_C_wanted)
            A.append(8.0 * (K1 / K2) * (n_O - n_C_wanted) ** 2.0 + 1.0 + 2.0 * K1 * (n_O - n_C_wanted))
            A.append(8.0 * (K1 / K2) * (n_O - n_C_wanted) + 2.0 * K3 + K1)
            A.append(2.0 * (K1 / K2) * (1.0 + 8.0 * K3 * (n_O - n_C_wanted)) + 2.0 * K1 * K3)
            A.append(8.0 * (K1 * K3 / K2))
            A.append(8.0 * (K1 * (K3 ** 2.0) / K2))
            for i in np.arange(26):
                file_roots_calculation.write("# Coefficients (A)")
                file_roots_calculation.write(
                    '\n' + str(A[0]) + '|' + str(A[1][i]) + '|' + str(A[2][i]) + '|' + str(A[3][i]) + '|' + str(
                        A[4][i]) + '|' + str(A[5][i]))
                roots = poly.polyroots([A[0], A[1][i], A[2][i], A[3][i], A[4][i], A[5][i]])
                file_roots_calculation.write("\n# Roots for T=" + str(T[i]))
                for w in np.arange(len(roots)):
                    two_n_C = float(-A[0])  # We know that n_CH4 should be approx. 2*n_C
                    possible_n_CH4 = roots[w]
                    file_roots_calculation.write("\nPossible ñ_CH4 (root " + str(w + 1) + "): " + str(possible_n_CH4))
                    file_roots_calculation.write("\n2 ñ_C (We know that ñ_CH4 should be approx. this): " + str(two_n_C))
                    if ('{0:.1f}'.format(possible_n_CH4.real) == '{0:.1f}'.format(two_n_C) or (
                                    possible_n_CH4.real <= two_n_C and possible_n_CH4.real > 0.0)) and possible_n_CH4.imag == 0.0:
                        file_roots_calculation.write(
                            "\n~~~~~ MATCH 1313 ~~~~~(T=" + str(T[i]) + ") | Possible_n_CH4.real=" + str(
                                possible_n_CH4.real) + " | Two_n_C=" + '{0:.6f}'.format(two_n_C) + '\n')
                        if j == 0:
                            n_CH4_1bar.append(possible_n_CH4.real)
                            count_1bar += 1  # This number should be 26 at the end of the loops (because we have 26 temperatures)
                        if j == 1:
                            n_CH4_P_wanted.append(possible_n_CH4.real)
                            count_P_wanted += 1

        print("***count_1bar*** = ", count_1bar)
        file_roots_calculation.write("\n***count_1bar*** = " + str(count_1bar))
        print("***count_P_wanted*** = ", count_P_wanted)
        file_roots_calculation.write("\n***count_P_wanted*** = " + str(count_P_wanted))
        print("(These two numbers should be 26 at the end of the loops (because we have 26 temperatures))")
        file_roots_calculation.write(
            "\n(These two numbers should be 26 at the end of the loops (because we have 26 temperatures))\n")

        file_roots_calculation.close()

        # METHANE
        file_abundance_CH4 = open('abundance_CH4_CO=' + str(int(CO_wanted)) + '.dat', "w")

        file_abundance_CH4.write("# METHANE")

        n_CH4_1bar = np.array(n_CH4_1bar)
        n_CH4_P_wanted = np.array(n_CH4_P_wanted)

        file_abundance_CH4.write(str(n_CH4_1bar) + "|" + str(len(n_CH4_1bar)) + "| P = 1bar")
        file_abundance_CH4.write(str(n_CH4_P_wanted) + "|" + str(len(n_CH4_P_wanted)) + "| P = " + P_wanted_str)

        file_abundance_CH4.close()

        # ACETYLENE
        file_abundance_C2H2 = open('abundance_C2H2_CO=' + str(int(CO_wanted)) + '.dat', "w")

        file_abundance_C2H2.write("# ACETYLENE")

        n_C2H2_1bar = K3_1bar * n_CH4_1bar ** 2.0
        n_C2H2_P_wanted = K3_P_wanted * n_CH4_P_wanted ** 2.0

        file_abundance_C2H2.write(str(n_C2H2_1bar) + "|" + str(len(n_C2H2_1bar)) + "| P = 1bar")
        file_abundance_C2H2.write(str(n_C2H2_P_wanted) + "|" + str(len(n_C2H2_P_wanted)) + "| P = " + P_wanted_str)

        file_abundance_C2H2.close()

        # WATER
        file_abundance_H2O = open('abundance_H2O_CO=' + str(int(CO_wanted)) + '.dat', "w")

        file_abundance_H2O.write("# WATER")

        n_H2O_1bar = (2.0 * K3_1bar * n_CH4_1bar ** 2.0) + n_CH4_1bar + (2.0 * (n_O - n_C_wanted))
        n_H2O_P_wanted = (2.0 * K3_P_wanted * n_CH4_P_wanted ** 2.0) + n_CH4_P_wanted + (2.0 * (n_O - n_C_wanted))

        file_abundance_H2O.write(str(n_H2O_1bar) + "|" + str(len(n_H2O_1bar)) + "| P = 1bar")
        file_abundance_H2O.write(str(n_H2O_P_wanted) + "|" + str(len(n_H2O_P_wanted)) + "| P = " + P_wanted_str)

        file_abundance_H2O.close()

        # CARBON MONOXIDE
        file_abundance_CO = open('abundance_CO_CO=' + str(int(CO_wanted)) + '.dat', "w")

        file_abundance_CO.write("# CARBON MONOXIDE")

        n_CO_1bar = K1_1bar * n_CH4_1bar * n_H2O_1bar
        n_CO_P_wanted = K1_P_wanted * n_CH4_P_wanted * n_H2O_P_wanted

        file_abundance_CO.write(str(n_CO_1bar) + "|" + str(len(n_CO_1bar)) + "| P = 1bar")
        file_abundance_CO.write(str(n_CO_P_wanted) + "|" + str(len(n_CO_P_wanted)) + "| P = " + P_wanted_str)

        file_abundance_CO.close()

        # CARBON DIOXIDE
        file_abundance_CO2 = open('abundance_CO2_CO=' + str(int(CO_wanted)) + '.dat', "w")

        file_abundance_CO2.write("# CARBON DIOXIDE")

        n_CO2_1bar = n_CO_1bar * n_H2O_1bar / K2_sin_presion
        n_CO2_P_wanted = n_CO_P_wanted * n_H2O_P_wanted / K2_sin_presion

        file_abundance_CO2.write(str(n_CO2_1bar) + "|" + str(len(n_CO2_1bar)) + "| P = 1bar")
        file_abundance_CO2.write(str(n_CO2_P_wanted) + "|" + str(len(n_CO2_P_wanted)) + "| P = " + P_wanted_str)

        file_abundance_CO2.close()

        return n_CH4_1bar, n_CH4_P_wanted, n_C2H2_1bar, n_C2H2_P_wanted, n_H2O_1bar, n_H2O_P_wanted, n_CO_1bar, n_CO_P_wanted, n_CO2_1bar, n_CO2_P_wanted

    if selection == '2':

        file_roots_calculation_selection2 = open('roots_calculation_selection2.dat', 'w')

        # METHANE

        K1 = K1_P_wanted
        K2 = K2_sin_presion
        K3 = K3_P_wanted

        print("K1=", K1)
        print("K2=", K2)
        print("K3=", K3)

        A0 = -2.0 * n_C_wanted
        A1 = 8.0 * (K1 / K2) * (n_O - n_C_wanted) ** 2.0 + 1.0 + 2.0 * K1 * (n_O - n_C_wanted)
        A2 = 8.0 * (K1 / K2) * (n_O - n_C_wanted) + 2.0 * K3 + K1
        A3 = 2.0 * (K1 / K2) * (1.0 + 8.0 * K3 * (n_O - n_C_wanted)) + 2.0 * K1 * K3
        A4 = 8.0 * (K1 * K3 / K2)
        A5 = 8.0 * (K1 * (K3 ** 2.0) / K2)

        file_roots_calculation_selection2.write("# Coefficients (A)")
        file_roots_calculation_selection2.write(
            '\n' + str(A0) + '|' + str(A1) + '|' + str(A2) + '|' + str(A3) + '|' + str(A4) + '|' + str(A5))
        roots = poly.polyroots([A0, A1, A2, A3, A4, A5])
        file_roots_calculation_selection2.write("\n# Roots for T=" + T_wanted_str)

        for w in np.arange(len(roots)):
            two_n_C = float(- A0)  # We know that n_CH4 should be approx. 2*n_C
            possible_n_CH4 = roots[w]
            file_roots_calculation_selection2.write("\nPossible ñ_CH4 (root " + str(w + 1) + "): " + str(possible_n_CH4))
            file_roots_calculation_selection2.write("\n2 ñ_C (We know that ñ_CH4 should be approx. this): " + str(two_n_C))
            if ('{0:.1f}'.format(possible_n_CH4.real) == '{0:.1f}'.format(two_n_C) or (
                    possible_n_CH4.real <= two_n_C and possible_n_CH4.real > 0.0)) and possible_n_CH4.imag == 0.0:
                file_roots_calculation_selection2.write(
                    "\n~~~~~ MATCH 1313 ~~~~~(T=" + T_wanted_str + ") | Possible_n_CH4.real=" + str(
                        possible_n_CH4.real) + " | Two_n_C=" + '{0:.6f}'.format(two_n_C) + '\n')
                n_CH4_P_wanted = possible_n_CH4.real

        file_roots_calculation_selection2.close()

        file_abundances_final_selection2 = open(
            'abundances_final_selection2_CO=' + str(int(CO_wanted)) + '_P_' + P_wanted_str + '_T_' + T_wanted_str + '.dat', "w")

        # METHANE

        file_abundances_final_selection2.write("# METHANE \n")
        file_abundances_final_selection2.write(
            str(n_CH4_P_wanted) + "| T = " + T_wanted_str + "| P = " + P_wanted_str + " | CO = " + CO_wanted_str + "\n")

        # ACETYLENE

        file_abundances_final_selection2.write("# ACETYLENE \n")

        n_C2H2_P_wanted = K3_P_wanted * n_CH4_P_wanted ** 2.0

        file_abundances_final_selection2.write(
            str(n_C2H2_P_wanted) + "| T = " + T_wanted_str + "| P = " + P_wanted_str + " | CO = " + CO_wanted_str + "\n")

        # WATER

        file_abundances_final_selection2.write("# WATER \n")

        n_H2O_P_wanted = (2.0 * K3_P_wanted * n_CH4_P_wanted ** 2.0) + n_CH4_P_wanted + (2.0 * (n_O - n_C_wanted))

        file_abundances_final_selection2.write(
            str(n_H2O_P_wanted) + "| T = " + T_wanted_str + "| P = " + P_wanted_str + " | CO = " + CO_wanted_str + "\n")

        # CARBON MONOXIDE

        file_abundances_final_selection2.write("# CARBON MONOXIDE \n")

        n_CO_P_wanted = K1_P_wanted * n_CH4_P_wanted * n_H2O_P_wanted

        file_abundances_final_selection2.write(
            str(n_CO_P_wanted) + "| T = " + T_wanted_str + "| P = " + P_wanted_str + " | CO = " + CO_wanted_str + "\n")

        # CARBON DIOXIDE

        file_abundances_final_selection2.write("# CARBON DIOXIDE \n")

        n_CO2_P_wanted = n_CO_P_wanted * n_H2O_P_wanted / K2_sin_presion

        file_abundances_final_selection2.write(
            str(n_CO2_P_wanted) + "| T = " + T_wanted_str + "| P = " + P_wanted_str + " | CO = " + CO_wanted_str + "\n")

        file_abundances_final_selection2.close()

        print("\nThe abundances for T="+T_wanted_str+", P="+P_wanted_str+" and C/O="+CO_wanted_str+" are:\n")
        print("ñ_CH4 = ", n_CH4_P_wanted)
        print("ñ_C2H2 = ", n_C2H2_P_wanted)
        print("ñ_H2O = ", n_H2O_P_wanted)
        print("ñ_CO = ", n_CO_P_wanted)
        print("ñ_CO2 = ", n_CO2_P_wanted)
        print("\nYou can find them in the file 'abundances_final_selection2_CO=" + str(int(CO_wanted)) + "_P_" + P_wanted_str + "_T_" + T_wanted_str + ".dat' \n")

        return n_CH4_P_wanted, n_C2H2_P_wanted, n_H2O_P_wanted, n_CO_P_wanted, n_CO2_P_wanted


def graphics_CO(P, figN, CO_wanted):  # P = "1bar" o P_wanted_str o ""

    T = np.linspace(500, 3000, 26)

    plt.figure()
    if P == "1bar":
        plt.plot(T, n_CH4_1bar, label='$CH_{4} (1bar)$', color='black', ls=':', linewidth=3)
        plt.plot(T, n_C2H2_1bar, label='$C_{2}H_{2} (1bar)$', color='y', ls=':', linewidth=3)
        plt.plot(T, n_H2O_1bar, label='$H_{2}O (1bar)$', color='m', ls=':', linewidth=3)
        plt.plot(T, n_CO_1bar, label='$CO (1bar)$', color='c', ls=':', linewidth=3)
        plt.plot(T, n_CO2_1bar, label='$CO_{2} (1bar)$', color='g', ls=':', linewidth=3)
    if P == P_wanted_str:
        plt.plot(T, n_CH4_P_wanted, label='$CH_{4} (P_wanted)$', color='black', ls=':', linewidth=1)
        plt.plot(T, n_C2H2_P_wanted, label='$C_{2}H_{2} (P_wanted)$', color='y', ls=':', linewidth=1)
        plt.plot(T, n_H2O_P_wanted, label='$H_{2}O (P_wanted)$', color='m', ls=':', linewidth=1)
        plt.plot(T, n_CO_P_wanted, label='$CO (P_wanted)$', color='c', ls=':', linewidth=1)
        plt.plot(T, n_CO2_P_wanted, label='$CO_{2} (P_wanted)$', color='g', ls=':', linewidth=1)
    if P == "":
        plt.plot(T, n_CH4_1bar, label='$CH_{4} (1bar)$', color='black', ls=':', linewidth=3)
        plt.plot(T, n_C2H2_1bar, label='$C_{2}H_{2} (1bar)$', color='y', ls=':', linewidth=3)
        plt.plot(T, n_H2O_1bar, label='$H_{2}O (1bar)$', color='m', ls=':', linewidth=3)
        plt.plot(T, n_CO_1bar, label='$CO (1bar)$', color='c', ls=':', linewidth=3)
        plt.plot(T, n_CO2_1bar, label='$CO_{2} (1bar)$', color='g', ls=':', linewidth=3)
        plt.plot(T, n_CH4_P_wanted, label='$CH_{4} (' + P_wanted_str + 'bar)$', color='black', ls='-', linewidth=1)
        plt.plot(T, n_C2H2_P_wanted, label='$C_{2}H_{2} (' + P_wanted_str + 'bar)$', color='y', ls='-', linewidth=1)
        plt.plot(T, n_H2O_P_wanted, label='$H_{2}O (' + P_wanted_str + 'bar)$', color='m', ls='-', linewidth=1)
        plt.plot(T, n_CO_P_wanted, label='$CO (' + P_wanted_str + 'bar)$', color='c', ls='-', linewidth=1)
        plt.plot(T, n_CO2_P_wanted, label='$CO_{2} (' + P_wanted_str + 'bar)$', color='g', ls='-', linewidth=1)
    plt.xlim([500, 3000])
    plt.ylim([10 ** -22, 10 ** 0])
    plt.yscale('log')
    plt.xticks(np.linspace(500, 3000, 6, endpoint=True))
    plt.yticks(10 ** (np.linspace(-22, 0, 12)))
    plt.tick_params(direction='in', bottom='on', top='on', left='on', right='on')
    plt.xlabel('Temperature (K)')
    plt.ylabel('$ñ_{x}$')
    plt.title('C/O = ' + str(CO_wanted))
    plt.legend(loc=4, prop={'size': 7}, ncol=2)
    plt.savefig('abundances-nist-' + str(figN) + '.png')
    plt.show()

# Now inputs

# Enter your Selection

selection = input("What do you want to do? \n [1] Generate a graphic of abundances on a range of temperatures for a specific C/O and pressure \n [2] Obtain the specific abundances of molecules on a specific temperature, C/O and pressure \n Answer ('1' or '2'): ")

# Enter your Pressure

P_wanted_str = input("Enter the pressure wanted (bar) : ")
P_wanted = float(P_wanted_str)

# Enter your C/O

CO_wanted_str = input("Enter the C/O wanted : ")
CO_wanted = float(CO_wanted_str)

if selection == '1':

    K1_1bar, K1_P_wanted, K2_sin_presion, K3_1bar, K3_P_wanted = classic_calculation()

    # NOW WE GENERATES MORE GRAPHICS

    n_CH4_1bar, n_CH4_P_wanted, n_C2H2_1bar, n_C2H2_P_wanted, n_H2O_1bar, n_H2O_P_wanted, n_CO_1bar, n_CO_P_wanted, n_CO2_1bar, n_CO2_P_wanted = abundances_norm_H(CO_wanted, P_wanted_str, selection)

    graphics_CO("", "fig_CO_" + CO_wanted_str + "_P_" + P_wanted_str, '{0:.2f}'.format(CO_wanted))

    print("The files of abundances and graphics for C/O=" + CO_wanted_str + " and P=" + P_wanted_str + " are ready! :D \n")


elif selection == '2':

    T_wanted_str = input("Enter the temperature wanted (between 500K and 6000K): ")
    T_wanted = float(T_wanted_str)

    K1_P_wanted, K2_sin_presion, K3_P_wanted = insert_one_temp(T_wanted, P_wanted)

    n_CH4_P_wanted, n_C2H2_P_wanted, n_H2O_P_wanted, n_CO_P_wanted, n_CO2_P_wanted = abundances_norm_H(CO_wanted, P_wanted_str, selection)

    print("The files of abundances for C/O=" + CO_wanted_str + " P=" + P_wanted_str + " and T=" + T_wanted_str + " are ready! :D \n")




