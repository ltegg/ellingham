# %% Compatability for Python 2.7.
from __future__ import division
from __future__ import unicode_literals
# %%
# ################
# ellingham.py
# ################
#
# by l t e g g @ l i v e . c o m  
# 2017-12-01
# Version 0.1
#
# A Python 3.7.0 script which generates Ellingham diagrams for some oxides,
# carbides, nitrides and fluorides, chlorides. Produces a .pdf which can be
# printed, and which can also be found in the github repository
# (https://github.com/ltegg/ellingham).
#
# As the change in Gibbs free energy is approximately linear with temperature,
# data is stored as a set of line segments, coded as pairs of co-ordinates.
# Line segments are stored separately for when the metal or the compound are in
# different phases.
#
# As the code is used only to generate a graph, no efforts have been made to 
# improve performance.
#
# Thermodynamic data for the oxides, nitrides, fluorides and chlorides comes
# from:
#  [1] Reed, T.B., 1971, Free Energy of Formation of Binary Compounds: An Atlas
#      of Charts for High-Temperature Chemical Calculations. MIT Press, 
#      Cambridge, Mass.
# Thermodynamic data for the carbides comes from
#  [2] Coltters, R.G., 1985, Thermodyanmics of binary metallic carbides: A
#      review. Materials Science and Engineering 76, 1-50.


# %% Import libraries
# Numpy for number handling
import numpy as np
# Patches to draw a rectangle
import matplotlib.patches as patches
# Pyplot to make plots
import matplotlib.pyplot as pl
# Use Arial as the font, and use regular text instead of math text for equations.
pl.rcParams.update({'mathtext.default': 'regular',
                    'mathtext.fontset':'custom',
                    'mathtext.it':'Arial:italic',
                    'mathtext.rm':'Arial',                 
                    'font.family':'Arial',})



    
# %% Metal Oxides

# Temperature values are in Kelvin
# DeltaG values is in kcal/moleO2
# Both are converted to other units later in the script.

# Metal: solid, oxide: solid
oxss = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                    [   0, 932, -266.6, -220.0, r'$\frac{4}{3} Al + O_2 = \frac{2}{3} Al_2O_3$', -13],
                    [   0, 904, -111.0,  -74.0, r'$\frac{4}{3} Sb + O_2 = \frac{2}{3} Sb_2O_3$', 3],
                    [   0, 983, -265.0, -222.0, '$2Ba + O_2 = 2BaO $', 0],
                    [   0, 544,  -92.0,  -69.0, r'$\frac{4}{3} Bi + O_2 = \frac{2}{3} Bi_2O_3$', 6],
                    [   0, 723, -200.5, -171.5, r'$\frac{4}{3} B + O_2 = \frac{2}{3} B_2O_3$', -3],
                    [   0,1123, -303.0, -249.0, '$2Ca + O_2 = 2CaO $', 0],
                    [   0,   0,  -55.6,  -55.6, '$2C + O_2 = 2CO $', 0],
                    [   0,   0,  -94.5,  -94.5, '$C + O_2 = CO_2 $', -8], 
                    [   0, 302, -151.8, -125.0, '$4Cs + O_2 = 2Cs_2O$', -14],
                    [   0,1357,  -80.0,  -33.0, '$4Cu + O_2 = 2Cu_2O $', 0],
                    [   0,1357,  -74.5,  -16.0, '$2Cu + O_2 = 2CuO $', 0],                   
                    [   0,   0, -119.3, -119.3, '$4H + O_2 = 2H_2O$', 12],
                    [   0,1642, -124.1,  -75.0, '$2Fe + O_2 = 2FeO$', -9],
                    [   0,1809, -129.2,  -55.5, r'$\frac{4}{3} Fe + O_2 = \frac{2}{3} Fe_2O_3$', -5],
                    [   0, 453, -286.0, -258.0, '$4Li + O_2 = 2Li_2O $', -12],
                    [   0,1068, -120.0,  -77.0, r'$ \frac{2}{3}Mo + O_2 = \frac{2}{3}MoO_3 $', -3],
                    [   0, 923, -286.0, -240.0, '$2Mg + O_2 = 2MgO $', 8],
                    [   0,   0,  -44.0,  -44.0, '$2Hg + O_2 = 2HgO$', 0],
                    [   0,1764, -181.0, -112.0, r'$\frac{4}{5}Nb + O_2 = \frac{2}{5}Nb_2O_5$', 0],
                    [   0, 734,  -32.0,    0.0, r'$\frac{3}{2} Pt + O_2 = \frac{1}{2}Pt_3O_4 $', 0],
                    [   0, 336, -172.0, -151.0, '$4K + O_2 = 2K_2O $', -14],
                    [   0, 312, -157.8, -138.0, '$4Rb + O_2 = 2Rb_2O$', -8],
                    [   0,1685, -216.5, -145.8, '$Si + O_2 = SiO_2 $', 0],
                    [   0, 480,  -14.0,    0.0, '$4Ag + O_2 = 2Ag_2O $', 0],
                    [   0, 371, -197.0, -176.0, '$4Na + O_2 = 2Na_2O $', 3],
                    [   0,1043, -281.0, -233.0, '$2Sr + O_2 = 2SrO$', 5],
                    [   0, 505, -138.8, -114.0, '$Sn + O_2 = SnO_2 $', -11],
                    [   0,1940, -225.5, -142.5, '$Ti + O_2 = TiO_2 $', 0],
                    [   0,1940, -247.5, -161.0, '$2Ti + O_2 = 2TiO$', 0],
                    [   0,1818, -168.0, -100.0, '$V + O_2 = VO_2$', -9],
                    [   0, 943, -149.5, -110.0, r'$\frac{4}{5}V + O_2 = \frac{2}{5}V_2O_5$', 2],                   
                    [   0,1743, -133.0,  -67.0, r'$ \frac{2}{3}W + O_2 = \frac{2}{3}WO_3 $', -10],
                    [   0, 693, -166.0, -134.0, '$2Zn + O_2 = 2ZnO $', 4],
                    [   0,2125, -262.0, -166.0, '$Zr + O_2 = ZrO_2 $', 9],
                ])

# Metal: liquid, oxide: solid
oxls = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                    [ 932,2345, -220.0, -147.6, r'$\frac{4}{3} Al + O = \frac{2}{3} Al_2O_3$', 0],
                    [ 904, 928,  -74.0,  -73.0, r'$\frac{4}{3} Sb + O_2 = \frac{2}{3} Sb_2O_3$', 0],                    
                    [ 983,1895, -222.0, -183.0, '$2Ba + O_2 = 2BaO $', 0],
                    [ 544,1098,  -69.0,  -44.0, r'$\frac{4}{3} Bi + O_2 = \frac{2}{3} Bi_2O_3$', 0],
                    [1123,1756, -249.0, -217.0, '$2Ca + O_2 = 2CaO $', 0],
                    [ 302, 763, -125.0,  -84.0, '$4Cs + O_2 = 2Cs_2O$', -16],                   
                    [1357,1509,  -33.0,  -28.0, '$4Cu + O_2 = 2Cu_2O $', 0],
                    [1357,1609,  -16.0,   -9.5, '$2Cu + O_2 = 2CuO $', 0], 
                    [ 453,1597, -258.0, -173.0, '$4Li + O_2 = 2Li_2O $', 0],
                    [ 336, 980, -151.0, -107.0, '$4K + O_2 = 2K_2O $', 0],
                    [   0, 630,  -44.0,  -10.0, '$2Hg + O_2 = 2HgO$', 0],
                    [ 312, 910, -138.0,  -96.0, '$4Rb + O_2 = 2Rb_2O$', -8],                   
                    [1685,1696, -145.8, -145.4, '$Si + O_2 = SiO_2 $', 0],
                    [ 371,1156, -176.0, -122.0, '$4Na + O_2 = 2Na_2O $', 0],
                    [1940,2128, -142.5, -134.5, '$Ti + O_2 = TiO_2 $', 0],  
                    [1940,2033, -161.0, -159.0, '$2Ti + O_2 = 2TiO$', 0],
                    [ 693,1180, -134.0, -109.0, '$2Zn + O_2 = 2ZnO $', 0],                   
                    [2125,2980, -166.0, -130.0, '$Zr + O_2 = ZrO_2 $', 0],                   
                ])

# Metal: gas, oxide: solid
oxgs = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                    [1895,2191, -183.0, -159.0, '$2Ba + O_2 = 2BaO $', 0],
                    [1756,2887, -217.0, -117.0, '$2Ca + O_2 = 2CaO $', 0],
                    [1597,2000, -173.0, -128.0, '$4Li + O_2 = 2Li_2O $', 0],
                    [ 923,1376, -240.0, -214.0, '$2Mg + O_2 = 2MgO $', 0],
                    [ 630, 740,  -10.0,    0.0, '$2Hg + O_2 = 2HgO$', 0],                    
                    [1156,1193, -122.0, -119.0, '$4Na + O_2 = 2Na_2O $', 0],
                    [1180,2240, -109.0,   -9.0, '$2Zn + O_2 = 2ZnO $', 0],                   
                    
                ])

# Metal: solid, oxide: liquid
oxsl = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                    [ 723,2313, -171.5, -112.0, r'$\frac{4}{3} B + O = \frac{2}{3} B_2O_3$', 0],
                    [1642,1809,  -75.0,  -71.9, '$2Fe + O_2 = 2FeO$', 0],
                    [1068,1530,  -77.0,  -64.0, '$ something Mo $', 0],
                    [1818,2190, -100.0,  -96.0, '$V + O_2 = VO_2$', 0],                    
                    [1743,2100,  -67.0,  -57.0, r'$ \frac{2}{3}W + O_2 = \frac{2}{3}WO_3 $', 0],                   
                ])              

# Metal: liquid, oxide: liquid
oxll = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                    [2345,2736, -147.6, -128.5, r'$\frac{4}{3} Al + O = \frac{2}{3} Al_2O_3$', 0],
                    [ 928,1698,  -73.0,  -45.0, r'$\frac{4}{3} Sb + O_2 = \frac{2}{3} Sb_2O_3$', 0],                                        
                    [1098,1852, -44.0,   -12.0, r'$\frac{4}{3} Bi + O_2 = \frac{2}{3} Bi_2O_3$', 0],
                    [1809,2000, -71.9,   -67.9, '$2Fe + O_2 = 2FeO$', 0],  
                    [1376,3125, -214.0,   52.0, '$2Mg + O_2 = 2MgO $', 0],
                    [ 763, 915,  -84.0,  -73.0, '$4Cs + O_2 = 2Cs_2O$', -16],                                       
                    [1509,2500,  -28.0,   -9.5, '$4Cu + O_2 = 2Cu_2O $', 0],
                    [1609,1870,   -9.5,      0, '$2Cu + O_2 = 2CuO $', 0],  
                    [ 980,1031, -107.0, -104.0, '$4K + O_2 = 2K_2O $', 0],
                    [ 910, 952,  -96.0,  -95.0, '$4Rb + O_2 = 2Rb_2O$', -8],                                       
                    [1696,2500, -145.4, -107.8, '$Si + O_2 = SiO_2 $', 0],
                    [2128,2500, -134.5, -121.5, '$Ti + O_2 = TiO_2 $', 0],
                    [2033,2500, -159.0, -142.5, '$2Ti + O_2 = 2TiO$', 0],
                    [2190,2500,  -96.0,  -81.0, '$V + O_2 = VO_2$', 0],                    
                    
                ])                    

# Metal: gas, oxide: liquid
oxgl = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                    [2191,2500, -159.0, -131.0, '$2Ba + O_2 = 2BaO $', 0],
                    [1031,1325, -104.0,  -71.0, '$4K + O_2 = 2K_2O $', 0],
                    [1193,1600, -119.0,  -62.0, '$4Na + O_2 = 2Na_2O $', 0], 
                    [2240,2340,   -9.0,    0.0, '$2Zn + O_2 = 2ZnO $', 0],                   
                    
                ])

# Metal: solid, oxide: gas
oxsg = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [   0,3400,  -55.6, -191.9, '$2C + O_2 = 2CO $', 0],
                [   0,3400,  -94.5,  -94.5, '$C + O_2 = CO_2 $', 0],
                [1530,2500,  -64.0,  -52.0, '$ something Mo $', 0],
                [2100,2500,  -57.0,  -52.0, r'$ \frac{2}{3}W + O_2 = \frac{2}{3}WO_3 $', 0],                
        ])

# Metal: liquid, oxide: gas
oxlg = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [ 915, 955,  -73.0,  -72.0, '$4Cs + O_2 = 2Cs_2O$', -16], 
                [1698,1908,  -45.0,  -32.0, r'$\frac{4}{3} Sb + O_2 = \frac{2}{3} Sb_2O_3$', 0],                                        
                
        ])

# Metal: gas, oxide: gas
oxgg = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [   0,3400, -119.3,  -26.6, '$4H + O_2 = 2H_2O$', 0],
                [1325,2160,  -71.0,    0.0, '$4K + O_2 = 2K_2O $', 0],
                [1600,2250,  -62.0,    0.0, '$4Na + O_2 = 2Na_2O $', 0],
                [1908,2380,  -32.0,    0.0, r'$\frac{4}{3} Sb + O_2 = \frac{2}{3} Sb_2O_3$', 0],                                        
                
                ])

    
# %% Metal Carbides

# Temperature values are in degrees Celcius
# DeltaG values is in kJ/moleC

# Metal: solid, carbide: solid
cass = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                    [   0, 1414,   -57,    -49, '$ Si + C = SiC $', -9],   
                    [   0, 1750,  -160, -150.5, '$ Ti + C = TiC $', -5],
                    [   0,  723,    23,     -1, '$3Fe + C = Fe_3C $', 0],
                    [   0, 1290,   -31,    -34, '$2W + C = W_2C$', -1],
                    [   0,  800, -39.5,    -45, '$W + C = WC$', -10],
                    [   0, 1000,   -70,    -59, '$2Mo + C = Mo_2C$', -12],
                    [   0,  720,  -183,   -175, '$Zr + C = ZrC$', 1],
                    
                    
                ])

# Metal: liquid, carbide: solid
cals = np.array([  #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                 [1414,2000,    -49,    -30, '$ Si + C = SiC $', 0],                    
                 ])

# Metal: gas, carbide: solid
cags = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form

# Metal: solid, carbide: liquid
casl = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form          

# Metal: liquid, carbide: liquid
call = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form                 

# Metal: gas, carbide: liquid
cagl = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form

# Metal: solid, carbide: gas
casg = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form

# Metal: liquid, carbide: gas
calg = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form

# Metal: gas, carbide: gas
cagg = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form


# %% Metal Nitrides

# Temperature values are in Kelvin
# DeltaG values is in kcal/moleO2
# Both are converted later in the script.

# Metal: solid, nitride: solid
niss = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                    [   0, 932, -144.3, -101.0, '$2Al + N_2 = 2AlN $', 0],
                    [   0,2300, -121.4,  -20.8, '$2B + N_2 = 2BN$', 0],
                    [   0,1809,   -5.8,   38.5, '$8Fe + N_2 = 2Fe_4N $', 12],
                    [   0, 923, -109.6,  -65.8, '$3Mg + N_2 = Mg_3N_2 $', 8],
                    [   0,1150,  -31.9,    0.0, '$4Mo + N_2 = Mo_2N $', 10],
                    [   0,0, -24.1,  -24.1, '$6H + N_2 = 2NH_3$', 0],
                    [   0,1680,  -90.0,  -22.5, r'$\frac{3}{2}Si + N_2 = \frac{1}{2}Si_3N_4 $', -6],
                    [   0,1940, -160.5,  -73.4, '$2Ti + N_2 = 2TiN $', 1],
                    [   0,2190,  -83.3,    3.6, '$2V + N_2 = 2VN$', -13],
                    [   0,2128, -163.8,  -67.2, '$2Zr + N_2 = 2ZrN $', -2],
                    
                ])

# Metal: liquid, nitride: solid
nils = np.array([     #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                    [2300,2500, -20.8,       0, '$2B + N_2 = 2BN$', 0],                    
                    [ 923,1376,  -65.8,  -41.3, '$3Mg + N_2 = Mg_3N_2', 0],
                    [1680,2130,  -22.5,    0.0, r'$\frac{3}{2}Si + N_2 = \frac{1}{2}Si_3N_4 $', 0],
          
                ])

# Metal: gas, nitride: solid
nigs = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form        

# Metal: solid, nitride: liquid
nisl = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form                     

# Metal: liquid, nitride: liquid
nill = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form                         

# Metal: gas, nitride: liquid
nigl = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form        

# Metal: solid, nitride: gas
nisg = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form        

# Metal: liquid, nitride: gas
nilg = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form        

# Metal: gas, nitride: gas
nigg = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [   0,2000,  -24.1,   85.2, '$6H + N_2 = 2NH_3$', 0],
                ])

    
# %% Metal Fluorides

# Temperature values are in Kelvin
# DeltaG values is in kcal/moleO2
# Both are converted later in the script.

# Metal: solid, fluoride: solid
flss = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [   0, 932, -215.3, -181.0, r'$\frac{2}{3}Al + F_2 = \frac{2}{3}AlF_3 $', 0],
                [   0,1123, -288.0, -245.0, '$ Ca + F_2 = CaF_2 $', 5],
                [   0,   0, -129.8, -129.8, '$2H + F_2 = 2HF $', 0],
                [   0,   0,  -81.2,  -81.2, r'$\frac{1}{2}C + F_2 = \frac{1}{2}CF_4 $', 2],
                [   0, 453, -290.0, -271.0, '$2Li + F_2 = 2LiF $', -8],
                [   0, 336, -270.0, -253.0, '$2K + F_2 = 2KF $', +0.3],
                [   0, 371, -274.0, -255.0, '$2Na + F_2 = 2NaF $', -0.3],             
                 ])

# Metal: liquid, fluoride: solid
flls = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [ 932,1545, -181.0, -156.0, r'$\frac{2}{3}Al + F_2 = \frac{2}{3}AlF_3 $', 0],
                [1123,1691, -245.0, -224.0, '$ Ca + F_2 = CaF_2 $', 0],
                [ 453,1120, -271.0, -240.0, '$2Li + F_2 = 2LiF $', 0],
                [ 336,1031, -253.0, -214.0, '$2K + F_2 = 2KF $', 0],
                [ 371,1187, -255.0, -214.0, '$2Na + F_2 = 2NaF', 0],                
                ])

# Metal: gas, fluoride: solid
flgs = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form        

# Metal: solid, fluoride: liquid
flsl = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form                     

# Metal: liquid, fluoride: liquid
flll = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [1545,2500, -156.0, -157.0, r'$\frac{2}{3}Al + F_2 = \frac{2}{3}AlF_3 $', 0],
                [1691,1955, -224.0, -222.0, '$ Ca + F_2 = CaF_2 $', 0],
                [1120,1597, -240.0, -216.0, '$2Li + F_2 = 2LiF $', 0],
                [1031,1130, -214.0, -209.0, '$2K + F_2 = 2KF $', 0],
                [1187,1268, -214.0, -209.0, '$2Na + F_2 = 2NaF', 0],                
                ])                 

# Metal: gas, fluoride: liquid
flgl = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [1955,2500, -222.0, -186.0, '$ Ca + F_2 = CaF_2 $', 0],
                [1597,1954, -216.0, -194.0, '$2Li + F_2 = 2LiF $', 0],
                [1130,1775, -209.0, -166.0, '$2K + F_2 = 2KF $', 0],
                [1268,1977, -209.0, -156.0, '$2Na + F_2 = 2NaF', 0],               
                ])

# Metal: solid, fluoride: gas
flsg = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form        

# Metal: liquid, fluoride: gas
fllg = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form        

# Metal: gas, fluoride: gas
flgg = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [   0,2500,  -81.2,  -36.0, r'$\frac{1}{2}C + F_2 = \frac{1}{2}CF_4 $', 0],
                [   0,   0, -129.8, -129.8, '$2H + F_2 = 2HF $', 0],               
                [1954,2500, -194.0, -179.0, '$2Li + F_2 = 2LiF $', 0],
                [1775,2500, -166.0, -150.0, '$2K + F_2 = 2KF $', 0],
                [1977,2500, -156.0, -150.0, '$2Na + F_2 = 2NaF', 0],
                [   0,1287, -129.8, -134.1, '$2H + F_2 = 2HF $', 0],                
                ])


# %% Metal Chlorides

# Temperature values are in Kelvin
# DeltaG values is in kcal/moleO2
# Both are converted later in the script.

# Metal: solid, chloride: solid
clss = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [   0, 465, -110.9,  -92.9, r'$\frac{2}{3}Al + Cl_2 = \frac{2}{3}AlCl_3 $', -5],
                [   0,1055, -188.0, -154.0, '$ Ca + Cl_2 = CaCl_2 $', 0],
                [   0,   0,  -12.3,  -12.3, r'$ \frac{1}{2}C + Cl_2 = \frac{1}{2}CCl_4 $  ', 0],
                [   0,   0,  -45.0,  -45.0, '$2H + Cl_2 = 2HCl $', -11],
                [   0, 459, -193.6, -177.6, '$2Li + Cl_2 = 2LiCl $', 0],
                [   0, 336, -209.4, -193.2, '$2K + Cl_2 = 2KCl $', 0],
                [   0, 371, -196.8, -180.0, '$2Na + Cl_2 = 2NaCl $', -5],
                [   0,   0,  -36.1,  -36.1, r'$\frac{1}{3} W + Cl_2 = \frac{1}{3} WCl_6 $', 8],
                 ])

# Metal: liquid, chloride: solid
clls = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [ 459, 887, -177.6, -161.0, '$2Li + Cl_2 = 2LiCl $', 0],
                [ 336,1031, -193.2, -161.0, '$2K + Cl_2 = 2KCl $', 0],
                [ 371,1073, -180.0, -149.4, '$2Na + Cl_2 = 2NaCl $', 0],
                
                ])

# Metal: gas, chloride: solid
clgs = np.array([[0, 0, 0, 0, ' ', 0]]) # no data for compounds in this form        

# Metal: solid, chloride: liquid
clsl = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [ 465, 500,  -92.9,  -91.7, r'$\frac{2}{3}Al + Cl_2 = \frac{2}{3}AlCl_3 $', 0],
                [1055,1123, -154.0, -152.0, '$ Ca + Cl_2 = CaCl_2 $', 0],
                [   0, 548,  -36.1,  -15.0, r'$\frac{1}{3} W + Cl_2 = \frac{1}{3} WCl_6 $', 0],
                
                ])             

# Metal: liquid, chloride: liquid
clll = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [1123,1755, -152.0, -136.0, '$ Ca + Cl_2 = CaCl_2 $', 0],
                [ 887,1597, -161.0, -141.2, '$2Li + Cl_2 = 2LiCl $', 0],
                [1031,1043, -161.0, -160.0, '$2K + Cl_2 = 2KCl $', 0],
                [1073,1156, -149.4, -145.6, '$2Na + Cl_2 = 2NaCl $', 0],
               
                ])                 

# Metal: gas, chloride: liquid
clgl = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [1755,1900, -136.0, -128.0, '$ Ca + Cl_2 = CaCl_2 $', 0],
                [1597,1655, -141.2, -138.4, '$2Li + Cl_2 = 2LiCl $', 0],
                [1043,1680, -160.0, -122.4, '$2K + Cl_2 = 2KCl $', 0],
                [1156,1738, -145.6, -110.0, '$2Na + Cl_2 = 2NaCl $', 0],
                
                ])

# Metal: solid, chloride: gas
clsg = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [ 500, 932,  -91.7,  -84.6, r'$\frac{2}{3}Al + Cl_2 = \frac{2}{3}AlCl_3 $', 0], 
                [ 548,1500,  -15.0,   -0.8, r'$\frac{1}{3} W + Cl_2 = \frac{1}{3} WCl_6 $', 0],
                
                ])

# Metal: liquid, chloride: gas
cllg = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [ 932,2273,  -84.6,  -70.2, r'$\frac{2}{3}Al + Cl_2 = \frac{2}{3}AlCl_3 $', 0],               
                
                ])

# Metal: gas, chloride: gas
clgg = np.array([ #T0,  T1,     G0,     G1,  chemical equation, label offset in vertical direction
                [2273,2500,  -70.2,  -71.6, r'$\frac{2}{3}Al + Cl_2 = \frac{2}{3}AlCl_3 $', 0],    
                [1900,2500, -128.0, -114.0, '$ Ca + Cl_2 = CaCl_2 $', 0],
                [   0,2500,  -12.3,   27.4, r'$ \frac{1}{2}C + Cl_2 = \frac{1}{2}CCl_4 $', 0],
                [   0,2500,  -45.0,  -53.3, '$2H + Cl_2 = 2HCl $', 0],
                [1655,2500, -138.4, -118.4, '$2Li + Cl_2 = 2LiCl $', 0],
                [1680,2500, -122.4, -110.4, '$2K + Cl_2 = 2KCl $', 0],
                [1738,2500, -110.0,  -96.8, '$2Na + Cl_2 = 2NaCl $', 0],         
                ])


# %% Convert temperatures to Celcius
# Convert temperatures for the oxide data
oxss[:,0:2] = oxss[:,0:2].astype(float) - 273.15
oxls[:,0:2] = oxls[:,0:2].astype(float) - 273.15
oxgs[:,0:2] = oxgs[:,0:2].astype(float) - 273.15
oxsl[:,0:2] = oxsl[:,0:2].astype(float) - 273.15
oxll[:,0:2] = oxll[:,0:2].astype(float) - 273.15
oxgl[:,0:2] = oxgl[:,0:2].astype(float) - 273.15
oxsg[:,0:2] = oxsg[:,0:2].astype(float) - 273.15
oxlg[:,0:2] = oxlg[:,0:2].astype(float) - 273.15
oxgg[:,0:2] = oxgg[:,0:2].astype(float) - 273.15

# Convert temperatures for the nitride data
niss[:,0:2] = niss[:,0:2].astype(float) - 273.15
nils[:,0:2] = nils[:,0:2].astype(float) - 273.15
nigs[:,0:2] = nigs[:,0:2].astype(float) - 273.15
nisl[:,0:2] = nisl[:,0:2].astype(float) - 273.15
nill[:,0:2] = nill[:,0:2].astype(float) - 273.15
nigl[:,0:2] = nigl[:,0:2].astype(float) - 273.15
nisg[:,0:2] = nisg[:,0:2].astype(float) - 273.15
nilg[:,0:2] = nilg[:,0:2].astype(float) - 273.15
nigg[:,0:2] = nigg[:,0:2].astype(float) - 273.15

# Convert temperatures for the fluoride data
flss[:,0:2] = flss[:,0:2].astype(float) - 273.15
flls[:,0:2] = flls[:,0:2].astype(float) - 273.15
flgs[:,0:2] = flgs[:,0:2].astype(float) - 273.15
flsl[:,0:2] = flsl[:,0:2].astype(float) - 273.15
flll[:,0:2] = flll[:,0:2].astype(float) - 273.15
flgl[:,0:2] = flgl[:,0:2].astype(float) - 273.15
flsg[:,0:2] = flsg[:,0:2].astype(float) - 273.15
fllg[:,0:2] = fllg[:,0:2].astype(float) - 273.15
flgg[:,0:2] = flgg[:,0:2].astype(float) - 273.15

# Convert temperatures for the chloride data
clss[:,0:2] = clss[:,0:2].astype(float) - 273.15
clls[:,0:2] = clls[:,0:2].astype(float) - 273.15
clgs[:,0:2] = clgs[:,0:2].astype(float) - 273.15
clsl[:,0:2] = clsl[:,0:2].astype(float) - 273.15
clll[:,0:2] = clll[:,0:2].astype(float) - 273.15
clgl[:,0:2] = clgl[:,0:2].astype(float) - 273.15
clsg[:,0:2] = clsg[:,0:2].astype(float) - 273.15
cllg[:,0:2] = cllg[:,0:2].astype(float) - 273.15
clgg[:,0:2] = clgg[:,0:2].astype(float) - 273.15


# %% Convert free energy change to kJ
# Convert free energy for oxide data
oxss[:,2:4] = oxss[:,2:4].astype(float) * 4.184
oxls[:,2:4] = oxls[:,2:4].astype(float) * 4.184
oxgs[:,2:4] = oxgs[:,2:4].astype(float) * 4.184
oxsl[:,2:4] = oxsl[:,2:4].astype(float) * 4.184
oxll[:,2:4] = oxll[:,2:4].astype(float) * 4.184
oxgl[:,2:4] = oxgl[:,2:4].astype(float) * 4.184
oxsg[:,2:4] = oxsg[:,2:4].astype(float) * 4.184
oxlg[:,2:4] = oxlg[:,2:4].astype(float) * 4.184
oxgg[:,2:4] = oxgg[:,2:4].astype(float) * 4.184

# Convert free energy for nitride data
niss[:,2:4] = niss[:,2:4].astype(float) * 4.184
nils[:,2:4] = nils[:,2:4].astype(float) * 4.184
nigs[:,2:4] = nigs[:,2:4].astype(float) * 4.184
nisl[:,2:4] = nisl[:,2:4].astype(float) * 4.184
nill[:,2:4] = nill[:,2:4].astype(float) * 4.184
nigl[:,2:4] = nigl[:,2:4].astype(float) * 4.184
nisg[:,2:4] = nisg[:,2:4].astype(float) * 4.184
nilg[:,2:4] = nilg[:,2:4].astype(float) * 4.184
nigg[:,2:4] = nigg[:,2:4].astype(float) * 4.184

# Convert free energy for fluoride data
flss[:,2:4] = flss[:,2:4].astype(float) * 4.184
flls[:,2:4] = flls[:,2:4].astype(float) * 4.184
flgs[:,2:4] = flgs[:,2:4].astype(float) * 4.184
flsl[:,2:4] = flsl[:,2:4].astype(float) * 4.184
flll[:,2:4] = flll[:,2:4].astype(float) * 4.184
flgl[:,2:4] = flgl[:,2:4].astype(float) * 4.184
flsg[:,2:4] = flsg[:,2:4].astype(float) * 4.184
fllg[:,2:4] = fllg[:,2:4].astype(float) * 4.184
flgg[:,2:4] = flgg[:,2:4].astype(float) * 4.184

# Convert free energy for chloride data
clss[:,2:4] = clss[:,2:4].astype(float) * 4.184
clls[:,2:4] = clls[:,2:4].astype(float) * 4.184
clgs[:,2:4] = clgs[:,2:4].astype(float) * 4.184
clsl[:,2:4] = clsl[:,2:4].astype(float) * 4.184
clll[:,2:4] = clll[:,2:4].astype(float) * 4.184
clgl[:,2:4] = clgl[:,2:4].astype(float) * 4.184
clsg[:,2:4] = clsg[:,2:4].astype(float) * 4.184
cllg[:,2:4] = cllg[:,2:4].astype(float) * 4.184
clgg[:,2:4] = clgg[:,2:4].astype(float) * 4.184


# %% Make the oxide diagram
# Create the figure, use approximate A3 landscape sizes
fig, (a1, a2) = pl.subplots(1,2,figsize=(15.5, 10.5))

# Use dummy data to include a legend
a1.plot(-1000, 1000, color='r', ls='-',  alpha=1,   label='Metal solid, oxide solid')
a1.plot(-1000, 1000, color='r', ls='--', alpha=1,   label='Metal liquid, oxide solid')
a1.plot(-1000, 1000, color='r', ls=':',  alpha=1,   label='Metal gas, oxide solid')
a1.plot(-1000, 1000, color='r', ls='-',  alpha=0.6, label='Metal solid, oxide liquid')
a1.plot(-1000, 1000, color='r', ls='--', alpha=0.6, label='Metal liquid, oxide liquid')
a1.plot(-1000, 1000, color='r', ls=':',  alpha=0.6, label='Metal gas, oxide liquid')
a1.plot(-1000, 1000, color='r', ls='-',  alpha=0.3, label='Metal solid, oxide gas')
a1.plot(-1000, 1000, color='r', ls='--', alpha=0.3, label='Metal liquid, oxide gas')
a1.plot(-1000, 1000, color='r', ls=':',  alpha=0.3, label='Metal gas, oxide gas')


linedict = {'marker':'.', 'markersize':2.25}

for row in oxss:
    a1.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='r', ls='-',  alpha=1,)
for row in oxls:
    a1.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='r', ls='--',  alpha=1,)
for row in oxgs:
    a1.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='r', ls=':',  alpha=1,)
for row in oxsl:
    a1.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='r', ls='-',  alpha=0.6,)
for row in oxll:
    a1.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='r', ls='--', alpha=0.6,) 
for row in oxgl:
    a1.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='r', ls=':',  alpha=0.6,)
for row in oxsg:
    a1.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='r', ls='-',  alpha=0.3,)
for row in oxgg:
    a1.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='r', ls=':',  alpha=0.3,)
for row in oxlg:
    a1.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='r', ls='--',  alpha=0.3,)

# Add text labels for each reaction    
for row in oxss:
    a1.text(float(row[0]) - 25,float(row[2]) + float(row[5]), 
            row[4],
            horizontalalignment='right',
            verticalalignment='center',
            fontsize=8)

# Set ticks and axis limits
xticks = [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
yticks = np.arange(-1300, 100, 100)
a1.set_xlim([-800,2000])
a1.set_xticks(xticks)
a1.set_ylim([-1300, 50])
a1.set_yticks(yticks)

# Manually add vertical grid lines only
for line in xticks:
    a1.axvline(line,color='0.5',alpha=0.5,zorder=-9)

# Add T = 0 and G = 0 lines
a1.axvline(0,color='k')
a1.axhline(0,color='k')

# Title and axis titles
a1.set_title('Oxides',fontweight='bold')
a1.set_xlabel('Temperature (°C)',x=0.64)#fontweight='bold')
a1.set_ylabel(r'Standard free energy of formation ($\Delta \mathit{G}_f \! \degree$) kJ/$mol_{O_2}$',)#fontweight='bold')


# A long hack to draw out a legend
#  Add a rectangle
rectpos = [900,1970, -1290, -1060]
a1.add_patch(patches.Rectangle(
        (rectpos[0], rectpos[2]),
        rectpos[1]-rectpos[0],
        rectpos[3]-rectpos[2],
        facecolor='#ffffff',
        fill=True, edgecolor='k', linewidth=1
        )
        )
#  Add a series of text labels
a1.text((rectpos[1]-rectpos[0])/2 + rectpos[0] + 155,
        rectpos[3]-30,
        'Metal',
        horizontalalignment='center',
        fontsize=9,
        fontweight='bold',
        )
a1.text((rectpos[1]-rectpos[0])/4 + rectpos[0] + 170,
        rectpos[3]-65,
        'Solid',
        horizontalalignment='center',
        fontsize=9,
        )
a1.text((rectpos[1]-rectpos[0])/2 + rectpos[0] + 155,
        rectpos[3]-65,
        'Liquid',
        horizontalalignment='center',
        fontsize=9,
        )
a1.text(3*(rectpos[1]-rectpos[0])/4 + rectpos[0] + 140,
        rectpos[3]-65,
        'Gas',
        horizontalalignment='center',
        fontsize=9,        
        )
a1.text(rectpos[0]+ 70,
        rectpos[3]-110,
        'Compound',
        horizontalalignment='center',
        fontsize=9,
        rotation=90,
        fontweight='bold',
        )
a1.text(rectpos[0]+ 290,
        rectpos[3]-110,
        'Solid',
        horizontalalignment='right',
        fontsize=9,
        )
a1.text(rectpos[0]+ 290,
        rectpos[3]-155,
        'Liquid',
        horizontalalignment='right',
        fontsize=9,
        )
a1.text(rectpos[0]+ 290,
        rectpos[3]-200,
        'Gas',
        horizontalalignment='right',
        fontsize=9,        
        )

#  Add short lines to indicate the line style
a1.plot([1260, 1400], [-1160, -1160], color='r', ls='-',  alpha=1,   label='Metal solid, oxide solid')
a1.plot([1520, 1660], [-1160, -1160], color='r', ls='--', alpha=1,   label='Metal liquid, oxide solid')
a1.plot([1780, 1920], [-1160, -1160], color='r', ls=':',  alpha=1,   label='Metal gas, oxide solid')
a1.plot([1260, 1400], [-1208, -1208], color='r', ls='-',  alpha=0.6, label='Metal solid, oxide liquid')
a1.plot([1520, 1660], [-1208, -1208], color='r', ls='--', alpha=0.6, label='Metal liquid, oxide liquid')
a1.plot([1780, 1920], [-1208, -1208], color='r', ls=':',  alpha=0.6, label='Metal gas, oxide liquid')
a1.plot([1260, 1400], [-1255, -1255], color='r', ls='-',  alpha=0.3, label='Metal solid, oxide gas')
a1.plot([1520, 1660], [-1255, -1255], color='r', ls='--', alpha=0.3, label='Metal liquid, oxide gas')
a1.plot([1780, 1920], [-1255, -1255], color='r', ls=':',  alpha=0.3, label='Metal gas, oxide gas')



# %% Make the carbide, nitride, fluoride and chloride diagram

# Add in the line segments from the carbide data
for row in cass:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='0.4', ls='-',  alpha=1)
for row in cals:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='0.4', ls='--',  alpha=1,)
#for row in cags:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='0.4', ls=':',  alpha=1,)
#for row in casl:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='0.4', ls='-',  alpha=0.6)
#for row in call:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='0.4', ls='--', alpha=0.6,)  
#for row in cagl:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='0.4', ls=':',  alpha=0.6,)
#for row in casg:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='0.4', ls='-',  alpha=0.3,)    
#for row in casg:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='0.4', ls='--',  alpha=0.3,)
#for row in cagg:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='0.4', ls=':',  alpha=0.3,)

# Add text labels for each reaction 
for row in cass:
    a2.text(float(row[0]) - 25,float(row[2]) + float(row[5]), 
            row[4],
            horizontalalignment='right',
            verticalalignment='center',
            fontsize=8)

# Add in the line segments from the nitride data
for row in niss:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='b', ls='-',  alpha=1)
for row in nils:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='b', ls='--',  alpha=1,)
#for row in nigs:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='b', ls=':',  alpha=1,)
#for row in nisl:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='b', ls='-',  alpha=0.6)
#for row in nill:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='b', ls='--', alpha=0.6,)    
#for row in nigl:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='b', ls=':',  alpha=0.6,)
#for row in nisg:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='b', ls='-',  alpha=0.3,)    
#for row in nisg:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='b', ls='--',  alpha=0.3,)
#for row in nigg:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='b', ls=':',  alpha=0.3,)

# Add text labels for each reaction 
for row in niss:
    a2.text(float(row[0]) - 25,float(row[2]) + float(row[5]), 
            row[4],
            horizontalalignment='right',
            verticalalignment='center',
            fontsize=8)
    
# Add in the line segments from the fluoride data
for row in flss:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color=[0, 1, 0], ls='-',  alpha=1)
for row in flls:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color=[0, 1, 0], ls='--',  alpha=1,)
#for row in flgs:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color=[0, 1, 0], ls=':',  alpha=1,)
#for row in flsl:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color=[0, 1, 0], ls='-',  alpha=0.6)
for row in flll:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color=[0, 1, 0], ls='--', alpha=0.6,)    
for row in flgl:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color=[0, 1, 0], ls=':',  alpha=0.6,)
#for row in flsg:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color=[0, 1, 0], ls='-',  alpha=0.3,)    
for row in flsg:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color=[0, 1, 0], ls='--',  alpha=0.3,)
for row in flgg:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color=[0, 1, 0], ls=':',  alpha=0.3,)

# Add text labels for each reaction 
for row in flss:
    a2.text(float(row[0]) - 25,float(row[2]) + float(row[5]), 
            row[4],
            horizontalalignment='right',
            verticalalignment='center',
            fontsize=8)

# Add in the line segments from the chloride data
for row in clss:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='g', ls='-',  alpha=1)
for row in clls:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='g', ls='--',  alpha=1,)
#for row in clgs:
#    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='g', ls=':',  alpha=1,)
for row in clsl:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='g', ls='-',  alpha=0.6)
for row in clll:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='g', ls='--', alpha=0.6,)   
for row in clgl:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='g', ls=':',  alpha=0.6,)
for row in clsg:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='g', ls='-',  alpha=0.3,)    
for row in cllg:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='g', ls='--',  alpha=0.3,)
for row in clgg:
    a2.plot([float(row[0]), float(row[1])],[float(row[2]), float(row[3])], **linedict, color='g', ls=':',  alpha=0.3,)

# Add text labels for each reaction 
for row in clss:
    a2.text(float(row[0]) - 25,float(row[2]) + float(row[5]), 
            row[4],
            horizontalalignment='right',
            verticalalignment='center',
            fontsize=8)



# Set ticks and axis limits
a2.set_xlim([-800,2000])
a2.set_xticks(xticks)
a2.set_ylim([-1300, 50])
a2.set_yticks(yticks)

# Manually add vertical grid lines only
for line in xticks:
    a2.axvline(line,color='0.5',alpha=0.5,zorder=-9)

# Add T = 0 and G = 0 lines
a2.axvline(0,color='k')
a2.axhline(0,color='k')

# Title and axis titles
a2.set_title('Carbides, nitrides, fluorides and chlorides',fontweight='bold')
a2.set_xlabel('Temperature (°C)',x=0.64)#fontweight='bold')
a2.set_ylabel('Standard free energy of formation ($\Delta \mathit{G}_f \! \degree$) kJ/$mol_{C/N_2/F_2/Cl_2}$',)#fontweight='bold')

# Another long hack to draw out a legend
#  Add a rectangle
rectpos = [900,1970, -1290, -1060]
a2.add_patch(patches.Rectangle(
        (rectpos[0], rectpos[2]),
        rectpos[1]-rectpos[0],
        rectpos[3]-rectpos[2],
        facecolor='#ffffff',
        fill=True, edgecolor='k', linewidth=1
        #zorder=-3
        )
        )
#  Add a series of text labels    
a2.text(rectpos[0] + 30,
        rectpos[3]-25,
        'Sources',
        fontsize=9,
        fontweight='bold',
        )

a2.text(rectpos[0] + 30,
        rectpos[3]-30,
        '$O_2$, $N_2$, $F_2$ and $Cl_2$ data from:',
        fontsize=9,
        verticalalignment='top',
        )
a2.text(rectpos[0] + 30,
        rectpos[3]-35,
        '\nReed, T.B., 1971. Free energy of \nformation of binary compounds. \nMIT Press, Cambridge, Mass.',
        fontsize=8,
        verticalalignment='top',
        fontstyle='italic',
        ) 
a2.text(rectpos[0] + 30,
        rectpos[3]-30,
        '\n\n\n\nC data from:',
        fontsize=9,
        verticalalignment='top',
        )
a2.text(rectpos[0] + 30,
        rectpos[3]-44,
        '\n\n\n\n\nColtters, R.G., 1985. Thermodynamics \nof binary metallic carbides: A review. \nMaterials Science and Engineering \n76, 1–50.',
        fontsize=8,
        verticalalignment='top',
        fontstyle='italic',        
        ) 

# Tight layout
#pl.tight_layout()

# Save the figures
pl.savefig('ellingham.eps',dpi=400,bbox_inches='tight')
pl.savefig('ellingham.png',dpi=400,bbox_inches='tight')
pl.savefig('ellingham.pdf',dpi=400,bbox_inches='tight')