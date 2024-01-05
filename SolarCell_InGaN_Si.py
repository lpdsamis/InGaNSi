import numpy as np
import matplotlib.pyplot as plt
import os



from solcore import siUnits, material, si
from solcore.solar_cell import SolarCell
from solcore.structure import Junction, Layer
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
from solcore.absorption_calculator import search_db

wavelengths = np.linspace(300, 1850, 700)* 1e-9 


light_source = LightSource(source_type='standard', x=wavelengths, version='AM1.5g')

MgF2_pageid = search_db(os.path.join("MgF2", "Rodriguez-de Marcos"))[0][0];
ZnS_pageid = search_db(os.path.join("ZnS", "Querry"))[0][0];
MgF2 = material(str(MgF2_pageid), nk_db=True)();
ZnS = material(str(ZnS_pageid), nk_db=True)();

Si = material("Si")
SiN = material("Si3N4")()
InGaN = material("InGaN")

### To minimize front surface reflection, we use a four-layer anti-reflection coating (ARC): ######

ARC = [Layer(si("100nm"), MgF2), Layer(si("15nm"), ZnS), Layer(si("15nm"), MgF2), Layer(si("50nm"), ZnS)]

############## TOP CELL - InGaN  ###########################

InGaN_n = InGaN(Nd=si("5e19cm-3"), hole_diffusion_length=si("10um"), electron_mobility=0.2740, hole_mobility=0.0289, relative_permittivity=11.7)
InGaN_p = InGaN(Na=si("1e15cm-3"), electron_diffusion_length=si("400um"), electron_mobility=0.2635, hole_mobility=0.0289, relative_permittivity=11.7)

########## BOTTOM CELL - Si ###########################

Si_n = Si(Nd=si("1e21cm-3"), hole_diffusion_length=si("10um"), electron_mobility=0.1, hole_mobility=0.04,relative_permittivity=11.7)
Si_p = Si(Na=si("1e16cm-3"), electron_diffusion_length=si("400um"), electron_mobility=0.1, hole_mobility=0.04,relative_permittivity=11.7)

##################### solar_cell ###############

solar_cell = SolarCell(
    ARC +
    [
        Junction([Layer(width=si("500nm"), material=InGaN_n, role='emitter'),
        	  Layer(width=si("1560nm"), material=InGaN_p, role='base'),
                  ], sn=1, sp=1, kind='DA'),
        Junction([Layer(width=si("1000nm"), material=Si_n, role='emitter'),
                  Layer(width=si("3000nm"), material=Si_p, role='base'),
                  ], sn=1, sp=1, kind='DA'),
       
      ], shading=0.08, R_series=2e-6)
      
################### Setting the depth spacing ###############################

position = len(solar_cell) * [0.1e-9] 
# position[-1] = 10e-9
      
##################### solar_cell_solver #############      

solar_cell_solver(solar_cell, 'qe', user_options={'wavelength': wavelengths, 'optics_method': "TMM",
                                                  'position': position, 'recalculate_absorption': True})

plt.figure(1)

plt.plot(wavelengths* 1e9, solar_cell[4].eqe(wavelengths) * 100, 'b', label='InGaN (TMM)')
plt.plot(wavelengths* 1e9, solar_cell[5].eqe(wavelengths) * 100, 'g', label='Si (TMM)')
plt.plot(wavelengths* 1e9, 100 * (1 - solar_cell.reflected), 'k--', label='1-R (TMM)')

plt.legend()
plt.ylim(0, 100)
plt.ylabel('EQE (%)')
plt.xlabel('Wavelength (nm)')
plt.tight_layout()
plt.title("EQE and absorption for 2J cell using TMM optical methods")


V = np.linspace(0, 3, 300)

solar_cell_solver(solar_cell, 'iv', user_options={'voltages': V, 'light_iv': True,
                                                  'wavelength': wavelengths, 'mpp': True,
                                                 'light_source': light_source,
                                                  'recalculate_absorption': True,
                                                  'optics_method': "TMM"})
plt.figure(2)

plt.plot(V, solar_cell.iv['IV'][1], 'k', linewidth=3, label='Total')
plt.plot(V, -solar_cell[4].iv(V), 'b', label='InGaN')
plt.plot(V, -solar_cell[5].iv(V), 'g', label='Si')

plt.text(1.4, 30, 'Efficieny (%): ' + str(np.round(solar_cell.iv['Eta'] * 100, 1)))
plt.text(1.4, 28, 'FF (%): ' + str(np.round(solar_cell.iv['FF'] * 100, 1)))
plt.text(1.4, 26, r'V$_{oc}$ (V): ' + str(np.round(solar_cell.iv["Voc"], 2)))
plt.text(1.4, 24, r'I$_{sc}$ (A/m$^2$): ' + str(np.round(solar_cell.iv["Isc"], 2)))

plt.legend()
plt.ylim(0, 40)
plt.xlim(0, 2)
plt.ylabel('Current (mA/cm$^2$)')
plt.xlabel('Voltage (V)')
plt.title("(2) IV characteristics of 2J cell")

plt.show()


