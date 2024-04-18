import numpy as np
import os
import matplotlib.pyplot as plt

from solcore import siUnits, material, si
from solcore.solar_cell import SolarCell
from solcore.structure import Junction, Layer
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
from solcore.absorption_calculator import search_db

wl = np.linspace(300, 1850, 700) * 1e-9

light_source = LightSource(source_type='standard', x=wl, version='AM1.5g')

MgF2_pageid = search_db(os.path.join("MgF2", "Rodriguez-de Marcos"))[0][0];
ZnS_pageid = search_db(os.path.join("ZnS", "Querry"))[0][0];
MgF2 = material(str(MgF2_pageid), nk_db=True)();
ZnS = material(str(ZnS_pageid), nk_db=True)();

ARC = [Layer(si("100nm"), MgF2), Layer(si("15nm"), ZnS), Layer(si("15nm"), MgF2), Layer(si("50nm"), ZnS)]

AlInP = material("AlInP")
InGaP = material("GaInP")
window_material = AlInP(Al=0.52)

top_cell_n_material = InGaP(In=0.49, Nd=siUnits(2e18, "cm-3"), hole_diffusion_length=si("200nm"))
top_cell_p_material = InGaP(In=0.49, Na=siUnits(1e17, "cm-3"), electron_diffusion_length=si("1um"))

GaAs = material("GaAs")

mid_cell_n_material = GaAs(Nd=siUnits(3e18, "cm-3"), hole_diffusion_length=si("500nm"))
mid_cell_p_material = GaAs(Na=siUnits(1e17, "cm-3"), electron_diffusion_length=si("5um"))

Ge = material("Ge")

bot_cell_n_material = Ge(Nd=siUnits(2e18, "cm-3"), hole_diffusion_length=si("800nm"), hole_mobility=0.01)
bot_cell_p_material = Ge(Na=siUnits(1e17, "cm-3"), electron_diffusion_length=si("50um"), electron_mobility=0.1)

solar_cell = SolarCell(
    ARC +
    [
        Junction([Layer(si("400nm"), material=bot_cell_n_material, role='emitter'),
                  Layer(si("100um"), material=bot_cell_p_material, role='base'),
                  ], sn=1, sp=1, kind='DA'),
    ], shading=0.05, R_series=2e-6)
    
position = len(solar_cell) * [0.1e-9]
position[-1] = 10e-9 # Indexing with -1 in a Python list/array gives you the last element

plt.figure(1)

# First calculate with TMM optical method
solar_cell_solver(solar_cell, 'qe', user_options={'wavelength': wl, 'optics_method': "TMM",
                                                  'position': position, 'recalculate_absorption': True})

#plt.plot(wl * 1e9, solar_cell[4].eqe(wl) * 100, 'b', label='GaInP (TMM)')
#plt.plot(wl * 1e9, solar_cell[4].eqe(wl) * 100, 'g', label='InGaAs (TMM)')
plt.plot(wl * 1e9, solar_cell[4].eqe(wl) * 100, 'r', label='Ge (TMM)')
plt.plot(wl * 1e9, 100 * (1 - solar_cell.reflected), 'k--', label='1-R (TMM)')

# Recalculate with simple Beer-Lambert (BL) law absorption to compare
solar_cell_solver(solar_cell, 'qe', user_options={'wavelength': wl, 'optics_method': "BL",
                                                  'position': position, 'recalculate_absorption': True})

#plt.plot(wl * 1e9, solar_cell[4].eqe(wl) * 100, 'b--', alpha=0.5, label='GaInP (BL)')
#plt.plot(wl * 1e9, solar_cell[4].eqe(wl) * 100, 'g--', alpha=0.5, label='InGaAs (BL)')
plt.plot(wl * 1e9, solar_cell[4].eqe(wl) * 100, 'r--', alpha=0.5, label='Ge (BL)')
plt.legend()
plt.ylim(0, 100)
plt.ylabel('EQE (%)')
plt.xlabel('Wavelength (nm)')
plt.tight_layout()
plt.title("(1) EQE and absorption for 3J cell using TMM and BL optical methods")
plt.show()

V = np.linspace(0, 3, 300)
solar_cell_solver(solar_cell, 'iv', user_options={'voltages': V, 'light_iv': True,
                                                  'wavelength': wl, 'mpp': True,
                                                  'light_source': light_source,
                                                  'recalculate_absorption': True,
                                                  'optics_method': "TMM"})

plt.figure(2)
plt.plot(V, solar_cell.iv['IV'][1], 'k', linewidth=3, label='Total')
#plt.plot(V, -solar_cell[4].iv(V), 'b', label='GaInP')
#plt.plot(V, -solar_cell[4].iv(V), 'g', label='InGaAs')
plt.plot(V, -solar_cell[4].iv(V), 'r', label='Ge')
plt.text(1.4, 220, 'Efficieny (%): ' + str(np.round(solar_cell.iv['Eta'] * 100, 1)))
plt.text(1.4, 200, 'FF (%): ' + str(np.round(solar_cell.iv['FF'] * 100, 1)))
plt.text(1.4, 180, r'V$_{oc}$ (V): ' + str(np.round(solar_cell.iv["Voc"], 2)))
plt.text(1.4, 160, r'I$_{sc}$ (A/m$^2$): ' + str(np.round(solar_cell.iv["Isc"], 2)))

plt.legend()
plt.ylim(0, 250)
plt.xlim(0, 3)
plt.ylabel('Current (A/m$^2$)')
plt.xlabel('Voltage (V)')
plt.title("(2) IV characteristics of 3J cell")

plt.show()
    
