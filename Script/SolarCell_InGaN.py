import numpy as np
import matplotlib.pyplot as plt

from solcore import material, si
from solcore.solar_cell import SolarCell, Layer, Junction
from solcore.solar_cell_solver import solar_cell_solver
from solcore.interpolate import interp1d
from solcore.light_source import LightSource

wavelengths = si(np.linspace(250, 1800, 200), "nm")
light_source = LightSource(source_type='standard', x=wavelengths, version='AM1.5g')

InGaN = material("InGaN")

InGaN_n = InGaN(Nd=si("5e19cm-3"), hole_diffusion_length=si("10um"), electron_mobility=0.2740, hole_mobility=0.0289, relative_permittivity=11.7)
InGaN_p = InGaN(Na=si("1e15cm-3"), electron_diffusion_length=si("400um"), electron_mobility=0.2635, hole_mobility=0.0289, relative_permittivity=11.7)

emitter_layer = Layer(width=si("0.6um"), material=InGaN_n, role='emitter')
base_layer = Layer(width=si("1.8um"), material=InGaN_p, role='base')

InGaN_junction = Junction([emitter_layer, base_layer], kind="DA")

solar_cell= SolarCell([InGaN_junction])

solar_cell_solver(solar_cell, 'qe', user_options={'wavelength': wavelengths, 'optics_method': "TMM",
                                                  'recalculate_absorption': True})

plt.figure(1)
plt.plot(wavelengths, solar_cell[0].eqe(wavelengths) * 100, 'k-', label="InGaN (TMM)")
plt.plot(wavelengths, 100*solar_cell[0].layer_absorption, label='Absorption')
plt.legend()
plt.ylim(0, 100)
plt.ylabel('EQE (%)')
plt.xlabel('Wavelength (m)')
plt.tight_layout()
plt.title("EQE and absorption for single InGaN cell using TMM optical methods")


V = np.linspace(0, 2, 300)

solar_cell_solver(solar_cell, 'iv', user_options={'voltages': V, 'light_iv': True,
                                                  'wavelength': wavelengths, 'mpp': True,
                                                 'light_source': light_source,
                                                  'recalculate_absorption': True,
                                                  'optics_method': "TMM"})
plt.figure(2)


plt.plot(V, -solar_cell[0].iv(V)/10, 'b', label='InGaN')

plt.text(1.04, 33, 'Efficieny (%): ' + str(np.round(solar_cell.iv['Eta'] * 100, 1)))
plt.text(1.04, 31, 'FF (%): ' + str(np.round(solar_cell.iv['FF'] * 100, 1)))
plt.text(1.04, 29, r'V$_{oc}$ (V): ' + str(np.round(solar_cell.iv["Voc"], 2)))
plt.text(1.04, 27, r'I$_{sc}$ (A/m$^2$): ' + str(np.round(solar_cell.iv["Isc"]/10, 2)))

plt.legend()
plt.ylim(0, 40)
plt.xlim(0, 1.4)
plt.ylabel('Current (A/m$^2$)')
plt.xlabel('Voltage (V)')
plt.title("IV characteristics of Single InGaN cell")

plt.show()


