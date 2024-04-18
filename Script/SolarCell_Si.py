import numpy as np
import matplotlib.pyplot as plt

from solcore import material, si
from solcore.solar_cell import SolarCell, Layer, Junction
from solcore.solar_cell_solver import solar_cell_solver
from solcore.interpolate import interp1d
from solcore.light_source import LightSource

wavelengths = si(np.linspace(250, 1800, 200), "nm")
light_source = LightSource(source_type='standard', x=wavelengths, version='AM1.5g')

Si = material("Si")
SiN = material("Si3N4")()
GaAs = material("GaAs")

Si_n = Si(Nd=si("1e21cm-3"), hole_diffusion_length=si("10um"),electron_mobility=0.1, hole_mobility=0.04, relative_permittivity=11.7)
Si_p = Si(Na=si("1e16cm-3"), electron_diffusion_length=si("400um"), electron_mobility=0.1, hole_mobility=0.04,relative_permittivity=11.7)

emitter_layer = Layer(width=si("1.8um"), material=Si_n, role='emitter') #2um
base_layer = Layer(width=si("100um"), material=Si_p, role='base') # 100um

Si_junction = Junction([emitter_layer, base_layer], kind="DA")

solar_cell_TMM_ARC = SolarCell([Layer(width=si(75, "nm"), material=SiN), Si_junction])

solar_cell_solver(solar_cell_TMM_ARC, 'qe', user_options={'wavelength': wavelengths, 'optics_method': "TMM",
                                                  'recalculate_absorption': True})

plt.figure(1)
plt.plot(wavelengths, 100*solar_cell_TMM_ARC[1].eqe(wavelengths), 'k-', label="Si (TMM)")
plt.plot(wavelengths, 100*solar_cell_TMM_ARC[1].layer_absorption, label='Absorption')
plt.legend()
plt.ylim(0, 100)
plt.ylabel('EQE (%)')
plt.xlabel('Wavelength (m)')
plt.tight_layout()
plt.title("EQE and absorption for single Si cell using TMM optical methods")


V = np.linspace(0, 2, 300)

solar_cell_solver(solar_cell_TMM_ARC, 'iv', user_options={'voltages': V, 'light_iv': True,
                                                  'wavelength': wavelengths, 'mpp': True,
                                                 'light_source': light_source,
                                                  'recalculate_absorption': True,
                                                  'optics_method': "TMM"})
plt.figure(2)


plt.plot(V, -solar_cell_TMM_ARC[1].iv(V)/10, 'b', label='Si')

plt.text(0.72, 34, 'Efficieny (%): ' + str(np.round(solar_cell_TMM_ARC.iv['Eta'] * 100, 1)))
plt.text(0.72, 32, 'FF (%): ' + str(np.round(solar_cell_TMM_ARC.iv['FF'] * 100, 1)))
plt.text(0.72, 30, r'V$_{oc}$ (V): ' + str(np.round(solar_cell_TMM_ARC.iv["Voc"], 2)))
plt.text(0.72, 28, r'I$_{sc}$ (A/m$^2$): ' + str(np.round(solar_cell_TMM_ARC.iv["Isc"]/10, 2)))

plt.legend()
plt.ylim(0, 40)
plt.xlim(0, 1)
plt.ylabel('Current (A/m$^2$)')
plt.xlabel('Voltage (V)')
plt.title("IV characteristics of Single Si cell")

plt.show()

