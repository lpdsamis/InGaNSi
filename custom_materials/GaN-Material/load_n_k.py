import numpy as np
import pandas as pd


data_txt = np.loadtxt('/home/amis/.solcore/custom_materials/GaN-Material/GaN-zb.refraction')

data_txtDF = pd.DataFrame(data_txt)

df2= data_txtDF.rename({0: 'Energy [eV]', 1: 'ref_ind_xx', 2:'extinct_xx'}, axis='columns')

df2=df2[(df2['Energy [eV]']<5)]

df2["Energy [eV]"]=1.239837538868846e-06/df2["Energy [eV]"]

df3=df2.rename({'Energy [eV]': 'wavelength [m]'}, axis='columns')

df3=df3.sort_values(by=['wavelength [m]'])

GaN_n_xx=df3[['wavelength [m]', 'ref_ind_xx']]

GaN_k_xx=df3[['wavelength [m]', 'extinct_xx']]


GaN_n_xx.to_csv(r'/home/amis/.solcore/custom_materials/GaN-Material/GaN_n.txt', header=None, index=None, sep=' ')

GaN_k_xx.to_csv(r'/home/amis/.solcore/custom_materials/GaN-Material/GaN_k.txt', header=None, index=None, sep=' ')

  


    
