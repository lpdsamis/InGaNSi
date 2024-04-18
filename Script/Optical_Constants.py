#!/usr/bin/env python
# coding: utf-8

import solcore
from solcore.absorption_calculator.nk_db import download_db, search_db, create_nk_txt
from solcore.absorption_calculator import calculate_rat, OptiStack
from solcore.material_system import create_new_material
from solcore import material
from solcore import si
from solcore.structure import Layer
import os
import numpy as np


import matplotlib.pyplot as plt
# pip install SciencePlots
import scienceplots

plt.style.use(['science','ieee'])

wl = si(np.arange(100, 1500, 5), 'nm')

GaN=material(name='365', nk_db=True)()

# create_new_material('InGaN', '../data/InGaN_n.txt', '../data/InGaN_k.txt', '../data/InGaN_params.txt')
# create_new_material('GaInN', '../data/InGaN_n.txt', '../data/InGaN_k.txt')

# GaInN = material('GaInN')(In=0.54)
InGaN = material('InGaN')()
GaInP = material("GaInP")(In=0.5)
Si=material('Si')()

param_n = dict(xlabel='Wavelength (nm)', ylabel=r' n ')
param_k = dict(xlabel='Wavelength (nm)', ylabel=r' k ')
with plt.style.context(['science', 'ieee']):
    fig, ax = plt.subplots()

    ax.plot(wl * 1e9, GaN.n(wl), label='GaN')
    ax.plot(wl * 1e9, InGaN.n(wl), label='InGaN')
    ax.plot(wl * 1e9, GaInP.n(wl), label='GaInP')
    ax.plot(wl * 1e9, Si.n(wl), label='Si')

    ax.legend()
    ax.autoscale(tight=True)
    ax.set(**param_n)
    fig.savefig('../figures/fig_n.jpg', dpi=300)
    plt.close()
with plt.style.context(['science', 'ieee']):
    fig, ax = plt.subplots()

    ax.plot(wl * 1e9, GaN.k(wl), label='GaN')
    ax.plot(wl * 1e9, InGaN.k(wl), label='InGaN')
    ax.plot(wl * 1e9, GaInP.k(wl), label='GaInP')
    ax.plot(wl * 1e9, Si.k(wl), label='Si')

    ax.legend()
    ax.autoscale(tight=True)
    ax.set(**param_k)
    fig.savefig('../figures/fig_k.jpg', dpi=300)
    plt.close()