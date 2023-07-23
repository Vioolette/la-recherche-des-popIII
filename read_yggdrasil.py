import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from os import walk
from astropy.table import Table

import re
plot_True     = True
plot_False    = False
verbose_True  = True
verbose_False = False

def read_yggdrasil(filename):
    """Read uggdrasil data file."""
    total_mass = None
    ages, gas_mass, stellar_mass, wl_grid, spectra = [], [], [], [], []
    current_wl, current_spectrum = [], []

    with open(filename, 'r') as infile:
        for line in infile.readlines():
            line = line.strip()
            if line.startswith('Mass available'):
                total_mass = float(line.split(':')[1])
            elif line.startswith('Age'):
                ages.append(float(line.split(':')[1]))
            elif line.startswith(('Gas')):
                gas_mass.append(float(line.split(':')[1])/total_mass)
            elif line.startswith(('Stellar')):
                stellar_mass.append(float(line.split(':')[1])/total_mass)
            elif re.match('^[0-9]', line):
                current_wl.append(float(line.split()[0]))
                current_spectrum.append(float(line.split()[1]) / total_mass)
            elif line == '':
                if current_wl:
                    if not wl_grid:
                        wl_grid = current_wl
                    if wl_grid and current_wl != wl_grid:
                        raise ValueError('The wavelength grid changed!')
                    spectra.append(current_spectrum)
                    current_wl, current_spectrum = [], []

    return ages, gas_mass, stellar_mass, wl_grid, spectra
    

# We read all filenames in the Yggdrasil directory
filenames = next(walk('./Yggdrasil'), (None, None, []))[2]  # [] if no file
if verbose_False:
    print('The filenames are:', filenames)

# We loop over the file names, and read it.
for filename in filenames:
    filename = './Yggdrasil/' + filename
    if filename.endswith('Spectra'):

        ages, gas_mass, stellar_mass, wl_grid, spectra = read_yggdrasil(filename)
        wl_grid = [float(wave)*1e-4 for wave in wl_grid]
        n_pixels = len(wl_grid)
        n_spectra = np.shape(spectra)[0]

# We read each of the individual models, plot them, and save them in separate files.
        data = np.zeros((n_pixels, 2))
        for i_spec in range(n_spectra):
            label = 'model_age:'+str(ages[i_spec]) + '_m_gas:' + str(gas_mass[i_spec]) + '_m_star:' + str(stellar_mass[i_spec])
            filename_out = './Yggdrasil/POPIII_models/Yggdrasil'+'_age'+str(ages[i_spec]) + '_m_gas' + str(gas_mass[i_spec]) + '_m_star' + str(stellar_mass[i_spec]) + '.fits'

# We save the models to a new file
            hdulist = fits.BinTableHDU.from_columns(  
               [fits.Column(name='lambda', format='E', array=wl_grid),
                fits.Column(name='Fnu', format='E', array=spectra[:][i_spec])])
            hdulist.writeto(filename_out, overwrite=True)

# We plot the models
            if plot_False:
                plt.plot(wl_grid, spectra[:][i_spec], label=label)
                plt.xscale('log')
                plt.yscale('log')
                plt.legend()
        if plot_False:
            plt.show()

