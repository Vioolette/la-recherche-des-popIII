from read import read_yggdrasil
from astropy.io import fits
from scipy import stats

import numpy as np
import matplotlib.pyplot as plt
import os


# Do we need to print stuff? Remplacer par True selon ce qu'on veut afficher
verbose_False   = False
verbose_True    = False
plot_test_False = False
plot_test_True  = False
plot_loglog     = False

def check_interp(x_before, y_before, x_after, y_after, plot_loglog):
# We check here whether the interpolation is OK
    plt.plot(x_before, y_before, 'ro', label='obs before')
    plt.plot(x_after, y_after, 'bx', label='obs after')
    if plot_loglog:
        plt.xscale('log')
        plt.yscale('log')
    plt.legend()
    plt.show()

def plot_line(x, y, label, plot_loglog):
    plt.plot(x, y, 'go-', label=label, ms=0, linewidth=1)
    if plot_loglog:
        plt.xscale('log')
        plt.yscale('log')
    plt.legend()
    plt.show()

# Le fichier doit être un spectre 1D
def redshift_dusty(fits_spec_filename):

    with fits.open(fits_spec_filename) as hdu_spec:
        data_spec = hdu_spec[1].data

    mask_obs = np.isfinite(data_spec['FLUX'])

    tab_wave = [x for x, y, z, a, b ,c, d, e, f, g, h, i, j, k, l, m, n, o in data_spec]
    spec_wave = np.asarray(tab_wave)[mask_obs]*1e4  # [um] -> [A]
    tab_flux = [y for x, y, z, a, b ,c, d, e, f, g, h, i, j, k, l, m, n, o in data_spec]
    spec_flux = np.asarray(tab_flux)[mask_obs]*1e3  # [Jy] -> [mJy]

    n_data_spec = spec_wave.size

    min_spec_wave = np.min(spec_wave)
    max_spec_wave = np.max(spec_wave)
    delta_spec_lambda_mean = (max_spec_wave - min_spec_wave) / (n_data_spec-1)
    delta_spec_lambda_min = min(el2 - el1 for el1, el2 in zip(spec_wave[:-1], spec_wave[1:]))
    n_spec_resampled = int((max_spec_wave - min_spec_wave) / delta_spec_lambda_min) + 1
# `Because this is a test, we have no uncertainties on the flux densities.
# So, let's assume we have 10% of uncertainties on the flux densities
    spec_error = data_spec['FLUX_ERROR']*1e3     # [mJy]
# We have the observed arrays (wave, flux) from about 0.6 to 6 microns.

    if verbose_False:
        print('Observed data: ', spec_wave, spec_flux, spec_error)
    if verbose_True:
        print('There are ', n_data_spec, 'datapoints from ', min_spec_wave, ' to ', max_spec_wave, ' microns')
        print('With a mean Delta_lambda = ', np.round(delta_spec_lambda_mean, 2), 'A', ' and a min Delta_lambda = ', delta_spec_lambda_min, 'A')

# Récuperation of models / templates
    fits_modelfilename = 'model_wmetals_dusty.fits'
    with fits.open(fits_modelfilename) as hdu_model:
        data_model = hdu_model[1].data
        n_data_model = data_model.size

    tab_wave = [x for x, y in data_model]
    model_wave = np.asarray(tab_wave)*1e4  # [um] -> [A]
    tab_flux = [y for x, y in data_model]
    model_flux = np.asarray(tab_flux)*1e3       # [Jy] -> [mJy]
    min_model_wave = np.min(model_wave)
    max_model_wave = np.max(model_wave)
    if verbose_False:
        print('Models: ', model_wave, model_flux)
    if verbose_True:
        print('There are ', n_data_model, 'datapoints from ', min_model_wave, ' to ', max_model_wave, ' microns')
# We have the models (wave, flux) rest-frame from below the Lyman break to 6 microns.

# Now, we proceed to a resampling in two steps:
#  1. we resample the observations to a constant step: delta_spec_lambda_mean

#Rééchantillonage du spectre observé#delta_lamda = 0.0193
    lambda_newobs  = np.linspace(min_spec_wave, max_spec_wave, n_spec_resampled)
    spec_newobs = np.interp(lambda_newobs, spec_wave, spec_flux)
# Vérification graphique que l'interpolation est ok
    if plot_test_False:
        check_interp(spec_wave, spec_flux, lambda_newobs, spec_newobs, plot_loglog)
# We have the resampled observed arrays (wave, flux) from about 0.6 to 6 microns.
# resampled to the minimum Delta_lambda

#  2. We resample the models to the same constant step: delta_speTruec_lambda_mean
# We redshift the model with n_redshifts redshift values in the range z = [z_min, z_max]
    z_min =  0.000
    z_max = 20.000
    n_redshifts = 201
    correlation = np.empty(n_redshifts)
    correlation[:] = np.NaN  #  set all correlation coefficients to NaN
    redshift_range = np.linspace(z_min, z_max, n_redshifts)

    for i_z_new, z_new in enumerate(redshift_range):
        if verbose_False:
            print ('redshift: ', np.round(z_new, 3))#, end='\r')
        lambda_break_obs = 912. * (1+z_new)
    # We redshift the wavelengths of the models for each of the explored redshifts
        lambda_newmodel = model_wave * (1+z_new)
    # We create a mask for the wavelengths outside the observed wavelength range
        mask_spec_wave = np.logical_and(lambda_newmodel>=min_spec_wave, lambda_newmodel<=max_spec_wave, np.isfinite(model_flux))
    # We resample the models to the resampled observation wavelength grid (lambda_newobs)
        spec_newmodel = np.interp(lambda_newobs, lambda_newmodel[mask_spec_wave], model_flux[mask_spec_wave])
# Vérification graphique que l'interpolation est ok
        if plot_test_False:
            check_interp(model_wave, model_flux, lambda_newobs, spec_newmodel, plot_loglog)
# We plot model and data
        if plot_test_False:
            check_interp(lambda_newobs[lambda_newobs>lambda_break_obs], spec_newmodel[lambda_newobs>lambda_break_obs]/np.mean(spec_newmodel[lambda_newobs>lambda_break_obs])*np.mean(spec_newobs[lambda_newobs>lambda_break_obs]), lambda_newobs[lambda_newobs>lambda_break_obs], spec_newobs[lambda_newobs>lambda_break_obs], plot_loglog)
# Calcul de la corrélation entre les deux
# We cross-correlate the observed and the modelled spectra to derive the redshift
        if len(spec_newmodel)>15:
            correlation[i_z_new] = stats.pearsonr(spec_newobs[lambda_newobs>lambda_break_obs], spec_newmodel[lambda_newobs>lambda_break_obs]/np.mean(spec_newmodel[lambda_newobs>lambda_break_obs])*np.mean(spec_newobs[lambda_newobs>lambda_break_obs]))[0]
        else:
            correlation[i_z_new] = np.nan
        if verbose_True:
            print ('Correlation for redshift: ', np.round(z_new, 3), '=', correlation[i_z_new])#, end='\r')

    best_ind = correlation.argmax()
    best_corr = correlation[best_ind]
    best_z = redshift_range[best_ind]
    info_best_z = 'Best z = ' + str(np.round(best_z, 3)) + ' and corr. coeff = ' + str(np.round(best_corr, 3))
    print('\n===>', info_best_z)

    if verbose_False:
        print(redshift_range, correlation)
    if plot_test_True:
        plt.scatter(1.+best_z, best_corr, c='red', s = 20, label=info_best_z)
        plot_line(1.+redshift_range, correlation, 'Corr. Coeff. (z)', plot_loglog)

    # plt.savefig("redshift_" + fits_spec_filename[:-5] + ".png")
    return(best_z, best_corr)


def redshift_hiz(fits_spec_filename):

    with fits.open(fits_spec_filename) as hdu_spec:
        data_spec = hdu_spec[1].data

    mask_obs = np.isfinite(data_spec['FLUX'])

    tab_wave = [x for x, y, z, a, b ,c, d, e, f, g, h, i, j, k, l, m, n, o in data_spec]
    spec_wave = np.asarray(tab_wave)[mask_obs]*1e4  # [um] -> [A]
    tab_flux = [y for x, y, z, a, b ,c, d, e, f, g, h, i, j, k, l, m, n, o in data_spec]
    spec_flux = np.asarray(tab_flux)[mask_obs]*1e3  # [Jy] -> [mJy]

    n_data_spec = spec_wave.size

    min_spec_wave = np.min(spec_wave)
    max_spec_wave = np.max(spec_wave)
    delta_spec_lambda_mean = (max_spec_wave - min_spec_wave) / (n_data_spec-1)
    delta_spec_lambda_min = min(el2 - el1 for el1, el2 in zip(spec_wave[:-1], spec_wave[1:]))
    n_spec_resampled = int((max_spec_wave - min_spec_wave) / delta_spec_lambda_min) + 1
# `Because this is a test, we have no uncertainties on the flux densities.
# So, let's assume we have 10% of uncertainties on the flux densities
    spec_error = data_spec['FLUX_ERROR']*1e3     # [mJy]
# We have the observed arrays (wave, flux) from about 0.6 to 6 microns.

    if verbose_False:
        print('Observed data: ', spec_wave, spec_flux, spec_error)
    if verbose_True:
        print('There are ', n_data_spec, 'datapoints from ', min_spec_wave, ' to ', max_spec_wave, ' microns')
        print('With a mean Delta_lambda = ', np.round(delta_spec_lambda_mean, 2), 'A', ' and a min Delta_lambda = ', delta_spec_lambda_min, 'A')

# Récuperation of models / templates
    fits_modelfilename = 'model_wmetals_hiz.fits'
    with fits.open(fits_modelfilename) as hdu_model:
        data_model = hdu_model[1].data
        n_data_model = data_model.size

    tab_wave = [x for x, y in data_model]
    model_wave = np.asarray(tab_wave)*1e4  # [um] -> [A]
    tab_flux = [y for x, y in data_model]
    model_flux = np.asarray(tab_flux)*1e3       # [Jy] -> [mJy]
    min_model_wave = np.min(model_wave)
    max_model_wave = np.max(model_wave)
    if verbose_False:
        print('Models: ', model_wave, model_flux)
    if verbose_True:
        print('There are ', n_data_model, 'datapoints from ', min_model_wave, ' to ', max_model_wave, ' microns')
# We have the models (wave, flux) rest-frame from below the Lyman break to 6 microns.

# Now, we proceed to a resampling in two steps:
#  1. we resample the observations to a constant step: delta_spec_lambda_mean

#Rééchantillonage du spectre observé#delta_lamda = 0.0193
    lambda_newobs  = np.linspace(min_spec_wave, max_spec_wave, n_spec_resampled)
    spec_newobs = np.interp(lambda_newobs, spec_wave, spec_flux)
# Vérification graphique que l'interpolation est ok
    if plot_test_False:
        check_interp(spec_wave, spec_flux, lambda_newobs, spec_newobs, plot_loglog)
# We have the resampled observed arrays (wave, flux) from about 0.6 to 6 microns.
# resampled to the minimum Delta_lambda

#  2. We resample the models to the same constant step: delta_spec_lambda_mean
# We redshift the model with n_redshifts redshift values in the range z = [z_min, z_max]
    z_min =  0.000
    z_max = 20.000
    n_redshifts = 201
    correlation = np.empty(n_redshifts)
    correlation[:] = np.NaN  #  set all correlation coefficients to NaN
    redshift_range = np.linspace(z_min, z_max, n_redshifts)

    for i_z_new, z_new in enumerate(redshift_range):
        if verbose_False:
            print ('redshift: ', np.round(z_new, 3))#, end='\r')
        lambda_break_obs = 912. * (1+z_new)
    # We redshift the wavelengths of the models for each of the explored redshifts
        lambda_newmodel = model_wave * (1+z_new)
    # We create a mask for the wavelengths outside the observed wavelength range

        mask_spec_wave = np.logical_and(lambda_newmodel>=min_spec_wave, lambda_newmodel<=max_spec_wave)
    # We resample the models to the resampled observation wavelength grid (lambda_newobs)
        spec_newmodel = np.interp(lambda_newobs, lambda_newmodel[mask_spec_wave], model_flux[mask_spec_wave])
# Vérification graphique que l'interpolation est ok
        if plot_test_False:
            check_interp(model_wave, model_flux, lambda_newobs, spec_newmodel, plot_loglog)
# We plot model and data
        if plot_test_False:
            check_interp(lambda_newobs[lambda_newobs>lambda_break_obs], spec_newmodel[lambda_newobs>lambda_break_obs]/np.mean(spec_newmodel[lambda_newobs>lambda_break_obs])*np.mean(spec_newobs[lambda_newobs>lambda_break_obs]), lambda_newobs[lambda_newobs>lambda_break_obs], spec_newobs[lambda_newobs>lambda_break_obs], plot_loglog)
# Calcul de la corrélation entre les deux
# We cross-correlate the observed and the modelled spectra to derive the redshift
        if len(spec_newmodel)>15:
            correlation[i_z_new] = stats.pearsonr(spec_newobs[lambda_newobs>lambda_break_obs], spec_newmodel[lambda_newobs>lambda_break_obs]/np.mean(spec_newmodel[lambda_newobs>lambda_break_obs])*np.mean(spec_newobs[lambda_newobs>lambda_break_obs]))[0]
        else:
            correlation[i_z_new] = np.nan
        if verbose_True:
            print ('Correlation for redshift: ', np.round(z_new, 3), '=', correlation[i_z_new])#, end='\r')

    best_ind = correlation.argmax()
    best_corr = correlation[best_ind]
    best_z = redshift_range[best_ind]
    info_best_z = 'Best z = ' + str(np.round(best_z, 3)) + ' and corr. coeff = ' + str(np.round(best_corr, 3))
    print('\n===>', info_best_z)

    if verbose_False:
        print(redshift_range, correlation)
    if plot_test_True:
        plt.scatter(1.+best_z, best_corr, c='red', s = 20, label=info_best_z)
        plot_line(1.+redshift_range, correlation, 'Corr. Coeff. (z)', plot_loglog)
    # plt.savefig("redshift_" + fits_spec_filename[:-5] + ".png")
    return(best_z, best_corr)


def redshift_pop3(fits_spec_filename, fits_modelfilename):

    with fits.open(fits_spec_filename) as hdu_spec:
        data_spec = hdu_spec[1].data

    mask_obs = np.isfinite(data_spec['FLUX'])

    tab_wave = [x for x, y, z, a, b ,c, d, e, f, g, h, i, j, k, l, m, n, o in data_spec]
    spec_wave = np.asarray(tab_wave)[mask_obs]*1e4  # [um] -> [A]
    tab_flux = [y for x, y, z, a, b ,c, d, e, f, g, h, i, j, k, l, m, n, o in data_spec]
    spec_flux = np.asarray(tab_flux)[mask_obs]*1e3  # [Jy] -> [mJy]

    n_data_spec = spec_wave.size

    min_spec_wave = np.min(spec_wave)
    max_spec_wave = np.max(spec_wave)
    delta_spec_lambda_mean = (max_spec_wave - min_spec_wave) / (n_data_spec-1)
    delta_spec_lambda_min = min(el2 - el1 for el1, el2 in zip(spec_wave[:-1], spec_wave[1:]))
    n_spec_resampled = int((max_spec_wave - min_spec_wave) / delta_spec_lambda_min) + 1
# `Because this is a test, we have no uncertainties on the flux densities.
# So, let's assume we have 10% of uncertainties on the flux densities
    spec_error = data_spec['FLUX_ERROR']*1e3     # [mJy]
# We have the observed arrays (wave, flux) from about 0.6 to 6 microns.

    if verbose_False:
        print('Observed data: ', spec_wave, spec_flux, spec_error)
    if verbose_True:
        print('There are ', n_data_spec, 'datapoints from ', min_spec_wave, ' to ', max_spec_wave, ' microns')
        print('With a mean Delta_lambda = ', np.round(delta_spec_lambda_mean, 2), 'A', ' and a min Delta_lambda = ', delta_spec_lambda_min, 'A')

# Récuperation of models / templates
    with fits.open(fits_modelfilename) as hdu_model:
        data_model = hdu_model[1].data
        n_data_model = data_model.size

    tab_wave = [x for x, y in data_model]
    model_wave = np.asarray(tab_wave)
    tab_flux = [y for x, y in data_model]
    for i in range(len(tab_flux)):               # [erg /s / A2] -> [erg / Hz / s / cm2]
        tab_flux[i] = (tab_flux[i] * 299792458 ) / ((tab_wave[i] * 1e-6)** 2)
    model_flux = np.asarray(tab_flux)*1e26       # [erg / Hz / s / cm2] -> [mJy]
    min_model_wave = np.min(model_wave)
    max_model_wave = np.max(model_wave)
    if verbose_False:
        print('Models: ', model_wave, model_flux)
    if verbose_True:
        print('There are ', n_data_model, 'datapoints from ', min_model_wave, ' to ', max_model_wave, ' microns')
# We have the models (wave, flux) rest-frame from below the Lyman break to 6 microns.

# Now, we proceed to a resampling in two steps:
#  1. we resample the observations to a constant step: delta_spec_lambda_mean

#Rééchantillonage du spectre observé#delta_lamda = 0.0193
    lambda_newobs  = np.linspace(min_spec_wave, max_spec_wave, n_spec_resampled)
    spec_newobs = np.interp(lambda_newobs, spec_wave, spec_flux)
# Vérification graphique que l'interpolation est ok
    if plot_test_False:
        check_interp(spec_wave, spec_flux, lambda_newobs, spec_newobs, plot_loglog)
# We have the resampled observed arrays (wave, flux) from about 0.6 to 6 microns.
# resampled to the minimum Delta_lambda

#  2. We resample the models to the same constant step: delta_spec_lambda_mean
# We redshift the model with n_redshifts redshift values in the range z = [z_min, z_max]
    z_min =  0.000
    z_max = 20.000
    n_redshifts = 201
    correlation = np.empty(n_redshifts)
    correlation[:] = np.NaN  #  set all correlation coefficients to NaN
    redshift_range = np.linspace(z_min, z_max, n_redshifts)

    for i_z_new, z_new in enumerate(redshift_range):
        if verbose_False:
            print ('redshift: ', np.round(z_new, 3))#, end='\r')
        lambda_break_obs = 912. * (1+z_new)
    # We redshift the wavelengths of the models for each of the explored redshifts
        lambda_newmodel = model_wave * (1+z_new) * 1e4
    # We create a mask for the wavelengths outside the observed wavelength range
        mask_spec_wave = np.logical_and(lambda_newmodel>=min_spec_wave, lambda_newmodel<=max_spec_wave)
    # We resample the models to the resampled observation wavelength grid (lambda_newobs)
        spec_newmodel = np.interp(lambda_newobs, lambda_newmodel[mask_spec_wave], model_flux[mask_spec_wave])
# Vérification graphique que l'interpolation est ok
        if plot_test_False:
            check_interp(model_wave, model_flux, lambda_newobs, spec_newmodel, plot_loglog)
# We plot model and data
        if plot_test_False:
            check_interp(lambda_newobs[lambda_newobs>lambda_break_obs], spec_newmodel[lambda_newobs>lambda_break_obs]/np.mean(spec_newmodel[lambda_newobs>lambda_break_obs])*np.mean(spec_newobs[lambda_newobs>lambda_break_obs]), lambda_newobs[lambda_newobs>lambda_break_obs], spec_newobs[lambda_newobs>lambda_break_obs], plot_loglog)
# Calcul de la corrélation entre les deux
# We cross-correlate the observed and the modelled spectra to derive the redshift
        if len(spec_newmodel)>15:
            correlation[i_z_new] = stats.pearsonr(spec_newobs[lambda_newobs>lambda_break_obs], spec_newmodel[lambda_newobs>lambda_break_obs]/np.mean(spec_newmodel[lambda_newobs>lambda_break_obs])*np.mean(spec_newobs[lambda_newobs>lambda_break_obs]))[0]
        else:
            correlation[i_z_new] = np.nan
        if verbose_True:
            print ('Correlation for redshift: ', np.round(z_new, 3), '=', correlation[i_z_new])#, end='\r')

    best_ind = correlation.argmax()
    best_corr = correlation[best_ind]
    best_z = redshift_range[best_ind]
    info_best_z = 'Best z = ' + str(np.round(best_z, 3)) + ' and corr. coeff = ' + str(np.round(best_corr, 3))
    print('\n===>', info_best_z)

    if verbose_False:
        print(redshift_range, correlation)
    if plot_test_True:
        plt.scatter(1.+best_z, best_corr, c='red', s = 20, label=info_best_z)
        plot_line(1.+redshift_range, correlation, 'Corr. Coeff. (z)', plot_loglog)

    # plt.savefig("redshift_" + fits_spec_filename[:-5] + ".png")
    return(best_z, best_corr)

red = {}
# Nom du dossier à modifier en fonction 
dossier_images = "S3_out_P5_v0.3_x1D"
for image in os.listdir(dossier_images):
    coeff = {}
    mod = {}
    print("./"+dossier_images+"/"+image)
    red_dusty = redshift_dusty("./"+dossier_images+"/"+image)
    coeff[red_dusty[0]]=red_dusty[1]
    mod[red_dusty[0]] = "dusty"

    red_hiz = redshift_hiz("./"+dossier_images+"/"+image)
    if red_hiz[0] in coeff:
        if max(coeff[red_hiz[0]], red_hiz[1]) == coeff[red_hiz[0]]:
            c = 1
        else:
            c = 0
        coeff[red_hiz[0]]=max(coeff[red_hiz[0]], red_hiz[1])
        if c == 0:
            mod[red_hiz[0]] = 'hiz'
    else:
        coeff[red_hiz[0]]=red_hiz[1]
        mod[red_hiz[0]]='hiz'

    for nom_modele in os.listdir("./Yggdrasil/POPIII_models"):
        red_pop3 = redshift_pop3("./"+dossier_images+"/"+image, "./Yggdrasil/POPIII_models/" + nom_modele)
        if red_pop3[0] in coeff:
            if max(coeff[red_pop3[0]], red_pop3[1]) == coeff[red_pop3[0]]:
                c = 1
            else:
                c = 0
            coeff[red_pop3[0]]=max(coeff[red_pop3[0]], red_pop3[1])
            if c == 0:
                mod[red_pop3[0]] = nom_modele
        else:
            coeff[red_pop3[0]]=red_pop3[1]
            mod[red_pop3[0]] = nom_modele
    
    red_max = max(coeff, key=coeff.get)
    red[image] = (red_max, coeff[red_max], mod[red_max])


# Création d'un fichier fits pour enregistrer les données

obj_data = []
for cle in red:
    obj_data.append(cle)

red_data = []
coef_corr_data = []
mod_data = []

for data in red.values():
    red_data.append(data[0])
    coef_corr_data.append(data[1])
    mod_data.append(data[2])


# Création des colonnes FITS
col1 = fits.Column(name='Objet', format='A100', array=obj_data)
col2 = fits.Column(name='redshift', format='D', array=red_data)
col3 = fits.Column(name='coefficient de correlation', format='D', array=coef_corr_data)
col4 = fits.Column(name='modele associe', format='A70', array=mod_data)

# Création du tableau FITS
cols = fits.ColDefs([col1, col2, col3, col4])
tbhdu = fits.BinTableHDU.from_columns(cols)

# Écriture du fichier FITS
tbhdu.writeto(dossier_images + '_res_redshifts.fits', overwrite=True)

# Afficher le contenu du fichier avec les résultats
im = fits.open("S3_out_P5_v0.3_x1D_res_redshifts.fits")
print(im[1].data)
