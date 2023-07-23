import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter

def spectre(nom_image, chemin_dossier):
    image = fits.open(chemin_dossier + "/" + nom_image)

    filtrage_x = 3.0
    filtrage_y = 3.0
    min_val = -1.0e-13
    max_val = 2.0e-13
    focus_y_min = 7
    focus_y_max = 13
    aspect_graphique = 13

    spectre_2D_data0 = image[1].data

    spectre_2D_data = gaussian_filter(spectre_2D_data0, [filtrage_x, filtrage_y])

    spectres, axes = plt.subplots(2)
    spectres.suptitle("Spectre 2D et 1D (du haut vers le bas)")

    axes[0].imshow(spectre_2D_data, cmap='viridis', vmin=min_val, vmax=max_val, aspect=aspect_graphique)
    liste_x = np.round(np.linspace(1.0, 5.0, 7), 2) # liste de 7 nonmbre entre 1 et 5 arrondis 2 chiffres après la virgule
    nb_x_2D = spectre_2D_data.shape[1]
    axes[0].set_xticks(np.linspace(0, nb_x_2D, 7))
    axes[0].set_xticklabels(liste_x)

    axes[0].axhline(y=focus_y_min, color='red')
    axes[0].axhline(y=focus_y_max, color='red')

    spectre_1D_data = np.sum(spectre_2D_data[focus_y_min:focus_y_max, :], axis=0)
    if "140" in nom_image:
        x_1D = np.linspace(0.97, 1.84, spectre_1D_data.shape[0])
    elif "235" in nom_image:
        x_1D = np.linspace(1.66, 3.07, spectre_1D_data.shape[0])
    else:
        x_1D = np.linspace(2.87, 5.1, spectre_1D_data.shape[0])
    axes[1].plot(x_1D, spectre_1D_data)

    spectres.savefig("spectres_" + nom_image[:-5] + ".png", format="png")
    # plt.show()
    return([x_1D, spectre_1D_data])


def spectre_global(nom_image1, chemin1, nom_image2, chemin2, nom_image3, chemin3):

    data_140M = spectre(nom_image1, chemin1)
    data_235M = spectre(nom_image2, chemin2)
    data_395M = spectre(nom_image3, chemin3)

    x_140M = data_140M[0]
    spectre_140M = data_140M[1]

    x_235M = data_235M[0]
    spectre_235M = data_235M[1]

    x_395M = data_395M[0]
    spectre_395M = data_395M[1]

# interpolation du 140M

    delta_lambda = 0.0038

    lambda_min_140M = 0.97
    lambda_max_140M = 1.84
    nb_points_140M = int((lambda_max_140M - lambda_min_140M) // delta_lambda)

    nouveaux_x_140M = np.linspace(lambda_min_140M, lambda_max_140M, nb_points_140M)
    spectre_140M = np.interp(nouveaux_x_140M, x_140M, spectre_140M)

# interpolation du 235M

    lambda_min_235M = 1.66
    lambda_max_235M = 3.07
    nb_points_235M = int((lambda_max_235M - lambda_min_235M) // delta_lambda)

    nouveaux_x_235M = np.linspace(lambda_min_235M, lambda_max_235M, nb_points_235M)
    spectre_235M = np.interp(nouveaux_x_235M, x_235M, spectre_235M)

# interpolation du 235M

    lambda_min_395M = 2.87
    lambda_max_395M = 5.10
    nb_points_395M = int((lambda_max_395M - lambda_min_395M) // delta_lambda)

    nouveaux_x_395M = np.linspace(lambda_min_395M, lambda_max_395M, nb_points_395M)
    spectre_395M = np.interp(nouveaux_x_395M, x_395M, spectre_395M)

# calcul du facteur entre 140M et 235M
    nb_recouvrement_140M = int(len(spectre_140M) * 0.21)
    nb_recouvrement_235M_1 = int(len(spectre_235M) * 0.13)

    zone_recouvrement_140M = spectre_140M[len(spectre_140M) - nb_recouvrement_140M:]
    moyenne_140M = sum(zone_recouvrement_140M) / len(zone_recouvrement_140M)

    zone_recouvrement_235M_1 = spectre_235M[:nb_recouvrement_235M_1]
    moyenne_235M_1 = sum(zone_recouvrement_235M_1)/len(zone_recouvrement_235M_1)

    facteur_235M = moyenne_235M_1 / moyenne_140M

# mise à niveau du 245M

    for el in spectre_235M:
        el *= facteur_235M


# calcul du facteur entre 235M et 395M
    nb_recouvrement_235M_2 = int(len(spectre_235M) * 0.14)
    nb_recouvrement_395M = int(len(spectre_395M) * 0.09)

    zone_recouvrement_235M_2 = spectre_235M[len(spectre_235M) - nb_recouvrement_235M_2:]
    moyenne_235M_2 = sum(zone_recouvrement_235M_2) / len(zone_recouvrement_235M_2)

    zone_recouvrement_395M = spectre_395M[:nb_recouvrement_395M]
    moyenne_395M = sum(zone_recouvrement_395M)/len(zone_recouvrement_395M)

    facteur_395M = moyenne_395M / moyenne_235M_2

# mise à niveau du 395M

    for el in spectre_395M:
        el *= facteur_395M

# tracé du spectre global

    zone_recouvrement_140M_g = int(len(spectre_140M)*0.21)
    zone_recouvrement_235M_g1 = int(len(spectre_235M)*0.13)
    zone_recouvrement_235M_g2 = int(len(spectre_235M)*0.14)
    zone_recouvrement_395M_g = int(len(spectre_395M)*0.09)

    ref_140M = len(spectre_140M) - zone_recouvrement_140M_g
    ref_235M = len(spectre_235M) - zone_recouvrement_235M_g2

    spectre_recouvrement_140M_235M = []
    for i in range(min(zone_recouvrement_140M_g, zone_recouvrement_235M_g1)):
        spectre_recouvrement_140M_235M.append((spectre_140M[ref_140M + i] + spectre_235M[i])/2)

    if zone_recouvrement_140M_g < zone_recouvrement_235M_g1:
        for el in spectre_235M[zone_recouvrement_140M_g:zone_recouvrement_235M_g1]:
            spectre_recouvrement_140M_235M.append(el)

    spectre_recouvrement_235M_395M = []
    for i in range(min(zone_recouvrement_235M_g2, zone_recouvrement_395M_g)):
        spectre_recouvrement_235M_395M.append((spectre_235M[ref_235M + i] + spectre_395M[i])/2)

    if zone_recouvrement_235M_g2 < zone_recouvrement_395M_g:
        for el in spectre_235M[zone_recouvrement_235M_g2:zone_recouvrement_395M_g]:
            spectre_recouvrement_235M_395M.append(el)

    spectre_data = list(spectre_140M[:ref_140M]) + spectre_recouvrement_140M_235M + list(spectre_235M[zone_recouvrement_235M_g1: ref_235M]) + spectre_recouvrement_235M_395M + list(spectre_395M[zone_recouvrement_395M_g:])
    x_spectre = np.linspace(0.97, 5.10, len(spectre_data))

# pour faire apparaître le spectre, décommenter
    plt.plot(x_spectre, spectre_data)
    plt.title("spectre global")

    plt.show()
    plt.savefig("spectre_" + chemin1[(len(chemin1) - 23) :] + ".png", format="png")
    return(x_spectre, spectre_data)
