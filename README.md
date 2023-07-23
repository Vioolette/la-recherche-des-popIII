# La recherche des galaxies popIII

Le fichier `spectre` contient le code permettant de tracer les spectres 2D et 1D de l'ensemble des données d'un fichier.

Le fichier `read_yggdrasil` contient le code permettant de créer des nouveaux fichiers fits en séparant les modèles Yggdrasil.

Le fichier `redshift` contient le code permettant de calculer le redshift et t de l'ensemble des objets d'un dossier. Selon les spectres et leur origine, certaines informations peuvent changer. Ainsi, le nombre de variable des données du fichier va changer. Il faut modifier les lignes 43, 45, 153, 155, 263 et 265 en concervant autant de lettres qu'il y a de varibles dans le fichier fits observé. Il a été écrit en collaboration avec mon encadrant. Le dossier sur lequel on lance le code doit être composer des fichiers fits du spectre 1D des objets célestes.
