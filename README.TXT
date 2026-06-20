Assembleur génomique : Dictionnaire vs Filtre de Bloom

Ce projet, réalisé dans le cadre du Master MISO à l'Université de Lille, compare deux approches pour l'assemblage génomique à partir de reads :

1. Approche par dictionnaire : utilisation d'une table de hachage (dictionnaire Python) pour stocker et manipuler les k-mers.
2. Approche par filtre de Bloom : utilisation d'une cascade de filtres de Bloom pour un stockage mémoire optimisé.

L'objectif est de comprendre les avantages et les limites de chaque méthode en termes de performances, d'utilisation mémoire et de qualité d'assemblage.


Partie Paul Lemmonier Dictionnaire 

Description:
Ce programme permet de faire l'assemblage d'un génome à partir de reads en utilisant un dictionnaire comme structure de données

il est nécessaire d'avoir installer le dossier data contenant les différents niveaux pour pouvoir faires les assemblages

Fichier : assembleur_dict.py

Paramètres a renseigner : 
-l = Niveau du dataset à assembler (0 à 10) Obligatoire.
-k = Taille des k-mers (défaut: 31).
-t = Seuil d'abondance pour filtrer le bruit (défaut: 5).
-d = Profondeur de recherche (lookahead) (défaut: 50).
-s = Taille cible max de l'assemblage en bp (défaut: 1500000).
-m = Longueur minimum des contigs à sauvegarder (défaut: 5000).

Exemple de ligne de commande : 
python3 assembler_dict.py -l 6 -k 31 -t 6 -d 50 -s 2000000 -m 15000

Partie Benoit Chaussiez Filtre de Bloom

Prerequis: 
mmh3
bitarray

Description:
Ce programme permet de faire l'assemblage d'un génome à partir de reads, ce dernier utilise un Bloom filter comme base afin d'assembler les k-mers pour en faire un génome

Fonctionnalités:
  Extraire tous les k-mers des reads

  Filtrer les k-mers erronés via une cascade de 6 filtres de Bloom

  Assembler le génome par extension à partir de seeds (k-mers fiables)

  Fusionner les contigs obtenus

Fichier utile au Bloom filter
parser.py # Parsing des différen paramètres souhaité
BLOOM.py #Fonction de Bloom Filter , Extension gauche et droite , fusion des contigs et construction du génome
decoupage.py #Fonction de decoupe des fichiers de reads pour le transformer en fichier avec des k-mers
main.py #Fonction d'affichage qui renvoie les resultats dans un fichier

Fichier du /Bloom:
decoupage.py
mains.py
BLOOM.py 
parser.py

Parser :
-r fichier de reads
-o fichier de k-mer crée
-k taille de k-mer souhaité
-t nombre d'essai/seed voulu

Exemple de ligne de commande : 
python3 mains.py -r reads.fasta -o kmer_30.fasta -k 30 -t 50

Projet réalisé dans le cadre du Master MISO - Méthodes Informatiques et Statistiques pour les Systèmes Omics - Université de Lille, Année 2025-2026.