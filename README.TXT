Partie Paul Dictionnaire 

Partie Benoit Filtre de Bloom

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
