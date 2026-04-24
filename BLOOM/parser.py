import argparse

def parse():
    parser = argparse.ArgumentParser(description='Bloom filter pour faire le génome')
    parser.add_argument('-r','--reads', required=True, help='Fichier FASTQ.gz des reads')
    parser.add_argument('-o','--out', required=True, help='Fichier de sortie')
    parser.add_argument('-t','--essai', required=True, help="Nombre d'essai")
    parser.add_argument('-k','--kmer', required=True, help='Taille de k-mer')
    return parser.parse_args()