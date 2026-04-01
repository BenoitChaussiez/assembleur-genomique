import gzip
import random

def lire_genome_fasta(fichier_fasta):
    """
    Lit un fichier FASTA (.fasta ou .fasta.gz) et retourne la séquence génomique complète.
    """
    genome = []

    if fichier_fasta.endswith(".gz"):
        f = gzip.open(fichier_fasta, 'rt')
    else:
        f = open(fichier_fasta, 'r')

    with f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            genome.append(line)

    return ''.join(genome)


def extraire_kmers_genome(genome, k):
    kmers = []
    for i in range(0, len(genome) - k + 1):
        kmers.append(genome[i:i+k])
    return kmers


def stock_kmers_dict(kmers):
    kmer_dict = {}
    for kmer in kmers:
        if kmer in kmer_dict:
            kmer_dict[kmer] += 1
        else:
            kmer_dict[kmer] = 1
    return kmer_dict


def sucessor(kmer, kmer_dict):
    prefix = kmer[1:]
    
    for base in ['A', 'T', 'C', 'G']:
        candidate = prefix + base
        if candidate in kmer_dict:
            kmer_dict[candidate] -= 1  # consomme
            return candidate
    return None

def predecessors(kmer, kmer_dict):
    suffix = kmer[:-1]
    
    for base in ['A', 'T', 'C', 'G']:
        candidate = base + suffix
        if candidate in kmer_dict:
            kmer_dict[candidate] -= 1  # consomme
            return candidate
    return None
    
def de_Bruijn_graph(kmer_dict):
    first_kmer = random.choice(list(kmer_dict.keys()))
    print(f"Premier k-mer choisi : {first_kmer}")
    assembly = first_kmer
    next = sucessor(first_kmer, kmer_dict)
    while next is not None:
        assembly += next[-1]
        next = sucessor(next, kmer_dict)
    previous = predecessors(first_kmer, kmer_dict)
    while previous is not None:
        assembly = previous[0] + assembly
        previous = predecessors(previous, kmer_dict)
    return assembly

def main():
    fichier_fasta = "./data/level00/reads.fasta.gz"
    genome = lire_genome_fasta(fichier_fasta)
    k = 25
    kmers = extraire_kmers_genome(genome, k)
    kmer_dict = stock_kmers_dict(kmers)
    print(kmer_dict)
    assembly = de_Bruijn_graph(kmer_dict)
    print("Assemblage de Bruijn :")
    print(assembly)
    print(f"Longueur de l'assemblage : {len(assembly)}")
    reference = lire_genome_fasta("./data/level00/reference.fasta")
    if assembly in reference:
        print("L'assemblage est présent dans le génome.")
    else:
        print("L'assemblage n'est pas présent dans le génome.")

if __name__ == "__main__":
    main()