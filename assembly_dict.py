import gzip
import random

def lire_genome_fasta(fichier_fasta):
    """Lit un fichier FASTA (.fasta ou .fasta.gz) et retourne la séquence génomique complète."""
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
    """Extrait tous les k-mers de longueur k du génome ou des reads."""
    kmers = []
    for i in range(len(genome) - k + 1):
        kmers.append(genome[i:i+k])
    return kmers

def stock_kmers_dict(kmers, threshold=0):
    """Compte les k-mers et supprime ceux dont la fréquence est <= threshold."""
    kmer_dict = {}
    for kmer in kmers:
        kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
    # Supprime les k-mers peu fréquents
    kmer_dict = {k: v for k, v in kmer_dict.items() if v > threshold}
    return kmer_dict

def sucessor(kmer, kmer_dict):
    """Retourne le successeur du k-mer et le consomme du dictionnaire."""
    prefix = kmer[1:]
    for base in ['A', 'T', 'C', 'G']:
        candidate = prefix + base
        if candidate in kmer_dict:
            del kmer_dict[candidate]  # supprime le k-mer utilisé
            return candidate
    return None

def predecessors(kmer, kmer_dict):
    """Retourne le prédécesseur du k-mer et le consomme du dictionnaire."""
    suffix = kmer[:-1]
    for base in ['A', 'T', 'C', 'G']:
        candidate = base + suffix
        if candidate in kmer_dict:
            del kmer_dict[candidate]  # supprime le k-mer utilisé
            return candidate
    return None

def de_Bruijn_graph(kmer_dict):
    """Construit un contig à partir des k-mers restants et les consomme."""
    if not kmer_dict:
        return None

    first_kmer = random.choice(list(kmer_dict.keys()))
    del kmer_dict[first_kmer]  # consomme le k-mer de départ

    assembly = first_kmer

    # extension vers l’avant
    next_kmer = sucessor(first_kmer, kmer_dict)
    while next_kmer:
        assembly += next_kmer[-1]
        next_kmer = sucessor(next_kmer, kmer_dict)

    # extension vers l’arrière
    prev_kmer = predecessors(first_kmer, kmer_dict)
    while prev_kmer:
        assembly = prev_kmer[0] + assembly
        prev_kmer = predecessors(prev_kmer, kmer_dict)

    return assembly

def main():
    fichier_fasta = "./data/level00/reads.fasta.gz"
    genome = lire_genome_fasta(fichier_fasta)
    k = 25
    kmers = extraire_kmers_genome(genome, k)
    kmer_dict = stock_kmers_dict(kmers, threshold=4)

    contigs = []
    while kmer_dict:
        contig = de_Bruijn_graph(kmer_dict)
        if contig:
            contigs.append(contig)

    # Affiche les contigs et leurs longueurs
    for i, c in enumerate(contigs, 1):
        print(f"Contig {i} (longueur {len(c)}): {c[:60]}{'...' if len(c)>60 else ''}")

    print(f"\nNombre total de contigs : {len(contigs)}")
    contigs_kmer_only = [c for c in contigs if len(c) == k]
    print(f"Nombre de contigs de longueur k : {len(contigs_kmer_only)}")

if __name__ == "__main__":
    main()