import gzip


class KmerIndex:
    def __init__(self, k):
        self.k = k
        self.index = {}  # dictionnaire k-mer -> compte

    def insert_sequence(self, seq):
        """Insère tous les k-mers d'une séquence"""
        for i in range(len(seq) - self.k + 1):
            kmer = seq[i:i+self.k]
            if kmer in self.index:
                self.index[kmer] += 1
            else:
                self.index[kmer] = 1

    def contains(self, kmer):
        """Vérifie si le k-mer est présent"""
        return kmer in self.index


def load_fasta_gz(file_path):
    sequences = []
    seq = ""
    # ouverture en mode texte avec gzip
    with gzip.open(file_path, 'rt') as f:  # 'rt' = read text
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line
        if seq:
            sequences.append(seq)
    return sequences

def filter_kmers_by_abundance(k_index, threshold):
    """
    Supprime du dictionnaire tous les k-mers dont l'abondance < threshold.
    Modifie k_index.index directement.
    """
    to_delete = [kmer for kmer, count in k_index.index.items() if count < threshold]
    for kmer in to_delete:
        del k_index.index[kmer]

def successors(kmer, k_index):
    """Retourne la liste des k-mers successeurs existants"""
    bases = ['A', 'C', 'G', 'T']
    succ = []
    for b in bases:
        next_kmer = kmer[1:] + b
        if k_index.contains(next_kmer):
            succ.append(next_kmer)
    return succ

def predecessors(kmer, k_index):
    """Retourne la liste des k-mers prédécesseurs existants"""
    bases = ['A', 'C', 'G', 'T']
    pred = []
    for b in bases:
        prev_kmer = b + kmer[:-1]
        if k_index.contains(prev_kmer):
            pred.append(prev_kmer)
    return pred

def build_contig(start_kmer, k_index):
    contig = start_kmer
    visited = set([start_kmer])

    # Étendre vers l’avant
    current = start_kmer
    while True:
        succ = [s for s in successors(current, k_index) if s not in visited]
        if not succ:
            break
        current = succ[0]  # on prend le premier successeur disponible
        contig += current[-1]  # ajouter la dernière base seulement
        visited.add(current)

    # Étendre vers l’arrière
    current = start_kmer
    prefix = ""
    while True:
        pred = [p for p in predecessors(current, k_index) if p not in visited]
        if not pred:
            break
        current = pred[0]  # on prend le premier prédécesseur
        prefix = current[0] + prefix  # ajouter la première base seulement
        visited.add(current)

    contig = prefix + contig
    return contig

# Paramètres
k = 25
fasta_file = "./data/level02/reads.fasta.gz"

# Chargement des séquences
sequences = load_fasta_gz(fasta_file)

# Création de l’index
k_index = KmerIndex(k)
for seq in sequences:
    k_index.insert_sequence(seq)

# Vérification d’un k-mer
print(k_index.contains(sequences[0][:k]))  # True
print(k_index.contains("A"*k))            # False probable

threshold = 3
filter_kmers_by_abundance(k_index, threshold)
# Assemblage d’un contig
import random
start_kmer = random.choice(list(k_index.index.keys()))
contig = build_contig(start_kmer, k_index)
print("Taille du contig assemblé :", len(contig))
with open("contig_assembled_lvl2.fasta", "w") as f:
    f.write(">Contig_1\n")
    f.write(contig + "\n")
