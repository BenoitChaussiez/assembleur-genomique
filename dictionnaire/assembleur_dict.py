import gzip
import argparse
import os
import matplotlib.pyplot as plt
from collections import Counter
import random
import sys
import time

# Augmenter la limite de récursion pour les graphes profonds
sys.setrecursionlimit(5000)
random.seed(42)

# ==========================================
# 1. Outils d'indexation et de charhement
# ==========================================

class KmerIndex:
    """
    Structure de données pour indexer et compter les k-mers d'une séquence.
    
    >>> idx = KmerIndex(3)
    >>> idx.insert_sequence("ATCGT")
    >>> idx.index == {'ATC': 1, 'TCG': 1, 'CGT': 1}
    True
    >>> idx.contains("ATC")
    True
    >>> idx.contains("AAA")
    False
    """
    def __init__(self, k):
        self.k = k
        self.index = {}

    def insert_sequence(self, seq):
        """Découpe une séquence en k-mers et met à jour leurs occurrences dans l'index."""
        for i in range(len(seq) - self.k + 1):
            kmer = seq[i:i+self.k]
            self.index[kmer] = self.index.get(kmer, 0) + 1

    def contains(self, kmer):
        """Vérifie si un k-mer spécifique existe dans l'index."""
        return kmer in self.index

def load_fasta_gz(file_path):
    """
    Charge les séquences d'un fichier FASTA compressé (.gz).
    
    Args:
        file_path (str): Le chemin vers le fichier FASTA.gz.
    
    Returns:
        list: Une liste contenant toutes les séquences ADN sous forme de chaînes de caractères.
    """
    sequences = []
    seq = ""
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq: sequences.append(seq)
                seq = ""
            else:
                seq += line
        if seq: sequences.append(seq)
    return sequences

def generer_spectre_kmers(k_index, k_size, max_freq=150, output_file="spectre_kmers.png"):
    """
    Génère et sauvegarde un graphique représentant le spectre d'abondance des k-mers.
    """
    abondances = list(k_index.index.values())
    compte_abondances = Counter(abondances)

    x = range(1, max_freq + 1)
    y = [compte_abondances.get(i, 0) for i in x]

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, marker='.', linestyle='-')
    plt.yscale('log')
    plt.title(f'Spectre des k-mers (k={k_size})')
    plt.xlabel('Abondance du k-mer (Fréquence d\'apparition)')
    plt.ylabel('Nombre de k-mers distincts (Échelle Log)')
    plt.grid(True, which="both", ls="--", alpha=0.5)

    plt.savefig(output_file)
    print(f"Graphique sauvegardé : {output_file}")

def filter_kmers_by_abundance(k_index, threshold):
    """
    Supprime de l'index les k-mers dont l'occurrence est strictement inférieure au seuil.
    
    >>> idx = KmerIndex(3)
    >>> idx.index = {'AAA': 5, 'CCC': 2, 'GGG': 1}
    >>> filter_kmers_by_abundance(idx, 3)
    Filtrage terminé : 1 k-mers restants.
    >>> idx.index == {'AAA': 5}
    True
    """
    to_delete = [kmer for kmer, count in k_index.index.items() if count < threshold]
    for kmer in to_delete:
        del k_index.index[kmer]
    print(f"Filtrage terminé : {len(k_index.index)} k-mers restants.")

# ==========================================
# 2. Implémentation du graph de Bruijn
# ==========================================

def successors(kmer, k_index):
    """
    Trouve tous les k-mers successeurs valides présents dans le dictionnaire.
    Un successeur est formé en supprimant la première base et en ajoutant A, C, G ou T à la fin.

    >>> idx = KmerIndex(3)
    >>> idx.index = {'TCG': 1, 'TCA': 1, 'TCC': 1, 'AAA': 1}
    >>> sorted(successors("ATC", idx))
    ['TCA', 'TCC', 'TCG']
    """
    suffix = kmer[1:]
    return [suffix + b for b in ['A', 'C', 'G', 'T'] if suffix + b in k_index.index]

def predecessors(kmer, k_index):
    """
    Trouve tous les k-mers prédécesseurs valides présents dans le dictionnaire.
    Un prédécesseur est formé en supprimant la dernière base et en ajoutant A, C, G ou T au début.

    >>> idx = KmerIndex(3)
    >>> idx.index = {'AAT': 1, 'CAT': 1, 'GAT': 1, 'CCC': 1}
    >>> sorted(predecessors("ATC", idx))
    ['AAT', 'CAT', 'GAT']
    """
    prefix = kmer[:-1]
    return [b + prefix for b in ['A', 'C', 'G', 'T'] if b + prefix in k_index.index]

def evaluate_path(current_node, k_index, global_visited, current_depth, max_depth):
    """
    Évalue récursivement la profondeur maximale atteignable depuis un nœud,
    afin d'aider l'algorithme à résoudre les bifurcations (Lookahead).
    """
    if current_depth >= max_depth:
        return current_depth
    
    neighbors = successors(current_node, k_index)
    valid_neighbors = [n for n in neighbors if n not in global_visited]
    
    if not valid_neighbors:
        return current_depth
    
    # On limite à l'exploration des 2 premiers voisins pour éviter l'explosion combinatoire
    return max(evaluate_path(n, k_index, global_visited, current_depth + 1, max_depth) 
               for n in valid_neighbors[:2])

def choose_best_next(next_nodes, k_index, visited, max_depth):
    """
    Sélectionne le meilleur nœud suivant parmi plusieurs choix possibles en utilisant
    l'évaluation de chemin (Lookahead). Privilégie un chemin aléatoire parmi ceux
    qui survivent jusqu'à la profondeur maximale.
    """
    if not next_nodes: return None
    if len(next_nodes) == 1: return next_nodes[0]
    
    path_scores = []
    for node in next_nodes:
        depth = evaluate_path(node, k_index, visited, 1, max_depth)
        path_scores.append((depth, node))
    
    path_scores.sort(key=lambda x: x[0], reverse=True)
    best_depth = path_scores[0][0]
    
    valid_paths = [node for depth, node in path_scores if depth >= max_depth]
    
    if valid_paths:
        return random.choice(valid_paths)
    return path_scores[0][1] if best_depth >= 10 else None

def build_contig(start_kmer, k_index, max_lookahead):
    """
    Construit un contig (chaîne ADN) en naviguant dans le graphe de De Bruijn
    à partir d'un k-mer initial.
    
    Returns:
        tuple: (Séquence du contig généré, Set des k-mers visités lors de la construction)
    """
    contig_kmers = [start_kmer]
    visited_in_contig = {start_kmer}
    
    current = start_kmer
    while True:
        neighbors = successors(current, k_index)
        next_nodes = [n for n in neighbors if n not in visited_in_contig]
        
        chosen = choose_best_next(next_nodes, k_index, visited_in_contig, max_lookahead)
        if chosen is None: break
            
        visited_in_contig.add(chosen)
        contig_kmers.append(chosen)
        current = chosen
        
    sequence = contig_kmers[0] + "".join([k[-1] for k in contig_kmers[1:]])
    return sequence, visited_in_contig

def get_n50(contigs):
    """
    Calcule les métriques N50 et L50 d'un assemblage.
    Le N50 est la longueur du plus petit contig de l'ensemble des plus grands contigs
    dont la somme couvre au moins 50% de la taille totale de l'assemblage.

    >>> contigs_1 = ["A"*100, "C"*50, "G"*30, "T"*20] # Total = 200, Moitié = 100
    >>> get_n50(contigs_1) # Le contig de 100 atteint la moitié à lui seul (L50 = 1)
    (100, 1)
    >>> contigs_2 = ["A"*40, "C"*30, "G"*20, "T"*10] # Total = 100, Moitié = 50
    >>> get_n50(contigs_2) # 40 + 30 = 70 (>= 50), donc N50 = 30, L50 = 2
    (30, 2)
    >>> get_n50([]) # Test liste vide
    (0, 0)
    """
    if not contigs:
        return 0, 0
        
    lengths = sorted([len(c) for c in contigs], reverse=True)
    total = sum(lengths)
    cumul = 0
    for i, l in enumerate(lengths):
        cumul += l
        if cumul >= total / 2: return l, i + 1
    return 0, 0

# ==========================================
# 3. Fonction main + parser pour préciser le niveau et les paramètres
# ==========================================

def main():
    """Point d'entrée principal du programme."""
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Assembleur de De Bruijn pour le projet de M1.")
    parser.add_argument('-l', '--level', type=int, required=True, help="Niveau du dataset à assembler (0 à 10).")
    parser.add_argument('-k', '--ksize', type=int, default=31, help="Taille des k-mers (défaut: 31).")
    parser.add_argument('-t', '--threshold', type=int, default=5, help="Seuil d'abondance pour filtrer le bruit (défaut: 5).")
    parser.add_argument('-d', '--lookahead', type=int, default=50, help="Profondeur de recherche (lookahead) (défaut: 50).")
    parser.add_argument('-s', '--target_size', type=int, default=1500000, help="Taille cible max de l'assemblage en bp (défaut: 1500000).")
    parser.add_argument('-m', '--min_len', type=int, default=5000, help="Longueur minimum des contigs à sauvegarder (défaut: 5000).")
    
    args = parser.parse_args()

    # Formatage des chemins d'accès des fichiers
    level_str = f"level{args.level:02d}"
    input_file = f"./data/{level_str}/reads.fasta.gz"
    output_fasta = f"FINAL_ASSEMBLY_{level_str.upper()}_k{args.ksize}_t{args.threshold}.fasta"
    output_png = f"spectre_{level_str}_k{args.ksize}.png"

    if not os.path.exists(input_file):
        print(f"❌ Erreur : Le fichier {input_file} n'existe pas.")
        return

    print(f"\nDÉMARRAGE ASSEMBLAGE - NIVEAU {args.level}")
    print(f"Paramètres : k={args.ksize} | seuil={args.threshold} | profondeur={args.lookahead} | filtre_taille={args.min_len}bp\n")

    #Indexation
    print(f"Chargement des données depuis {input_file}...")
    sequences = load_fasta_gz(input_file)
    k_index = KmerIndex(args.ksize)
    for seq in sequences: 
        k_index.insert_sequence(seq)

    #Diagnostic & Filtrage
    generer_spectre_kmers(k_index, args.ksize, max_freq=100, output_file=output_png)
    filter_kmers_by_abundance(k_index, args.threshold)

    #Assemblage
    global_visited = set()
    all_contigs = []
    sorted_kmers = sorted(k_index.index.keys(), key=lambda x: k_index.index[x], reverse=True)

    print("\nDébut de l'assemblage avec déduplication...")
    for start_kmer in sorted_kmers:
        if start_kmer in global_visited:
            continue
        
        seq, kmers_utilises = build_contig(start_kmer, k_index, args.lookahead)
        
        if len(seq) >= 500: # Seuil minimal interne avant le filtrage final
            all_contigs.append(seq)
            if len(all_contigs) % 5 == 0:
                print(f"  -> {len(all_contigs)} contigs trouvés. Somme actuelle : {sum(len(c) for c in all_contigs)} bp")

        global_visited.update(kmers_utilises)
        for kmer in kmers_utilises:
            for n in successors(kmer, k_index): global_visited.add(n)
            for p in predecessors(kmer, k_index): global_visited.add(p)
            
        if sum(len(c) for c in all_contigs) >= args.target_size:
            print("\n✅ Taille cible atteinte, arrêt de l'assemblage.")
            break

    #Filtrage Final & Statistiques
    all_contigs.sort(key=len, reverse=True)
    final_selection = [c for c in all_contigs if len(c) >= args.min_len]

    n50, l50 = get_n50(final_selection)

    print(f"\n--- RAPPORT FINAL NIVEAU {args.level} ---")
    print(f"Contigs retenus (>{args.min_len}bp) : {len(final_selection)}")
    print(f"Nmax : {len(final_selection[0]) if final_selection else 0} bp")
    print(f"Somme totale : {sum(len(c) for c in final_selection)} bp")
    print(f"N50 : {n50} bp | L50 : {l50}")

    #Sauvegarde
    with open(output_fasta, "w") as f:
        for i, seq in enumerate(final_selection):
            f.write(f">Contig_{i+1}_len_{len(seq)}\n")
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + "\n")

    print(f"\nFichier '{output_fasta}' généré avec succès.")

    end_time = time.time()
    elapsed_time = end_time - start_time
    
    #Formatage propre pour l'affichage (secondes ou minutes)
    if elapsed_time < 60:
        print(f"Temps d'exécution total : {elapsed_time:.2f} secondes\n")
    else:
        minutes = int(elapsed_time // 60)
        secondes = elapsed_time % 60
        print(f"Temps d'exécution total : {minutes} minutes et {secondes:.2f} secondes\n")

if __name__ == "__main__":
    main()