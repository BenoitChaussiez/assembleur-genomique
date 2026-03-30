import hashlib
import time
import sys

class BloomFilter:
    def __init__(self, size, hash_count):
        self.size = size
        self.hash_count = hash_count
        # Utiliser bytearray au lieu de list pour économiser la mémoire
        self.bit_array = bytearray(size)

    def _hashes(self, item):
        # Utiliser SHA256 pour meilleure distribution
        h = int(hashlib.sha256(item.encode()).hexdigest(), 16)
        for i in range(self.hash_count):
            yield (h + i * 0x9e3779b97f4a7c15) % self.size

    def add(self, item):
        for idx in self._hashes(item):
            self.bit_array[idx] = 1

    def check(self, item):
        return all(self.bit_array[idx] for idx in self._hashes(item))


def compter_lignes(fichier):
    """Compte rapidement le nombre de lignes dans un fichier."""
    with open(fichier, 'r') as f:
        return sum(1 for _ in f)


def comparer_fichiers(fichier1, fichier2, taille_filtre=500_000_000, nb_hash=7):
    """
    Compare deux fichiers de kmers.
    """
    print("="*60)
    print("COMPARAISON DE FICHIERS DE KMERS")
    print("="*60)
    
    print("\nAnalyse des fichiers...")
    n1 = compter_lignes(fichier1)
    n2 = compter_lignes(fichier2)
    print(f"Fichier 1 : {fichier1} ({n1:,} kmers)")
    print(f"Fichier 2 : {fichier2} ({n2:,} kmers)")
    
    print(f"\nInitialisation du filtre de Bloom (taille: {taille_filtre:,}, hash: {nb_hash})...")
    bf = BloomFilter(size=taille_filtre, hash_count=nb_hash)
    table = set()
    
    print(f"\nChargement de {fichier1}...")
    start_load = time.time()
    
    with open(fichier1, 'r') as f:
        for i, line in enumerate(f):
            kmer = line.strip()
            if kmer:
                bf.add(kmer)
                table.add(kmer)
            
            if i > 0 and i % 2_000_000 == 0:
                elapsed = time.time() - start_load
                print(f"  {i:,} kmers chargés ({elapsed:.1f}s)")
    
    end_load = time.time()
    print(f" {len(table):,} kmers uniques chargés en {end_load - start_load:.2f}s")
    
    print(f"\nRecherche dans {fichier2}...")
    commons = []
    start_search = time.time()
    
    with open(fichier2, 'r') as f:
        for i, line in enumerate(f):
            kmer = line.strip()
            if not kmer:
                continue
            
            if bf.check(kmer) and kmer in table:
                commons.append(kmer)
            
            if i > 0 and i % 2_000_000 == 0:
                elapsed = time.time() - start_search
                print(f"  {i:,} kmers analysés... {len(commons):,} trouvés ({elapsed:.1f}s)")
    
    end_search = time.time()
    
    print("\n" + "="*60)
    print("RÉSULTATS FINAUX")
    print("="*60)
    print(f" {fichier1} : {len(table):,} kmers uniques")
    print(f" {fichier2} : {n2:,} kmers")
    print(f" Kmers communs : {len(commons):,}")
    
    if len(commons) > 0:
        pourcentage = (len(commons) / min(len(table), n2)) * 100
        print(f" Pourcentage : {pourcentage:.2f}% du plus petit fichier")
    
    print(f"\n Temps de chargement : {end_load - start_load:.2f}s")
    print(f" Temps de recherche : {end_search - start_search:.2f}s")
    print(f"  Temps total : {end_search - start_load:.2f}s")
    
    if len(commons) > 0:
        print(f"\n 10 premiers kmers communs :")
        for i, kmer in enumerate(commons[:10]):
            print(f"   {i+1}. {kmer}")
    
    return commons

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python kmer.py <fichier_kmers1> <fichier_kmers2>")
        print("Exemple: python kmer.py fichier1.fa fichier2.fa")
        sys.exit(1)
    
    try:
        communs = comparer_fichiers(sys.argv[1], sys.argv[2])
    except FileNotFoundError as e:
        print(f" Erreur: Fichier non trouvé - {e}")
    except MemoryError:
        print("Erreur: Mémoire insuffisante. Essayez de réduire la taille du filtre.")
    except Exception as e:
        print(f" Erreur inattendue: {e}")