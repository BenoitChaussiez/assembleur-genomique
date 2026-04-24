import hashlib
import time
import sys
import mmh3
import random
from bitarray import bitarray




def extraire_kmers_25(fichier_entree, fichier_sortie, k=30):

    with open(fichier_entree, "r") as f_in, open(fichier_sortie, "w") as f_out:
        
        sequence = ""

        for line in f_in:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    for i in range(len(sequence) - k + 1):
                        kmer = sequence[i:i+k]
                        f_out.write(kmer + "\n")
                sequence = ""
            else:
                sequence += line

        if sequence:
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                f_out.write(kmer + "\n")

class BloomFilter:
    def __init__(self, size, hash_count):
        self.size = size
        self.hash_count = hash_count
        self.bit_array = bitarray(size)
        self.bit_array.setall(0)

    def _hashes(self, item):
        """Génère les indices en utilisant hash64 pour une meilleure distribution"""
        hash1, hash2 = mmh3.hash64(item, seed=0)
        
        for i in range(self.hash_count):
            idx = (hash1 + i * hash2) % self.size
            yield idx

    def add(self, item):
        for idx in self._hashes(item):
            self.bit_array[idx] = 1
    
    def check(self, item):
        return all(self.bit_array[idx] for idx in self._hashes(item))
    

def compter_lignes(fichier):
    """Compte rapidement le nombre de lignes dans un fichier."""
    with open(fichier, 'r') as f:
        return sum(1 for _ in f)
    

def extension_gauche(max_iterations, kmer_initial, bf):
    """Extension vers la droite avec choix du meilleur chemin quand il y a des branches"""
    
    kmer_actuel = kmer_initial
    sequence_ajoutee_gauche = ""
    iteration = 0
    bf_verif_gauche=BloomFilter(size=800_000_000, hash_count=7)
    while iteration < max_iterations:
        base = kmer_actuel[:-1]
        candidats_gauche = {
            'A': 'A' + base ,
            'T': 'T' + base ,
            'C': 'C' + base,
            'G': 'G' + base 
        }
        
        candidats_valides = []
        for nuc, kmer_test in candidats_gauche.items():
            if bf.check(kmer_test) and not bf_verif_gauche.check(kmer_test):
                candidats_valides.append((nuc, kmer_test))
        
        if not candidats_valides:
            print(f'Arrêt gauche: aucun nouveau kmer à iteration {iteration}')
            break

        if len(candidats_valides) > 1:
            scores = []
            for nuc, kmer_test in candidats_valides:
                base_suivante = kmer_test[1:]
                score = sum(1 for n in 'ATCG' if bf.check(base_suivante + n))
                scores.append(score)
            
            meilleur_score = max(scores)
            meilleurs = [candidats_valides[i][0] for i, s in enumerate(scores) if s == meilleur_score]
            nuc_choisi = random.choice(meilleurs)
        else:
            nuc_choisi = candidats_valides[0][0]
        
        kmer_actuel = candidats_gauche[nuc_choisi]
        bf_verif_gauche.add(kmer_actuel)
        sequence_ajoutee_gauche = nuc_choisi + sequence_ajoutee_gauche
        iteration += 1
    
    return sequence_ajoutee_gauche
    
def extension_droite(max_iterations, kmer_initial, bf):
    """Extension vers la droite avec choix du meilleur chemin quand il y a des branches"""
    kmer_actuel = kmer_initial
    sequence_ajoutee_droite = ""
    iteration = 0
    bf_verif_droite=BloomFilter(size=800_000_000, hash_count=7)
    while iteration < max_iterations:
        base = kmer_actuel[1:]
        candidats_droite = {
            'A': base + 'A',
            'T': base + 'T',
            'C': base + 'C',
            'G': base + 'G'
        }
        
        candidats_valides = []
        for nuc, kmer_test in candidats_droite.items():
            if bf.check(kmer_test) and not bf_verif_droite.check(kmer_test):
                candidats_valides.append((nuc, kmer_test))
        
        if not candidats_valides:
            print(f'Arrêt droite: aucun nouveau kmer à iteration {iteration}')
            break

        if len(candidats_valides) > 1:
            scores = []
            for nuc, kmer_test in candidats_valides:
                base_suivante = kmer_test[1:]
                score = sum(1 for n in 'ATCG' if bf.check(base_suivante + n))
                scores.append(score)
            
            meilleur_score = max(scores)
            meilleurs = [candidats_valides[i][0] for i, s in enumerate(scores) if s == meilleur_score]
            nuc_choisi = random.choice(meilleurs)
        else:
            nuc_choisi = candidats_valides[0][0]
        
        kmer_actuel = candidats_droite[nuc_choisi]
        bf_verif_droite.add(kmer_actuel)
        sequence_ajoutee_droite += nuc_choisi
        iteration += 1
    
    return sequence_ajoutee_droite

def fusionner_contigs(contigs, k=30):
    """Fusionne les contigs qui se chevauchent"""
    
    contigs_tries = sorted(contigs, key=len, reverse=True)
    fusionnes = []
    
    for contig in contigs_tries:
        ajoute = False
        for i, existing in enumerate(fusionnes):
            if contig in existing:
                ajoute = True
                break
            
            if existing in contig:
                fusionnes[i] = contig
                ajoute = True
                break
            
            for overlap in range(min(k, len(contig)), 20, -1):
                if existing.endswith(contig[:overlap]):
                    fusionnes[i] = existing + contig[overlap:]
                    ajoute = True
                    break
                elif contig.endswith(existing[:overlap]):
                    fusionnes[i] = contig + existing[overlap:]
                    ajoute = True
                    break
            
            if ajoute:
                break
        
        if not ajoute:
            fusionnes.append(contig)
    
    return fusionnes
    
def construction_BLOOM(fichier1,nb_essai, taille_filtre=0, nb_hash=7):


    compte = compter_lignes(fichier1)
    taille_filtre = 1_400_000_000
    start = time.time()
    
    bf = BloomFilter(size=taille_filtre, hash_count=nb_hash)
    bf2 = BloomFilter(size=taille_filtre, hash_count=nb_hash)
    bf3 = BloomFilter(size=taille_filtre, hash_count=nb_hash)
    bf4 = BloomFilter(size=taille_filtre, hash_count=nb_hash)
    bf5 = BloomFilter(size=taille_filtre, hash_count=nb_hash)
    bf6 = BloomFilter(size=taille_filtre, hash_count=nb_hash)
    with open(fichier1, 'r') as f:
        tous_les_kmers = [line.strip() for line in f]
        print(f"Nombre total de k-mers: {len(tous_les_kmers)}")
        
        for kmer in tous_les_kmers:
            if bf.check(kmer):
                if bf2.check(kmer):
                    if bf3.check(kmer):
                        if bf4.check(kmer):
                            if bf5.check(kmer):
                                bf6.add(kmer)
                            else:
                                bf5.add(kmer)  
                        else:
                            bf4.add(kmer)  
                    else:
                        bf3.add(kmer)  
                else:
                    bf2.add(kmer)  
            else:
                bf.add(kmer)  

    del (bf, bf2, bf3, bf4, bf5)
    

    seeds_utilisees = set()  
    tous_les_contigs = []
    
    for essai in range(nb_essai):

        kmer_aleatoire = None
        max_tentatives = 100
        for _ in range(max_tentatives):
            candidat = random.choice(tous_les_kmers)
            if bf6.check(candidat) and candidat not in seeds_utilisees:
                kmer_aleatoire = candidat
                break
        
        if kmer_aleatoire is None:
            break
        
        seeds_utilisees.add(kmer_aleatoire) 
        

        sequence_droite = extension_droite(compte, kmer_aleatoire, bf6)
        sequence_gauche = extension_gauche(compte, kmer_aleatoire, bf6)
        
        contig = sequence_gauche + kmer_aleatoire + sequence_droite
        
        if len(contig) >= 5000:
            tous_les_contigs.append(contig)

    

    contigs = tous_les_contigs
    for passe in range(3):
        contigs = fusionner_contigs(contigs, k=30)

    

    with open('fichier_sortie.fasta', 'w') as f_out:
        for i, contig in enumerate(contigs):
            f_out.write(f">contig_{i}_k=30\n")
            for j in range(0, len(contig), 80):
                f_out.write(contig[j:j+80] + "\n")
    
    taille_totale = sum(len(c) for c in contigs)
    print(f"Nombre de contigs: {len(contigs)}")
    print(f"Taille totale: {taille_totale} pb")
    
    if len(contigs) > 0:
        print(f"Contig max: {max(len(c) for c in contigs)} pb")
    
    end = time.time()
    print(f"TEMPS construction contig: {end-start:.2f} secondes")
    
    return contigs