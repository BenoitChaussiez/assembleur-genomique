def trie(kmers: list[str]) -> list[str]:
    """
    Renvoie une nouvelle liste de k-mers triés par ordre alphabétique.
    """
    return sorted(kmers)


def dichotomique(list1: list[str], list2: list[str]) -> list[str]:
    """
    Recherche les k-mers de list2 dans list1 (après tri de list1)
    et renvoie la liste des k-mers trouvés.
    """

    def recherche_dichotomique(kmers: list[str], cible: str) -> int:
        gauche = 0
        droite = len(kmers) - 1

        while gauche <= droite:
            milieu = (gauche + droite) // 2

            if kmers[milieu] == cible:
                return milieu
            elif kmers[milieu] < cible:
                gauche = milieu + 1
            else:
                droite = milieu - 1

        return -1

    # 1. on trie la première liste
    list1_triee = trie(list1)

    # 2. on cherche chaque k-mer de list2 dans list1
    resultats = []
    for kmer in list2:
        if recherche_dichotomique(list1_triee, kmer) != -1:
            resultats.append(kmer)

    return resultats

list1 = [
    "TTGACGATCGATCGATCGATCGATC",
    "AACGATCGATCGATCGATCGATCGA",
    "GGTACGATCGATCGATCGATCGATC",
    "CCGATCGATCGATCGATCGATCGAT"]

list2 = [
    "GGTACGATCGATCGATCGATCGATC",
    "AAAAAAAAAAAAAAAAAAAAAAAAA",
    "CCGATCGATCGATCGATCGATCGAT"]

print(dichotomique(list1, list2))