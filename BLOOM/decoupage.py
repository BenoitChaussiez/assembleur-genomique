import sys

def extraire_kmers(fichier_entree, fichier_sortie, k=30):

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


# ----- EXECUTION -----
if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python script.py input.fa output_kmers.txt")
        sys.exit(1)

    extraire_kmers(sys.argv[1], sys.argv[2])

    print("Extraction terminée.")