from parser import parse
from BLOOM import construction_BLOOM
from decoupage import extraire_kmers


if __name__ == "__main__":


    args = parse()
    
    extraire_kmers(args.reads,args.out,int(args.kmer))

    construction_BLOOM(args.out,int(args.essai))