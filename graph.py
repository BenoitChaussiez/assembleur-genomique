import argparse
def parser():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="Integrer 2 fichier")
    parser.add_argument('-f1', '--fichier1', required=True,dest='fichier1')
    parser.add_argument('-f2', '--fichier2',required=True,dest='fichier2')
    parser.add_argument('-o', '--output',required=True,dest='output')
    return parser.parse_args()

def intersection_to_file(fileA, fileB, output_edges):
    count = 0

    with open(fileA) as fa, open(fileB) as fb, open(output_edges, "w") as out:
        a = fa.readline().strip()
        b = fb.readline().strip()

        while a and b:
            if a == b:
                prefix = a[:-1]
                suffix = a[1:]
                out.write(f"{prefix}\t{suffix}\n")
                count += 1

                a = fa.readline().strip()
                b = fb.readline().strip()

            elif a < b:
                a = fa.readline().strip()
            else:
                b = fb.readline().strip()

    print(f"Intersection size: {count}")

def main():
    args = parser()
    intersection_to_file(args.fichier1,args.fichier2,args.output)
    
if __name__ == "__main__":
    main()