with open('level01/reference.fasta', 'r') as g:
    contenu = g.read()
    nombre_caracteres = len(contenu)
    print(f"Nombre de caractères: {nombre_caracteres}")