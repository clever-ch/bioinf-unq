from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna
#Dependencias para GUI
import tkinter as tk
from tkinter import filedialog
#Dependencia para abortar/terminar programa
import sys

#Constantes
CODONES_STOP = ('UAA', 'UAG', 'UGA')
proteinaLetras = {'E','F','I','J','L','O','P','Q','Z'}


proteina = "ACCGGUUUAACUUGE"
ARN = "ACCGGUUUAACUUGA"
ADN = "ACCGGTTTAACTTGA"

def isProtein(fasta):
    letrasUnicas = set(fasta)
    i = proteinaLetras.intersection(letrasUnicas)
    return (i != set())

def isADN(fasta):
    return 'T' in fasta

def isARN(fasta):
    return 'U' in fasta

# print(isProtein(proteina))
# print(isADN(ADN))
# print(isARN(ARN))
#
# messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG", generic_rna)
# print(messenger_rna.translate())
# coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", generic_dna)
# print(coding_dna.translate())


#Elimina los guines de la secuencia si es que los hubiese
def desalinear(secuencia):
    secuencia._set_seq(secuencia.seq.ungap('-')) #TODO medio raro porque tecnicamente un seqRecord es inmutable

def divisiblePor3(secuencia):
    print(f'Longitud de la secuencia: {len(secuencia)}')
    if len(secuencia) % 3 != 0:
        print("La secuencia no es divisible por 3")
        sys.exit(0)

def hayCodonesStop(secuencia):
    print(secuencia)
    #isStopCodon = lambda x: any(x[i:i + 3] in CODONES_STOP for i in range(0, len(x), 3))

    codonesEncontrados = 0
    for i in range(0, len(secuencia), 3):
        if secuencia.seq[i:i+3] in CODONES_STOP:
            codonesEncontrados += 1

    if codonesEncontrados > 1:
        print(f'Hay multiples({codonesEncontrados}) codones de stop en la secuencia')
    elif codonesEncontrados == 1 and secuencia.seq.endswith(CODONES_STOP):
        print('Hay un codon de stop pero esta al final asi que no pasa nada')
    else:
        print('No hay codones de stop')



if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()

    fastaFileDirection = filedialog.askopenfilename(title = 'Seleccione archivo FASTA')

    secuencias = list(SeqIO.parse(fastaFileDirection, "fasta"))

    if len(secuencias) > 0:
        print("Se encontraron varias secuencias en el archivo.")
        index = 0
        for seq_record in secuencias:
            print("ID: " + str(index))
            index += 1
            print(seq_record.id)
            print(repr(seq_record.seq))
            print(f'longitud: {len(seq_record)}')
        secuenciaNoSeleccionado = True
        options = list(map(lambda o: str(o), range(len(secuencias))))
        while secuenciaNoSeleccionado :
            idSeleccionado = input('Indique el ID de la secuencia que quiere utilizar:')
            if idSeleccionado not in options:
                print("No hay cadena con ese ID")
            else:
                secuenciaNoSeleccionado = False
                print("Cadena seleccionada:")
                seq_selected = secuencias[int(idSeleccionado)]
    else:
        seq_selected = secuencias[0]

    print(seq_selected.id)
    desalinear(seq_selected)
    quiereMarcoLectura = input("¿Desea definir un marco de lectura? Y/N")
    if quiereMarcoLectura == 'Y':
        marcoLectura = (int(input("Inicio: ")), int(input("Final: ")))
        seq_selected = seq_selected[marcoLectura[0]:marcoLectura[1]]
        print("Seccion seleccionada:")
        print(seq_selected.seq)
    #Chequeos de la secuencia:
    divisiblePor3(seq_selected)
    hayCodonesStop(seq_selected)
    print(seq_selected)
