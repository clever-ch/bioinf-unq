from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna
from Bio.Blast import NCBIWWW
from Bio.PDB import PDBList

#Dependencias para GUI
import tkinter as tk
from tkinter import filedialog, messagebox

#Dependencia para abortar/terminar programa
import sys

#BLAST
from Bio.Blast import NCBIXML

#Dependencias para mutacion manual
import os, subprocess, tempfile

#Dependencia para modelado
from modelado import *

#Pymol
import pymol
from pymol import cmd
import __main__

#Constantes
CODONS_STOP = ('UAA', 'UAG', 'UGA')
CODON_START = ['AUG']
proteinaLetras = {'E','F','I','J','L','O','P','Q','Z'}
ARNLetras = {'A','C','G','U'}


def isProtein(fasta):
    letrasUnicas = set(fasta)
    i = proteinaLetras.intersection(letrasUnicas)
    return (i != set())

def isADN(fasta):
    return 'T' in fasta

def isARN(fasta):
    return 'U' in fasta



#Elimina los guines de la secuencia si es que los hubiese
def desalinear(secuencia):
    secuencia._set_seq(secuencia.seq.ungap('-'))

def abortarYMensaje(mensaje):
    print(mensaje)
    print('Fin del programa')
    sys.exit(0)

def divisiblePor3(secuencia):
    print(f'Longitud de la secuencia: {len(secuencia)}')
    if len(secuencia) % 3 != 0:
        abortarYMensaje("La secuencia no es divisible por 3")


def countCodon(secuencia, codonesBuscados):
    codonesEncontrados = 0
    for i in range(0, len(secuencia), 3):
        if secuencia[i:i + 3] in codonesBuscados:
            codonesEncontrados += 1
            print(codonesBuscados)
            print(f'ubicacion: {i}')
            print(secuencia[i:i + 3])
    return codonesEncontrados

def hayCodonesStop(secuencia):
    print(secuencia)

    codonesEncontrados = countCodon(secuencia, CODONS_STOP)

    if codonesEncontrados > 1:
        abortarYMensaje(f'Hay multiples({codonesEncontrados}) codones de stop en la secuencia')
    elif codonesEncontrados == 1 and secuencia.endswith(CODONS_STOP):
        print('Hay un codon de stop pero esta al final asi que esta bien')
    else:
        abortarYMensaje('No hay codones de stop')

def hayCodonStart(secuencia):

    codonesEncontrados = countCodon(secuencia, CODON_START)

    if codonesEncontrados == 0:
        print('No hay condon de inicio')
    else:
        print('Hay codon de Inicio')


def identidadPorcenaje(align):
    hsp = align.hsps[0]
    return (hsp.identities / hsp.align_length)*100

def eValue(align):
    return align.hsps[0].expect

def score(align):
    return align.hsps[0].score

def bestAlign(aligns):
    #Se ordena segun los siguientes criterios:
    #porcentaje de identidad (decreciente), e value (creciente), score (decreciente)
    ordenado = sorted(aligns, key=lambda a: (identidadPorcenaje(a), -eValue(a), score(a)), reverse=True)
    best = ordenado[0]

    if identidadPorcenaje(best) < 40:
        abortarYMensaje("No se ha encontrado secuencia homologa")

    return best


def pdbId(align):
    title = align.title
    endCode = title.find(" ")
    code = title[0:endCode]
    pdbId = code.split('|')[3]
    print("codigo")
    print(code)
    print(pdbId)
    return pdbId


def soloLetrasARN(seq):
    letrasUnicas = set(seq)
    letrasUnicas.discard('\n')
    print(letrasUnicas)
    print(letrasUnicas == ARNLetras)
    #return letrasUnicas == ARNLetras
    if letrasUnicas != ARNLetras:
        abortarYMensaje("La secuencia debe ser un ARN valido")

def guardarFastaMutacion(seqOriginal,seqMutada):
    fastaMutacionFile = open('secuencia_original_y_mutada.fasta','w')
    fastaMutacionFile.write('>Mutacion\n' + str(seqMutada.translate()) + '\n')
    fastaMutacionFile.write('>Original\n' + str(seqOriginal.translate()) + '\n')
    fastaMutacionFile.close()

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

    desalinear(seq_selected)
    if isProtein(seq_selected.seq):
        print("La secuencia es una proteina")
        sys.exit(0)
    elif isADN(seq_selected.seq):
        print("La secuencia era ADN")
        seq_selected = seq_selected.seq.transcribe()
    elif isARN(seq_selected.seq):
        print("La secuencia era ARN")
        seq_selected = seq_selected.seq


    print(seq_selected)
    quiereMarcoLectura = input("¿Desea definir un marco de lectura? Y/N:  ")
    if quiereMarcoLectura == 'Y':
        marcoLectura = (int(input("Inicio: ")), int(input("Final: ")))
        seq_selected = seq_selected[marcoLectura[0]:marcoLectura[1]]
        print("Seccion seleccionada:")
        print(seq_selected)

    #Chequeos de la secuencia:
    divisiblePor3(seq_selected)
    hayCodonesStop(seq_selected)
    hayCodonStart(seq_selected)
    print(seq_selected)



    #BLAST
    print("Consultando contra BLAST. Aguarde un momento...")
    result_handle = NCBIWWW.qblast("blastx", "pdb", seq_selected)
    with open("my_blast.xml", "w") as save_file:
        blast_results = result_handle.read()
        save_file.write(blast_results)

    result_handle.close()
    
    
    #Parseo BLAST
    result_handle = open("my_blast.xml")
    blast_record = NCBIXML.read(result_handle)

    best = bestAlign(blast_record.alignments)

    print("sequence:", best.title)
    print("IDENTIDAD: ", identidadPorcenaje(best))
    print("e value:", best.hsps[0].expect)
    print("score: ", best.hsps[0].score)


    bestPDBId = pdbId(best)

    #Descarga PDB
    
    pdbl = PDBList()
    pdbWyldTypeDir = pdbl.retrieve_pdb_file(bestPDBId, file_format='pdb' ,pdir='PDB')
    print(pdbWyldTypeName)
    pdbWyldTypeName = pdbWyldTypeDir.split('/')[1].split('.')[0]
    print(pdbWyldTypeName)



    #Mutacion Manual de la secuencia
    print("")
    print("A continuacion se le solcita que haga alguna mutacion a la secuencia ARN")
    print("Para esto se abrira el editor Vi. Al finalizar sus cambios guarde salga del editor.")
    print("Recordatorio: para guardar y salir de Vi, precione Esc e introduzca el comando ':wq'")
    input("presione Enter para continuar")
    (fd, path) = tempfile.mkstemp()
    fp = os.fdopen(fd, 'w')
    fp.write(str(seq_selected))
    fp.close()

    editor = os.getenv('EDITOR', 'vi')
    print(editor, path)
    subprocess.call('%s %s' % (editor, path), shell=True)

    with open(path, 'r') as f:
        seq_mutated = f.read().rstrip()
        print(seq_mutated)


    os.unlink(path)

    # Chequeos de la secuencia mutada
    # todo: bucle de mutacion
    soloLetrasARN(seq_mutated)
    seq_mutated = Seq(seq_mutated, generic_rna)
    divisiblePor3(seq_mutated)
    hayCodonesStop(seq_mutated)
    hayCodonStart(seq_mutated)

    print("Secuencia Original")
    print(seq_selected)
    print("")
    print("Secuencia Mutada")
    print(seq_mutated)

    guardarFastaMutacion(seq_selected, seq_mutated)

    #modelado de proteina mutada
    pdbMutacionName = modelado(pdbWyldTypeName)


    pymol.finish_launching()
    __main__.pymol_argv = ['pymol', '-qc']  # Pymol: quiet and no GUI
    pymol.cmd.load(pdbWyldTypeDir)
    pymol.cmd.load(pdbMutacionName)