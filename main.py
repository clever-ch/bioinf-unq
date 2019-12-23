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

#Elimina los guiones de la secuencia si es que los hubiese
def desalinear(secuencia):
    secuencia._set_seq(secuencia.seq.ungap('-'))

def abortarYMensaje(mensaje):
    print(mensaje)
    print('Fin del programa')
    sys.exit(0)

#Chequea que la secuencia sea multiplo de 3
def divisiblePor3(secuencia):
    print(f'Longitud de la secuencia: {len(secuencia)}')
    ret = True
    if len(secuencia) % 3 != 0:
        print("La secuencia no es divisible por 3")
        ret = False
    print(f'Divisible por 3: {ret}')
    return ret


#Cuenta la cantidad de apariciones de un conjunto de codones en una secuencia
def countCodon(secuencia, codonesBuscados):
    codonesEncontrados = 0
    for i in range(0, len(secuencia), 3):
        if secuencia[i:i + 3] in codonesBuscados:
            codonesEncontrados += 1
    return codonesEncontrados

#Vefifica si la secuencia tiene uno o varios codones de stop;
# y si solo hubiese uno, si el codon de stop se encuentra al final de la secuencia
def hayCodonesStop(secuencia):
    print(secuencia)

    codonesEncontrados = countCodon(secuencia, CODONS_STOP)

    if codonesEncontrados > 1:
        print(f'Hay multiples({codonesEncontrados}) codones de stop en la secuencia')
        ret = False
    elif codonesEncontrados == 1 and secuencia.endswith(CODONS_STOP):
        print('Hay un codon de stop pero esta al final asi que esta bien')
        ret = True
    else:
        print('No hay codones de stop')
        ret = False
    print(f'Hay codones stop: {ret}')
    return ret

#Verifica que haya un codon de inicio
def hayCodonStart(secuencia):

    codonesEncontrados = countCodon(secuencia, CODON_START)

    if codonesEncontrados == 0:
        print('No hay condon de inicio')
        ret = False
    else:
        print('Hay codon de Inicio')
        ret = True
    print(f'Hay Codones start: {ret}')
    return ret

#Verifica que la secuencia sea valida
def chequeosSecuenciaSonCorrectos(secuencia):
    return divisiblePor3(secuencia) and hayCodonesStop(secuencia) and hayCodonStart(secuencia)


#Retorna el porcentaje de identidad de uno de los resultados de la consulta a BLAST
def identidadPorcenaje(align):
    hsp = align.hsps[0]
    return (hsp.identities / hsp.align_length)*100

#Retorna el e value de uno de los resultados de la consulta a BLAST
def eValue(align):
    return align.hsps[0].expect

#Retorna el score de uno de los resultados de la consulta a BLAST
def score(align):
    return align.hsps[0].score


#Retorna el mejor de los resultados de la consulta a BLAST
#Si ninguno supera el 40% de identidad se aborta el programa
def bestAlign(aligns):
    #Se ordena segun los siguientes criterios:
    #porcentaje de identidad (decreciente), e value (creciente), score (decreciente)
    ordenado = sorted(aligns, key=lambda a: (identidadPorcenaje(a), -eValue(a), score(a)), reverse=True)
    best = ordenado[0]

    if identidadPorcenaje(best) < 40:
        abortarYMensaje("No se ha encontrado secuencia homologa")

    return best

#Retorna el ID de PDB de uno de los resultados de la consulta a BLAST
def pdbId(align):
    title = align.title
    endCode = title.find(" ")
    code = title[0:endCode]
    pdbId = code.split('|')[3]
    print(f'\nPDB ID {pdbId}')
    return pdbId


#Verificacion de que la secuencia sea un ARN valido
def soloLetrasARN(seq):
    letrasUnicas = set(seq)
    letrasUnicas.discard('\n')
    ret = letrasUnicas == ARNLetras
    if not ret:
        print("La secuencia debe ser un ARN valido")
        caracteresInvalidos = letrasUnicas - ARNLetras
        print("Caracteres invalidos encontrados:")
        print(caracteresInvalidos)
    return ret

#Genera un archivo FASTA con las proteinas de la seqcuencia original y la mutada
def guardarFastaMutacion(seqOriginal,seqMutada):
    fastaMutacionFile = open('secuencia_original_y_mutada.fasta','w')
    fastaMutacionFile.write('>Mutacion\n' + str(seqMutada.translate()) + '\n')
    fastaMutacionFile.write('>Original\n' + str(seqOriginal.translate()) + '\n')
    fastaMutacionFile.close()



#*******************************
#MAIN
#*******************************
if __name__ == '__main__':

    root = tk.Tk()
    root.withdraw()

    #Solicito la ubicacion del archivo FASTA
    fastaFileDirection = filedialog.askopenfilename(title = 'Seleccione archivo FASTA')
    print("Archivo ingresado")
    print(fastaFileDirection)
    secuencias = list(SeqIO.parse(fastaFileDirection, "fasta"))

    #Si en el archivo FASTA habia mas de una secuencia consulto cual desea utilizar
    #Si solo habia una secuencia, directamente se utiliza esa
    if len(secuencias) > 1:
        print("Se encontraron varias secuencias en el archivo.")
        index = 0

        #Imprimo una breve descripcion de cada secuencia
        for seq_record in secuencias:
            print("ID: " + str(index))
            index += 1
            print(seq_record.id)
            print(repr(seq_record.seq))
            print(f'longitud: {len(seq_record)}')

        #Solicito que secuencia desea utilizar
        #Repito la solicitud hasta que se ingrese un ID valido
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

    #Elimino los alineamientos si es que los hubiese
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

    #Consulto si desea definir un marco de lectura de la secuencia.
    #Se realizan chequeos sobre la secuencia. Verifica si es divisible por 3, si tiene codon de inicio y un solo codon de parada.
    #Si no cumple con alguno de los cheques se vuelve a consultar por un marco de lectura.
    seq_es_incorrecta = True
    while seq_es_incorrecta:
        quiereMarcoLectura = input("¿Desea definir un marco de lectura? Y/N:  ").upper()
        if quiereMarcoLectura == 'Y':
            marcoLectura = (int(input("Inicio: ")), int(input("Final: ")))
            seleccion = seq_selected[marcoLectura[0]:marcoLectura[1]]
            print("\nSeccion seleccionada:")
            print(seleccion)
            seq_es_incorrecta = not chequeosSecuenciaSonCorrectos(seleccion)
            if not seq_es_incorrecta:
                seq_selected = seleccion
        else:
            print("\nChequeos de secuencia:")
            seq_es_incorrecta = not chequeosSecuenciaSonCorrectos(seq_selected)

    print(seq_selected)


    #*******************************
    #BLAST
    #*******************************
    #Se consulta contra BLAST y se guarda el resultado en un archivo 'my_blast.xml'
    
    print("Consultando contra BLAST. Puede llevar bastante tiempo. Aguarde un momento...")
    result_handle = NCBIWWW.qblast("blastx", "pdb", seq_selected)
    with open("my_blast.xml", "w") as save_file:
        blast_results = result_handle.read()
        save_file.write(blast_results)

    result_handle.close()
    

    # *******************************
    #Parseo BLAST
    # *******************************
    result_handle = open("my_blast.xml")
    blast_record = NCBIXML.read(result_handle)

    #se selecciona el mejor de los resultados devueltos por BLAST
    #siendo este aquel con mayor porcentaje de identidad,
    #                   si hubiese varios se seleccionaentre estos aquel on menor e value,
    #                   y si aun hubiese varios se seleecciona entre estos aquel con mayor score
    best = bestAlign(blast_record.alignments)

    #breve descripcion del mejor resultado de BLAST
    print("sequence:", best.title)
    print("IDENTIDAD: ", identidadPorcenaje(best))
    print("e value:", best.hsps[0].expect)
    print("score: ", best.hsps[0].score)


    bestPDBId = pdbId(best)

    # *******************************
    #Descarga PDB
    # *******************************
    pdbl = PDBList()
    pdbWyldTypeDir = pdbl.retrieve_pdb_file(bestPDBId, file_format='pdb' ,pdir='PDB')
    #Extraigo el nombre de la direccion absoluta del archivo. Luego se le es pasasda como parametro a Modeller
    pdbWyldTypeName = pdbWyldTypeDir.split('/')[1].split('.')[0]


    #*******************************
    #Mutacion Manual de la secuencia
    #*******************************
    #Se abre un editor Vi con la secuencia ARN para que el usuario pueda editarla como desee.
    #Se verifica que la cadena resultante solo utilice caracteres validos, en caso de que no, se notifica al usuario cuales son aquellos que caracteres que no puede utilizar.
    #Se verifica que la secuencia sea valida.
    #Si alguna verificacion falla se le consulta al usuario si desea seguir editando la mutacion que habia hecho o si prefiere comenzar con una mutacion nueva.
    print("\nA continuacion se le solcita que haga alguna mutacion a la secuencia ARN")
    print("Para esto se abrira el editor Vi. Al finalizar sus cambios guarde y salga del editor.")
    print("Recordatorio: para guardar y salir de Vi, precione Esc e introduzca el comando ':wq'")
    input("presione Enter para continuar")
    (fd, path) = tempfile.mkstemp()
    fp = os.fdopen(fd, 'w')
    fp.write(str(seq_selected))
    fp.close()

    editor = os.getenv('EDITOR', 'vi')
    print(editor, path)

    seqMutada_es_incorrecta = True
    while seqMutada_es_incorrecta:
        subprocess.call('%s %s' % (editor, path), shell=True)

        with open(path, 'r') as f:
            seq_mutated = f.read().rstrip()
            print(seq_mutated)



        # Chequeos de la secuencia mutada
        print("\nChequeos de secuencia:")
        letrasARN = soloLetrasARN(seq_mutated)
        seq_mutated = Seq(seq_mutated, generic_rna)

        seqMutada_es_incorrecta = not (letrasARN and chequeosSecuenciaSonCorrectos(seq_mutated))

        if seqMutada_es_incorrecta:
            print("¿Quiere continuar modificando la mutacion? Si no, comenzara una mutacion nueva")
            accion = input('Y/N:  ').upper()
            if accion != 'Y':
                (fd, path) = tempfile.mkstemp()
                fp = os.fdopen(fd, 'w')
                fp.write(str(seq_selected))
                fp.close()

    os.unlink(path)

    print("Secuencia Original")
    print(seq_selected)
    print("")
    print("Secuencia Mutada")
    print(seq_mutated)

    #Crea un archivo FASTA en el cual esta las proteinas de la secuencia original y de la mutada
    guardarFastaMutacion(seq_selected, seq_mutated)

    # *******************************
    #Modelado de proteina mutada
    # *******************************
    #Se le pide al usuario la cantidad de modelos que desea que se realisen
    #Se selecciona como mejor modelo aquel con mayor DOPE score
    pdbMutacionName = modelado(pdbWyldTypeName)

    # *******************************
    #Visualizacion PDBs con Pymol
    # *******************************
    pymol.finish_launching()
    __main__.pymol_argv = ['pymol', '-qc']  # Pymol: quiet and no GUI
    pymol.cmd.load(pdbWyldTypeDir)
    pymol.cmd.load(pdbMutacionName)