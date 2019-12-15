import tkinter as tk
from tkinter import filedialog
from modeller import *
from modeller.automodel import *  # Load the automodel class

# Para loguear lo que se imprime en pantalla
import sys

#Inicio seccion para convertir de FASTA to PIR
sys.stdout = open('FastaToPir-Log.txt', 'w')
print('--Comienzo log Fasta to Pir--')

e = environ()

root = tk.Tk()
root.withdraw()

fastaFile = filedialog.askopenfilename(title = 'Seleccione archivo FASTA')
print('Archivo ingresado: ' + fastaFile)

a = alignment(e, file=fastaFile, alignment_format='FASTA')
nombreFile = fastaFile.split('/')[-1].split('.')[0]

a.write(file='%s.pir'%nombreFile, alignment_format='PIR')       #Esto me retorna un archivo .ali
print('Archivo .ali generado: ' + nombreFile)

print('--Fin log Fasta to Pir--')
sys.stdout.close()

#Inicio seccion de alineamiento
sys.stdout = open('Alineamiento-Log.txt', 'w')
print('--Comienzo log Alineamiento--')
pdbFile = filedialog.askopenfilename(title='Seleccione archivo PDB')
target = nombreFile + '.pir'

env = environ()
aln = alignment(env)
mdl = model(env, file=pdbFile, model_segment=('FIRST:A', 'LAST:A'))
aln.append_model(mdl, align_codes=target, atom_files=pdbFile)
aln.append(file=target, align_codes='NM_134326.2|ALBU_RAT')
aln.align2d()
aln.write(file='alineamiento.ali', alignment_format='PIR')
aln.write(file='alineamiento.pap', alignment_format='PAP', alignment_features ="INDICES HELIX BETA")

print('--Fin log Alineamiento--')
sys.stdout.close()

#Inicio seccion para generar n cantidad de modelos
sys.stdout = open('Modelado-Log.txt', 'w')
print('--Comienzo log Modelado--')

log.verbose()  # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = './:../atom_files'
env.io.hetatm = True

a = automodel(env,
              alnfile=target,  # alignment filename
              knowns=('3V03'),  # codes of the templates
              sequence='NM_134326',
              assess_methods=(assess.DOPE, assess.GA341))  # code of the target
a.starting_model = 1  # index of the first model
a.ending_model = 2  # index of the last model
# (determines how many models to calculate)
a.make()  # do the actual homology modeling

print('--Fin log Modelado--')
sys.stdout.close()