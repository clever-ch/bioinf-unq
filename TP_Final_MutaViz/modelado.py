#Dependencias para GUI
import tkinter as tk
from tkinter import filedialog, messagebox

from modeller import * # Load standard Modeller classes
from modeller.automodel import *  # Load the automodel class

def modelado(pdb):

    env = environ()

    #FASTA to PIR

    #target='pdb_70'
    #target='aln_rat_3v03'

    root = tk.Tk()
    root.withdraw()

    #fastaFileDirection = filedialog.askopenfilename(title = 'Seleccione archivo FASTA')
    fastaFileDirection = 'secuencia_original_y_mutada.fasta'
    print(fastaFileDirection)
    #fastaFileDirection = fastaFileDirection+".algo.otra"
    #print(fastaFileDirection)

    #de la direccion absoluta extraigo el nombre del archivo, sin su extension
    #nameFile = '.'.join(fastaFileDirection.split('/')[-1][::-1].split('.')[1:])[::-1]
    nameFile = 'secuencia_original_y_mutada'
    print(nameFile)


    a = alignment(env, file=fastaFileDirection, alignment_format='FASTA')
    #a.write(file='%s.ali'%target, alignment_format='PIR')
    a.write(file='PIR/%s.pir'%nameFile, alignment_format='PIR')


    #Alineamiento templado

    #target='aln_pig_3v03'
    #template='3V03_A'
    #target='aln_rat_3v03'
    #target='aln_rat_3v03_auxiliar'
    #template='3V03'
    target=nameFile
    #template='pdb3lee' #todo Ver de donde sale el nombre------------------------------------------------------------
    template = pdb
    alignName='%s-%s' % (target,template)

    aln = alignment(env)
    mdl = model(env, file="PDB/%s.ent" % template, model_segment=('FIRST:A','LAST:A'))
    aln.append_model(mdl, align_codes=template, atom_files="PDB/%s.ent" % template)
    #aln.append(file='%s.ali'%target, align_codes=target)
    #aln.append(file='%s.pir'%target, align_codes="NM_134326")
    #aln.append(file='%s.pir'%target, align_codes="NM_134326.2")
    #aln.append(file='PIR/%s.pir'%target, align_codes='all')
    aln.append(file='PIR/%s.pir'%target, align_codes='Mutacion')
    aln.align2d()
    aln.write(file='ALI/%s.ali' % (alignName), alignment_format='PIR')
    #aln.write(file='%s-%s.pap' % (target,template), alignment_format='PAP', alignment_features ="INDICES HELIX BETA")


    #Generador modelos

    # Homology modeling with multiple templates

    # para logear lo que imprima en pantalla
    #import sys

    #sys.stdout = open('file', 'w')
    print('Comienza')

    log.verbose()  # request verbose output
    env = environ()  # create a new MODELLER environment to build this model in

    # directories for input atom files
    env.io.atom_files_directory = './:../atom_files'
    env.io.hetatm = True

    print('alignName:   '+alignName)
    print('template:    '+template)
    print('nameFile:    '+nameFile)

    cantModelos = input('¿Cuantos modelos desea que Modeller realice?')

    a = automodel(env,
                  alnfile='ALI/%s.ali' % (alignName),  # alignment filename
                  knowns=(template),  # codes of the templates
                  sequence='Mutacion',
                  assess_methods=(assess.DOPE, assess.GA341))  # code of the target
    a.starting_model = 1  # index of the first model
    a.ending_model = int(cantModelos)  # index of the last model
    # (determines how many models to calculate)
    a.make()  # do the actual homology modeling

    # Get a list of all successfully built models from a.outputs
    ok_models = [x for x in a.outputs if x['failure'] is None]

    # Rank the models by DOPE score
    key = 'DOPE score'

    ok_models.sort(key=lambda a: a[key])

    # Get top model
    m = ok_models[0]
    print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))

    pdbFinal = str(m['name'])

    #print('PDB Final')
    #print(pdbFinal)
    print('Termina Modelado')
    return pdbFinal