from modeller import *
from modeller.automodel import *

def modelado(pdb):

    env = environ()

    #FASTA to PIR

    fastaFileDirection = 'secuencia_original_y_mutada.fasta'
    nameFile = 'secuencia_original_y_mutada'

    a = alignment(env, file=fastaFileDirection, alignment_format='FASTA')
    a.write(file='PIR/%s.pir'%nameFile, alignment_format='PIR')


    #Alineamiento templado

    target=nameFile
    template = pdb
    alignName='%s-%s' % (target,template)

    aln = alignment(env)
    mdl = model(env, file="PDB/%s.ent" % template, model_segment=('FIRST:A','LAST:A'))
    aln.append_model(mdl, align_codes=template, atom_files="PDB/%s.ent" % template)
    aln.append(file='PIR/%s.pir'%target, align_codes='Mutacion')
    aln.align2d()
    aln.write(file='ALI/%s.ali' % (alignName), alignment_format='PIR')


    #Generador modelos

    # Homology modeling with multiple templates
    log.verbose()  # request verbose output
    env = environ()  # create a new MODELLER environment to build this model in

    # directories for input atom files
    env.io.atom_files_directory = './:../atom_files'
    env.io.hetatm = True

    print('alignName:   '+alignName)
    print('template:    '+template)
    print('nameFile:    '+nameFile)

    cantModelos = input('\n¿Cuantos modelos desea que Modeller realice?')

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

    print('Termina Modelado')
    return pdbFinal