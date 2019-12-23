# MutaViz

MutaViz es un software de experimentación y visualización de mutaciones sobre secuencias de ADN y ARN.

### Este software le va a permitir al usuario

- Cargar un archivo fasta con encabezado y secuencia ADN / ARN
- Hacer una mutación sobre la secuencia cargada
- Visualizarla

### Requisitos
- Python (_Versión utilizada 3.8_)
- Vi
- [BioPython](https://biopython.org/wiki/Download)
- [Modeller](https://salilab.org/modeller/download_installation.html)
- [PyMol](https://pymol.org/2/support.html?#installation)

### Parametros de funciones
- Busqueda BLAST con NCBIWWW.qblast con [parametros defaults](https://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html#qblast) con excepción de:
    - program = 'blastx'
    - database = 'pdb'

> De entre los resultados se elige como mejor aquel con mayor porcentaje de identidad, si hubiese varios se selecciona entre estos aquel con menor e value, y si aun hubiese varios se selecciona entre estos aquel con mayor score

- Descarga PDB con pdbl.retrieve_pdb_file con [parametros defaults](https://biopython.org/DIST/docs/api/Bio.PDB.PDBList%27.PDBList-class.html#retrieve_pdb_file) con excepcion de:
    - file_format = 'pdb'
    - pdir = 'PDB'

- Modelado con modeller.automodel con [parametros defaults](https://salilab.org/modeller/9v8/manual/node42.html) con excepción de:
    - assess_methods = (assess.DOPE, assess.GA341)
    
>De entre los modelos se elige como mejor aquel con menor DOPE score

### Secuencia uso
1. Ejecutar main.py
2. Seleccionar el archivo FASTA
    > 1. Si en el archivo FASTA hubiese varias secuencias se mostrara una breve descripción de estas, cada una con un ID.
    > 2. Indicar el ID de la secuencia que se quiere utilizar
3. Se consulta si desea definir un marco de lectura (si solo quiere utilizar una porcion de la secuencia).
    1. Si se quiere un marco de lectura, se debe indicar la posicion de inicio (es inclusive y comienza desde 0) y de fin (no inclusive)
    2. Si no se quiere un marco de lectura se utiliza la secuencia entera
> Se realizan chequeos sobre la secuencia (si es divisible por 3, si tiene codon de inicio y un solo codon de parada). Si no cumple con alguna condicion se volvera a consultar por un marco de lectura.
4. Se realiza la consulta a BLAST. Esto puede llevar varios minutos. Y se descarga el PDB.
   > Si no se consigue una secuencia con más de 40% de identidad se aborta el programa
5. Se abre un editor Vi con la secuencia original como ARN, para que el usuario pueda hacer las mutaciones que desee editando la secuencia.
> Recordatorio de Vi:
> - para entrar en modo edicion toque 'i' (para salir de este toque Esc)
> - para moverse: h(izquieda), j(abajo), k(arriba), l(derecha)
> - para guardar toque Esc y escriba el comando ':w'
> - para cerrar toque Esc y escriba el comando ':q'
> - para guardar y cerrar toque Esc y escriba el comando ':wq'

> Se verifica que solo haya caracteres validos y se realizan los chequeos sobre la secuencia (si es divisible por 3, si tiene codon de inicio y un solo codon de parada). Si alguna verificacion falla se le consulta al usuario si desea seguir editando la mutacion que habia hecho o si prefiere comenzar con una mutacion nueva.
6. Comienza el proceso de modelado. En un momento se solicita cuantos modelos queremos que se realicen.
7. Se imprime la ubicacion del PDB de la secuencia original y el de la mutada (la ubicacion es desde la posicion del archivo main.py)

### Autores
- Clever Chuquimia
- Julián Uriel Espinoza
