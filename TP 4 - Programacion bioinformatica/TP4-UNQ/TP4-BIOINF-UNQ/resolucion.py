# Definimos un diccionario K,V donde 'K' es la abreviación de una letra de cada aminoácido y 'V' es el Codón ARN

dic = {'A':'GCU', 'R':'AGA', 'N':'AAU' , 'D':'GAU', 'C':'UGC', 'F':'UUU', 'G':'GGC', 'E':'GAA', 'Q':'CAG', 'H':'CAU',
       'I':'AUU', 'L':'UUA', 'K':'AAA', 'P':'CCU', 'S':'UCU', 'Y':'UAU', 'T':'ACU', 'W':'UGG', 'V':'GUA', 'M':'AUG'}

# TODO: Falta el codón de Stop "UAA", que no tiene valor de una letra.
#items = dic.items()

# Constantes

CADENA = 'AGEKGKKIFVQKCSQCHTVCSQCHTVEKGGKHKTGPNEKGKKIFVQKCSQCHTVLHGLFGRKTGQA'


SEC_2 = 'GTTATAATATTGCTAAAATTATTCAGAGTAATATTGTGGATTAAAGCCACAATAAGATTTATAATCTTAAATGATGGGACTACCATCCTTACTCTCTCCATTTCAAGGCTGAC' \
        'GATAAGGAGACCTGCTTTGCCGAGGAGGTACTACAGTTCTCTTCACAAACAATTGTCTTACAAAATGAATAAAACAGCACTTTGTTTTTATCTCCTGCTTTTAATATGTCCAGTATTC' \
        'ATTTTTGCATGTTTGGTTAGGCTAGGGCTTAGGGATTTATATATCAAAGGAGGCTTTGTACATGTGGGACAGGGATCTTATTTTAGATTTATATATCAAAGGAGGCTTTGTACATGTGG' \
        'GACAGGGATCTTATTTTACAAACAATTGTCTTACAAAATGAATAAAACAGCACTTTGTTTTTATCTCCTGCTCTATTGTGCCATACTGTTGAATGTTTATAATGCATGTTCTGTTTCC' \
        'AAATTTCATGAAATCAAAACATTAATTTATTTAAACATTTACTTGAAATGTTCACAAACAATTGTCTTACAAAATGAATAAAACAGCACTTTGTTTTTATCTCCTGCTTTTAATATGT' \
        'CCAGTATTCATTTTTGCATGTTTGGTTAGGCTAGGGCTTAGGGATTTATATATCAAAGGAGGCTTTGTACATGTGGGACAGGGATCTTATTTTAGATTTATATATCAAAGGAGGCT'


SEC_3 = 'ACTTTGTTTTTATCTCCTGCTCTATTGTGCCATACTGTTGAATGTTTATACCCTACATGGTGCATGTTCTGTTTCCAAATTTCATGAAATCAAAACATTAATTTATTTAAACATTTACT' \
        'TGAAATGTTCACAAACAATTGTCTTACAAAATGAATAAAACAGCACTTTGTTTTTATCTCCTGCTTTTAATATGTCCAGTATTCATTTTTGCATGTTTGGTTAGGCTAGGGCTTAGGGAT' \
        'TTATATATCAAAGGAGGCTTTGTACATGTGGGACAGGGATCTTATTTTAGATTTATATATCAAAGGAGGCT'

# End Constantes


# Funciones auxiliares

def split(word):
    return [char for char in word]


def incrementarSiCorresponde(char):
    return (1 if ((char == "G") | (char == "C")) else 0)


def calcularPorcentaje(total, porcion):
    return ((porcion * 100) / total)

# End Funciones auxiliares

##########################################################################################################
# Ejercicio 1.a
def print_ARN_From_Peptido(strPeptido):
    arn = ''
    listChar = split(strPeptido)

    for c in listChar:
        codon = dic.get(c)

        if (codon != None):
            arn += codon
            ##print(arn)

    return arn

##########################################################################################################
# Ejercicio 1.b
def print_ARN_From_Peptido_Retorna_Tupla(strPeptido):
    arn = ''
    cantGC = 0
    listChar = split(strPeptido)

    for c in listChar:
        cantGC += incrementarSiCorresponde(c)
        codon = dic.get(c)

        if (codon != None):
            arn += codon

    porcentaje = calcularPorcentaje(len(listChar), cantGC)

    # Retorno ambos valores
    return arn, porcentaje

##########################################################################################################
# Ejercicio 2
def buscarPalabraEnSecuencia(palabra, secuencia):
    return secuencia.find(palabra) > 0

#
# MAIN
# Se prueba cada función definida
#
if __name__ == '__main__':
    print(buscarPalabraEnSecuencia('TATAAA', SEC_2))
