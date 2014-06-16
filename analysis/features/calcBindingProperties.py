# http://pubs.acs.org/doi/pdf/10.1021/bi00035a029

def getEnergyForSgRNA (seq):
    """Gives energy of binding this sgRNA to the DNA. converts input DNA to RNA and computes energy of binding to DNA

    :param seq: 20-mer DNA sequence
    :return Positive floating point representing energy of binding
    """
    # Mapping is RNA sequence to Gibbs free energy
    energyMap={'AA':-1.0,'AC':-5.9,'AG':-1.8,'AU':-0.9
    ,'CA':-0.9,'CC':-2.1,'CG':-1.7,'CU':-0.9,
           'GA':-1.3,'GC':-2.7,'GG':-2.9,'GU':-1.1,
           'UA':-0.6,'UC':-1.5,'UG':-1.6,'UU':-0.2}

    rna=seq.upper().replace('T','U')
    energySum=0.0
    for j in xrange(0,len(seq)-1):
        #grab this and next nucleotide
        energySum+=energyMap[rna[j:j+2]]
    return -1*energySum

if __name__ == "__main__":
    print str(getEnergyForSgRNA('cgccctgacgccgcatggcc'))