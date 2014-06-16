from subprocess import call
import os

#write this to temp.in
def readFile(inFile, seq):
    inFile = open(inFile, 'r')
    numSet = -1
    matches = [[0 for i in xrange(len(seq))] for j in xrange(len(seq))]
    for line in inFile:
        if line[0] == '%' or line == '\n':
            continue
        if numSet == -1:  #this line contains # of values
            numSet = int(line)
        else:
            (a, b, val) = line.split()
            a = int(a) - 1
            b = int(b) - 1
            if a >= len(seq) or b >= len(seq):
                continue  # ignore the values used to represent unpaired probability
            val = float(val)
            matches[a][b] = matches[b][a] = val
    return matches

HOME_DIR = '/Users/davidc/PycharmProjects/parseRegulonDB/'
def callNuPack (seq):
    """Returns 2D array representing probability of interacting
    Prob is a symmetrix matrix

    :param seq: String in upper or lowercase
    :return: matrix is 0-indexed
    """
    seq=seq.lower()
    cacheFile=HOME_DIR+'nupack3.0.3/cache/'+seq+'.ppairs'

    if os.path.exists(cacheFile):
        matches=readFile(cacheFile,seq)
    else:
        tempFile=open(HOME_DIR+'temp.in','w')
        tempFile.write(seq+'\n')
        tempFile.close()

        call([HOME_DIR+'nupack3.0.3/bin/pairs',HOME_DIR+'temp','-material',HOME_DIR+'nupack3.0.3/parameters/rna1995'])
        inFile=HOME_DIR+'temp.ppairs'
        matches = readFile(inFile, seq)
        os.unlink(HOME_DIR+'temp.in')
        os.rename(HOME_DIR+'temp.ppairs',cacheFile)
    return matches

if __name__ == '__main__':
    test = callNuPack('gcgagaaaccttattaacca')
    print test
