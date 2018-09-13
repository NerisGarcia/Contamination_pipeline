#NERIS GARCIA GONZALEZ
# USAGE :
# print "Argument 1: Name of GC_content fro R Script file.\n " \
#      "Argument 2: Name of the TV gene info file.\n " \
#      "Argument 3: Expected family name.\n "
#      "Argument 5: Fasta file.\n" \
#      "Argument 4: Output file name.\n" \


import sys


def select_bad_contigs(family):
    d_badcontigs = {}
    # Open GC_content file
    GC_file = open(sys.argv[1], "r")
    # out = open(sys.argv[3], "w")
    TV_geneinfo_file = open(sys.argv[2], "w")
    # Dictionary for bad contigs
    lFALSEContigs = []
    GC_file.readline()
    # CReate a dictionray with the FALSE contigs to remove
    for line in GC_file:
        lLine = line.split(",")
        #print lLine
        if "FALSE" in str(lLine[8]):
            lFALSEContigs.append(str(lLine[1])[2:-1])
    GC_file.close()
    # Create a dictiory with the contigs from TV file
    d_geneinfo = {}
    #TV_geneinfo_file.readline()
    # Read file
    for line in TV_geneinfo_file:
        lLine = line.split("\t")
       # If contig name in file is in teh FALSE dictionary
        if lLine[4][1:] in lFALSEContigs:
            if lLine[4][1:] in d_geneinfo:
                d_geneinfo[lLine[4][1:]] = d_geneinfo[lLine[4][1:]] + lLine[14]
            else:
                d_geneinfo[lLine[4][1:]] = lLine[14]
    for k, v in d_geneinfo.iteritems():
        lFALSEContigs.remove(k)
        fam = str(d_geneinfo[k])
        famset= set(fam.split(" ")[1:])
        if family not in famset:
            d_badcontigs[k] = str(fam)
    for key in lFALSEContigs:
        d_badcontigs[key] = str("short contig")
    return d_badcontigs



def length(d_bad_contigs, TV_geneinfo_file, fasta, out):
    open(TV_geneinfo_file, "r")
    out_len = open(out, "w")
    # Gene number
    d_ngenes = {}
    TV_geneinfo_file.readline()
    for line in TV_geneinfo_file:
        lline = line.split("\t")
        if lline[4][1:] in d_bad_contigs:
            if lline[4][1:] in d_ngenes:
                d_ngenes[lline[4][1:]] = d_ngenes[lline[4][1:]] + 1
            else:
                d_ngenes[lline[4][1:]] = 1
    #Length
    fasta = open(fasta, "r")
    cab = ""
    dFasta = {}
    seq = ""
    for line in fasta:
        if line.startswith(">"):
            dFasta[cab[1:-1]]=len(seq)
            cab = line
            seq = ""
        else:
            seq += line
    print dFasta
    dDEFFASTA = {}
    for k, v in dFasta.iteritems():
        if k in d_bad_contigs:
            dDEFFASTA[k] = dFasta[k]
    print dDEFFASTA
    for k, v in dDEFFASTA.iteritems():
        tax = ""
        if k in d_bad_contigs:
            tax = str(d_bad_contigs[k])
        if k in d_ngenes:
            out_len.write( k + "\t" +str(dDEFFASTA[k])+"\t"+ str(d_ngenes[k]) +"\t" + tax +"\n")
        else:
            out_len.write(k + "\t" + str(dDEFFASTA[k]) + "\t0\t" + tax + "\n")

def main(argv):
    family = str(sys.argv[3])
    fasta = str(sys.argv[4])
    out =   str(sys.argv[5])
    d_bad_contigs = select_bad_contigs(family)
    length(d_bad_contigs, TV_geneinfo_file, fasta, out)

if __name__ == '__main__':
    main(sys.argv[1:])


