#NERIS GARCIA GONZALEZ
# USAGE :
# print "Argument 1: Name of GC_content fro R Script file.\n " \
#      "Argument 2: Name of the TV gene info file.\n " \
#      "Argument 3: Output file name.\n" \
#      "Argument 4: Expected family name."

import sys


def select_bad_contigs(family):
    d_badcontigs = {}
    # Open GC_content file
    # GC_file = open(sys.argv[1], "r")
    GC_file = open("bj87_GCcontent.csv", "r")
    # out = open(sys.argv[3], "w")
    # TV_geneinfo_file = open(sys.argv[2], "w")
    TV_geneinfo_file = open(
        "/home/neris/Laboratorio/boris/4_anotaciones_TV/Prokka_MIRA_TV/2_blast_prokka_TV/Prokka_BLAST_w16_TV/2_TV_DEFFFF/bj87_gene_info",
        "r")
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
    print lFALSEContigs
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



def length(d_bad_contigs):
    TV_geneinfo_file = open("/home/neris/Laboratorio/boris/4_anotaciones_TV/Prokka_MIRA_TV/2_blast_prokka_TV/Prokka_BLAST_w16_TV/2_TV_DEFFFF/bj87_gene_info", "r")
    out_len = open("bj87_len", "w")
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
    fasta = open("/home/neris/Laboratorio/boris/3_ensamblado_y_mapeado/MIRA4/MyFirstAssembly_MIRA_bj87_sub1_out.unpadded.fasta", "r")
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
    # Comprobar el TV con el contenido en GC
    # GC_file_name = raw_input("Name of GC_content fro R Script file: ")
    # tv_file_name = raw_input("Name of the TV gene info file: ")
    # family = str(sys.argv[4])
    family = "Flavobacteriaceae"
    d_bad_contigs = select_bad_contigs(family)
    #print d_bad_contigs
    length(d_bad_contigs)

if __name__ == '__main__':
    main(sys.argv[1:])


