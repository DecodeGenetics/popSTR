#!/nfs_mount/bioinfo/users/florian/opt/Python-2.7.6/bin/python2.7

from sys import argv
from collections import Counter

if len(argv) < 4:
    print 'Usage: ' + argv[0] + ' <filtered markers file> <outputDir> <attributes file 1> [ ... <attributes file N> ]'
    exit(1)
    
filtered_markers_file = argv[1]
out_directory = argv[2]
attributes_files = argv[3:]

filtered_markers = []
with open(filtered_markers_file, 'r') as fmfile:
    for line in fmfile:
        #print line
        line = line.split()
        filtered_markers[len(filtered_markers):] = [(line[0],line[1],line[2])]
filtered_markers = set(filtered_markers)

for attrfile in attributes_files:
    with open(attrfile, 'r') as infile:
        pn = infile.readline().split()[0]
        with open(out_directory + pn + '', 'w') as outfile:
            outfile.write(pn + '\n')
            line = infile.readline()
            while line != '':
                if line[0:3] != 'chr':
                    print 'Something is wrong!'
                    print line
                    exit(1)
                chrom,begin,end,motif,refrep,nreads,refseq,a1,a2 = line.split()
                nreads = int(nreads)
                if (chrom,begin,end) in filtered_markers:
                    reads = []
                    alleles = []
                    for i in range(0, nreads):
                        l = infile.readline()
                        l = l.split()
                        l[0] = str(round(float(len(l[10]))/float(len(motif)),1))
                        reads.append(l)
                        alleles.append(l[0])
                    count = Counter(alleles)
                    a1 = count.most_common()[0][0]
                    if len(count.most_common())>1:
                        if float(count.most_common()[1][1]/count.most_common()[0][1])>0.2:
                            a2 = count.most_common()[1][0]
                        else:
                            a2=a1
                    else:
                        a2=a1
                    outfile.write(chrom + '\t' + begin +'\t'+ end +'\t'+ motif + '\t'+ refrep + '\t' + str(nreads) + '\t' + refseq + '\t' + a1 + '\t' + a2 + '\n')
                    for read in reads:
                       outfile.write('\t'.join(read)+'\n') 
                else:
                    for i in range(0, nreads):
                        infile.readline()
                line = infile.readline()
