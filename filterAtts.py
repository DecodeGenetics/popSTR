#!/nfs_mount/bioinfo/users/florian/opt/Python-2.7.6/bin/python2.7

from sys import argv

if len(argv) < 3:
    print 'Usage: ' + argv[0] + ' <filtered markers file> <attributes file 1> [ ... <attributes file N> ]'
    exit(1)
    
filtered_markers_file = argv[1]
attributes_files = argv[2:]

filtered_markers = []
with open(filtered_markers_file, 'r') as fmfile:
    for line in fmfile:
        line = line.split()
        filtered_markers[len(filtered_markers):] = [(line[0],line[1],line[2])]
filtered_markers = set(filtered_markers)

for attrfile in attributes_files:
    with open(attrfile, 'r') as infile:
        with open(attrfile + 'filtered', 'w') as outfile:
            pn = infile.readline()
            outfile.write(pn)
            line = infile.readline()
            while line != '':
                if line[0:3] != 'chr':
                    print 'Something is wrong!'
                    print line
                    exit(1)
                chrom,begin,end,motif,refrep,nreads,refseq,a1,a2 = line.split()
                nreads = int(nreads)
                if (chrom,begin,end) in filtered_markers:
                    outfile.write(line)
                    for i in range(0, nreads):
                        l = infile.readline()
                        outfile.write(l)
                else:
                    for i in range(0, nreads):
                        infile.readline()
                line = infile.readline()
