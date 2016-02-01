#!/nfs_mount/bioinfo/users/florian/opt/Python-2.7.6/bin/python2.7

from sys import argv
from collections import defaultdict

if len(argv) < 2:
    print 'Usage: ' + argv[0] + ' <old attribute file> <new attribute file>'
    exit(1)
    
old_attribute_file = argv[1]
new_attribute_file = argv[2]
posToAlleles = defaultdict(list)
nDiffAlleles = 0
newAlleles = 0
nMoreReads = 0
nLessReads = 0
nOldMarkers = 0
nNewMarkers = 0
nCompMarkers = 0
nLostAlleles = 0

with open(old_attribute_file, 'r') as old_atts_file:
    pn = old_atts_file.readline()
    line = old_atts_file.readline()
    while line != '':
        if line[0:3] == 'chr':
            nOldMarkers += 1
            chrom,begin,end,motif,refrep,nreads,refseq = line.split()
            for i in range(0, int(nreads)):
                l = old_atts_file.readline()
                posToAlleles[begin].append((l.split()[0],l.split()[6]))                
            line = old_atts_file.readline()
        else:
            print 'Something is wrong!'
            print line
            exit(1)
            
with open(new_attribute_file, 'r') as new_atts_file:
    pn = new_atts_file.readline()
    line = new_atts_file.readline()
    newAlleleSet = set()
    while line != '':
        if line[0:3] == 'chr':
            nNewMarkers += 1
            chrom,begin,end,motif,refrep,nreads,refseq,a1,a2 = line.split()
            if begin in posToAlleles: 
                nCompMarkers += 1
                oldAlleleSet = set([x[0] for x in posToAlleles[begin]])                  
                for i in range(0, int(nreads)):
                    l = new_atts_file.readline()
                    if l.split()[0] not in [x[0] for x in posToAlleles[begin]]:
                        newAlleles += 1
                    newAlleleSet.add(l.split()[0])                                                
                if int(nreads) > len(posToAlleles[begin]):
                    nMoreReads += 1
                    if len(newAlleleSet.difference(oldAlleleSet))>0:
                        print 'More alleles in new attributes file at: ' + str(begin)
                        print newAlleleSet.difference(oldAlleleSet)
                if int(nreads) < len(posToAlleles[begin]):
                    nLessReads += 1
                    if len(oldAlleleSet.difference(newAlleleSet))>0:
                        print 'More alleles in old attributes file at: ' + str(begin)
                        print oldAlleleSet.difference(newAlleleSet)
                        nLostAlleles += 1                    
                if newAlleles > 0:
                    nDiffAlleles += 1
                    newAlleles = 0
                newAlleleSet.clear()          
            else:
                for i in range(0, int(nreads)):
                    new_atts_file.readline()
            line = new_atts_file.readline()       
        else:
            print 'Something is wrong!'
            print line
            exit(1)

print 'Number of markers in old file:' + str(nOldMarkers) + ''
print 'Number of markers in new file:' + str(nNewMarkers) + ''
print 'Total number of compared markers:' + str(nCompMarkers) + ''
print 'Number of markers with more reads in new file: ' + str(nMoreReads) + ''
print 'Number of markers with more reads in old file: ' + str(nLessReads) + ''
print 'Number of markers where I lose alleles in new version: ' + str(nLostAlleles) + ''
print 'Number of markers with new alleles:' + str(nDiffAlleles) + ''
            
