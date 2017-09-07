#!/nfs_mount/bioinfo/users/florian/opt/Python-2.7.6/bin/python2.7

from sys import argv

def checkBf(refBf, repeatSeq):
    returnValue = 4
    flankSeq = refBf[-returnValue:]
    while flankSeq in repeatSeq:
        returnValue += 1
        flankSeq = refBf[-returnValue:]
    return returnValue
    
def checkAf(refAf, repeatSeq):
    returnValue = 4
    flankSeq = refAf[:returnValue]
    while flankSeq in repeatSeq:
        returnValue += 1
        flankSeq = refBf[:returnValue]
    return returnValue

if len(argv) < 3:
    print 'Usage: ' + argv[0] + ' <original markers file> <augmented markers file>'
    exit(1)
    
original_markers_file = argv[1]
augmented_markers_file = argv[2]

with open(original_markers_file, 'r') as om_file:
    with open(augmented_markers_file, 'w') as am_file:
        for line in om_file:
            line = line.split()
            refBf = line[5]
            refAf = line[6]
            repeatSeq = line[7]
            nBasesRequiredBf = checkBf(refBf, repeatSeq)
            #if nBasesRequiredBf > 4:
                #print 'Require ' + str(nBasesRequiredBf) + ' before marker at ' + line[1]
            nBasesRequiredAf = checkAf(refAf, repeatSeq)
            #if nBasesRequiredAf > 4:
                #print 'Require ' + str(nBasesRequiredAf) + ' after marker at ' + line[1]
            line.append(str(nBasesRequiredBf))
            line.append(str(nBasesRequiredAf))
            am_file.write('\t'.join(line)+'\n')
