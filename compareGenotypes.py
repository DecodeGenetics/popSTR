#!/nfs_mount/bioinfo/users/florian/opt/Python-2.7.6/bin/python2.7

def compareGenotypes( dict1, dict2, dict1def, dict2def):
    nCompared = 0
    nDiff = 0
    diffMarkers = set()
    sameMarkers = set()
    for marker in dict1:
        dict1genotype = dict1[marker]
        if marker in dict2:
            nCompared += 1
            dict2genotype = dict2[marker]
            if dict1genotype[2] != dict2genotype[2]:
                nDiff += 1
                diffMarkers.add(marker)
                print marker + ' ' + str(nameToCoords[marker])
                print dict1def + ' : ' + str(dict1genotype)
                print dict2def + ' : ' + str(dict2genotype)
            if dict1genotype[2] == dict2genotype[2]:
                sameMarkers.add(marker)                
    print dict1def + ' and ' + dict2def + ' disagree on ' + str(nDiff) + ' out of ' + str(nCompared) + ' genotypes.'
    return sameMarkers, diffMarkers             

from sys import argv

if len(argv) < 6:
    print 'Usage: ' + argv[0] + '<benchMarkGenotypes> <myGenotypes> <lobSTRgenotypes> <markerNameToCoordsFile> <successFile>'
    exit(1)

benchmark_genotype_file = argv[1]
my_genotype_file = argv[2]
lobSTR_genotype_file = argv[3]
name_to_coords_file = argv[4]
success_file = argv[5]
nComparedMeVsLobSTR = 0
nComparedMe = 0
nComparedLobSTR = 0
nDiffMe = 0
nDiffLobSTR = 0
nDiffMeVsLobSTR = 0

coordsToName = dict()
nameToCoords = dict()
with open (name_to_coords_file, 'r') as name_to_coords_file:
    for line in name_to_coords_file:
        line = line.split()
        coordsToName[(line[1],line[2])] = line[0]
        nameToCoords[line[0]] = (line[1],line[2])

markerToBenchmarkGenotype = dict()
with open(benchmark_genotype_file, 'r') as benchmark_genotype_file:
    for line in benchmark_genotype_file:
        markerName,pnId,A1,A2,delta = line.split()
        markerToBenchmarkGenotype[markerName] = (A1,A2,float(delta))

markerToLobSTRgenotype = dict()        
with open(lobSTR_genotype_file, 'r') as lobSTR_genotype_file:
    for i in range(0,23):
        lobSTR_genotype_file.readline()
    line = lobSTR_genotype_file.readline()
    while line != '':
        line = line.split()
        if (line[0],line[1]) not in coordsToName:
            print 'There is a marker in the lobSTR vcf file which is not in the matchedMarkers file: ' + line[0] + ' ' + line[1]
            exit(1)
        markerName = coordsToName[(line[0],line[1])]
        gtInfo = line[9]
        alleles = gtInfo.split(':')[5]
        A1 = alleles.split('/')[0]
        A2 = alleles.split('/')[1]
        markerToLobSTRgenotype[markerName] = (A1,A2,abs(int(A1)-int(A2)))
        line = lobSTR_genotype_file.readline()       
     
markerToMyGenotype = dict()
with open (my_genotype_file, 'r') as my_genotype_file:
    for i in range(0,13):
        my_genotype_file.readline()
    line = my_genotype_file.readline()
    while line != '':
        line = line.split()
        if (line[0],line[1]) not in coordsToName:
            print 'There is a marker in my vcf file which is not in the matchedMarkers file: ' + line[0] + ' ' + line[1]
            exit(1)
        markerName = coordsToName[(line[0],line[1])]
        motifLength = len(line[7].split(';')[1].split('=')[1])
        gtInfo = line[9]        
        alleles = gtInfo.split(':')[0]
        A1 = alleles.split('/')[0]
        A2 = alleles.split('/')[1]
        markerToMyGenotype[markerName] = (A1,A2,round(abs(float(A1)-float(A2))*motifLength))
        line = my_genotype_file.readline()

print 'Starting comparison between me and benchmark data:'
meBenchSame, meBenchDiff = compareGenotypes( markerToMyGenotype,  markerToBenchmarkGenotype, 'my genotypes', 'benchmark genotypes')
with open (success_file, 'a') as success_file:
    for marker in meBenchSame:
        success_file.write(marker + ' success \n')
    for marker in meBenchDiff:
        success_file.write(marker + ' fail \n')
print 'Starting comparison between lobSTR and benchmark data:'
lobSTRbencSame, lobSTRbenchDiff = compareGenotypes( markerToLobSTRgenotype,  markerToBenchmarkGenotype, 'lobSTR genotypes', 'benchmark genotypes')
print 'Starting comparison between me and lobSTR:'
meLobSTRsame, meLobSTRdiff = compareGenotypes( markerToLobSTRgenotype,  markerToBenchmarkGenotype, 'my genotypes', 'lobSTR genotypes')
bothWrong = meBenchDiff.intersection(lobSTRbenchDiff)
print 'There are ' + str(len(bothWrong)) + ' markers where both me and lobSTR are wrong'
for marker in bothWrong: 
    print marker + ' ' + str(nameToCoords[marker])
    print 'My genotype: ' +  str(markerToMyGenotype[marker])
    print 'LobSTR genotype: ' +  str(markerToLobSTRgenotype[marker])
    print 'Benchmark genotype ' + str(markerToBenchmarkGenotype[marker])
    
meWrongLobSTRright = meBenchDiff.difference(lobSTRbenchDiff)
lobSTRwrongMeRight = lobSTRbenchDiff.difference(meBenchDiff)
print 'Markers where I am wrong and lobSTR is right: '
for marker in meWrongLobSTRright:
    print marker + ' ' + str(nameToCoords[marker])
    print 'My genotype: ' +  str(markerToMyGenotype[marker])
    if marker in markerToLobSTRgenotype:
        print 'LobSTR genotype: ' +  str(markerToLobSTRgenotype[marker])
    else:
        print 'No lobSTR genotype available'
print 'Markers where lobSTR is wrong and I am right: ' + str(lobSTRwrongMeRight)
