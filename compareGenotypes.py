#!/nfs_mount/bioinfo/users/florian/opt/Python-2.7.6/bin/python2.7

def compareGenotypes(my, lobSTR, benchMark):
    readDepthToSuccess = dict()
    nCompared = 0
    nDiff = 0
    lobSTRnDiff = 0
    diffMarkers = set()
    sameMarkers = set()
    for marker in my:
        myGenotype = my[marker]
        if marker in lobSTR and marker in benchMark:
            nCompared += 1
            lobSTRgenotype = lobSTR[marker]
            benchMarkGenotype = benchMark[marker]
            if lobSTRgenotype[2] != benchMarkGenotype[2]:
                lobSTRnDiff += 1
                if lobSTR[marker][3] in readDepthToSuccess:
                    readDepthToSuccess[lobSTR[marker][3]][0] += 1                    
                else:
                    readDepthToSuccess[lobSTR[marker][3]] = [1,0,0,0]
            else:                
                if lobSTR[marker][3] in readDepthToSuccess:
                    readDepthToSuccess[lobSTR[marker][3]][1] += 1
                else:
                    readDepthToSuccess[lobSTR[marker][3]] = [0,1,0,0]
            if myGenotype[2] == benchMarkGenotype[2]:
                readDepthToSuccess[lobSTR[marker][3]][3] += 1
                sameMarkers.add(marker)
            else:
                readDepthToSuccess[lobSTR[marker][3]][2] += 1
                nDiff += 1
                diffMarkers.add(marker)                
    #print dict1def + ' and ' + dict2def + ' disagree on ' + str(nDiff) + ' out of ' + str(nCompared) + ' genotypes.'
    for depth in readDepthToSuccess:
        lobSTRerrors = readDepthToSuccess[depth][0]
        lobSTRrights = readDepthToSuccess[depth][1]
        myErrors = readDepthToSuccess[depth][2]
        myRights = readDepthToSuccess[depth][3]
        print str(depth) + ' ' + str(lobSTRerrors) + ' ' + str(lobSTRrights) + ' ' + str(myErrors) + ' ' + str(myRights)
    return sameMarkers, diffMarkers             

from sys import argv

if len(argv) != 8:
    print 'Usage: ' + argv[0] + '<benchMarkGenotypes> <myGenotypes> <markerNameToCoordsFile> <successFile> <readDepthFile> <minReadDepth> <lobSTRgenotypes>'
    exit(1)
one_allele_markers = ['D5S2498', 'D5S1468','D7S2248','D9S1793','D11S2368','D13S289','D14S611','D16S3034','D16S516','D22S423','D13S1315','D16S3129','D2S427']
benchmark_genotype_file = argv[1]
my_genotype_file = argv[2]
lobSTR_genotype_file = argv[7]
name_to_coords_file = argv[3]
success_file = argv[4]
read_depth_file = argv[5]
min_readDepth = int(argv[6])
#nComparedMeVsLobSTR = 0
nComparedMe = 0
#nComparedLobSTR = 0
nDiffMe = 0
#nDiffLobSTR = 0
#nDiffMeVsLobSTR = 0
ratioSum = 0
nRatios = 0

markerToReadDepth = dict()
with open (read_depth_file, 'r') as read_depth_file:
    for line in read_depth_file:
        line = line.split()
        markerToReadDepth[(line[0],line[1])] = line[3]

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
        if markerName not in one_allele_markers:
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
        gtInfo = line[9]        
        if int(gtInfo.split(':')[4]) >= min_readDepth:
            if (line[0],line[1]) in markerToReadDepth:
                ratioSum += float(markerToReadDepth[(line[0],line[1])])/float(gtInfo.split(':')[4])
                nRatios += 1
            alleles = gtInfo.split(':')[5]
            A1 = alleles.split('/')[0]
            A2 = alleles.split('/')[1]
            markerName = coordsToName[(line[0],line[1])]
            markerToLobSTRgenotype[markerName] = (A1,A2,abs(int(A1)-int(A2)),gtInfo.split(':')[4])
        line = lobSTR_genotype_file.readline()
if nRatios > 0:
    avRatio = ratioSum/nRatios
    #print 'Average ratio between my coverage and lobSTR coverage: ' + str(avRatio)
else:
    avRatio = 1
    #print 'Average ratio between my coverage and lobSTR coverage: 1'
     
markerToMyGenotype = dict()
default = 0
with open (my_genotype_file, 'r') as my_genotype_file:
    for i in range(0,13):
        my_genotype_file.readline()
    line = my_genotype_file.readline()
    while line != '':
        line = line.split()
        if (line[0],line[1]) not in coordsToName:
            print 'There is a marker in my vcf file which is not in the matchedMarkers file: ' + line[0] + ' ' + line[1]
            exit(1)
        if int(markerToReadDepth.get((line[0],line[1]), default)) >= 1:
            markerName = coordsToName[(line[0],line[1])]
            motifLength = len(line[7].split(';')[1].split('=')[1])
            gtInfo = line[9]        
            alleles = gtInfo.split(':')[0]
            A1 = alleles.split('/')[0]
            A2 = alleles.split('/')[1]
            markerToMyGenotype[markerName] = (A1,A2,round(abs(float(A1)-float(A2))*motifLength), markerToReadDepth.get((line[0],line[1]), default))
        line = my_genotype_file.readline()       

#print 'Starting comparison between me and benchmark data:'
meBenchSame, meBenchDiff = compareGenotypes( markerToMyGenotype, markerToLobSTRgenotype, markerToBenchmarkGenotype)
with open (success_file, 'a') as success_file:
    for marker in meBenchSame:
        success_file.write(marker + ' success \n')
    for marker in meBenchDiff:
        success_file.write(marker + ' fail \n')
#print 'Starting comparison between me and lobSTR:'
#meLobSTRsame, meLobSTRdiff = compareGenotypes( markerToLobSTRgenotype,  markerToMyGenotype, 'my genotypes', 'lobSTR genotypes')
#bothWrong = meBenchDiff.intersection(lobSTRbenchDiff)
#wrongAgree = 0
#for marker in bothWrong: 
    #print marker + ' ' + str(nameToCoords[marker])
    #print 'My genotype: ' +  str(markerToMyGenotype[marker])
    #print 'LobSTR genotype: ' +  str(markerToLobSTRgenotype[marker])
    #if markerToMyGenotype[marker][2] == markerToLobSTRgenotype[marker][2]:
        #wrongAgree += 1
    #print 'Benchmark genotype ' + str(markerToBenchmarkGenotype[marker])
#print 'There are ' + str(len(bothWrong)) + ' markers where both me and lobSTR are wrong and we agree on: ' + str(wrongAgree) + ' of them.'
    
#meWrongLobSTRright = meBenchDiff.difference(lobSTRbenchDiff)
#lobSTRwrongMeRight = lobSTRbenchDiff.difference(meBenchDiff)
#print 'There are ' + str(len(meWrongLobSTRright)) + ' markers where I am wrong and lobSTR is right'
#print 'There are ' + str(len(lobSTRwrongMeRight)) + ' markers where I am right and lobSTR is wrong'
#print 'Markers where I am wrong and lobSTR is right: '
#for marker in meWrongLobSTRright:
    #print marker + ' ' + str(nameToCoords[marker])
    #print 'My genotype: ' +  str(markerToMyGenotype[marker])
    #if marker in markerToLobSTRgenotype:
        #print 'LobSTR genotype: ' +  str(markerToLobSTRgenotype[marker])
    #else:
        #print 'No lobSTR genotype available'
#print 'Markers where lobSTR is wrong and I am right: ' + str(lobSTRwrongMeRight)
