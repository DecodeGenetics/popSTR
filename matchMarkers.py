#!/nfs_mount/bioinfo/users/florian/opt/Python-2.7.6/bin/python2.7

from sys import argv

if len(argv) < 5:
    print 'Usage: ' + argv[0] + '<markerNames> <markerCoordsAndNamesFile> <markerInfoFile> <outputFile>'
    exit(1)

names_file = argv[1]    
coords_and_names_file = argv[2]
marker_info_file = argv[3]
output_file = argv[4]
matchedMarkers = 0

markerNames = set()
with open(names_file, 'r') as names_file:
    for line in names_file:
        line = line.split()
        markerNames.add(line[0])

markerNamesAndCoords = []    
with open(coords_and_names_file, 'r') as coords_and_names_file:
    for line in coords_and_names_file:
        if line[0:3] == 'chr':
            line = line.split()
            chrom,begin,end = line[0:3]
            aliases = line[3:len(line)]
            if len(markerNames.intersection(set(aliases)))==1:
                markerNamesAndCoords.append((list(markerNames.intersection(set(aliases)))[0],chrom,begin,end))
            if len(markerNames.intersection(set(aliases)))>1: 
                print 'This is weird, two names for the same marker are in markersToCompare: ' + markerNames.intersection(set(aliases))
                print line
                exit(1)
        else:
            print 'Something is wrong!'
            print line
            exit(1)
print 'Number of markers in memory: ' + str(len(markerNamesAndCoords))

index = 0
with open(marker_info_file, 'r') as marker_info_file:
    with open(output_file, 'w') as outfile:
        line = marker_info_file.readline()
        line = line.split()
        chrom, begin, end = line[0:3]
        while line != '':
            if int(begin) > int(markerNamesAndCoords[index][2]) and int(end) < int(markerNamesAndCoords[index][3]):
                value = ' '.join(map(str, line))
                outfile.write(markerNamesAndCoords[index][0] + ' ' + value + '\n')
                if index > len(markerNamesAndCoords)-1:
                    break
                matchedMarkers +=1
            if int(begin) > int(markerNamesAndCoords[index][2]) and int(end) >= int(markerNamesAndCoords[index][3]):
                index +=1
                if index > len(markerNamesAndCoords)-1:
                    break
                continue           
            line = marker_info_file.readline()
            line = line.split()
            chrom, begin, end = line[0:3]
            
print 'Number of matched markers: ' + str(matchedMarkers)   

