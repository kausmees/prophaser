#!/usr/bin/python

# this script is based on the code in interpolate_maps.py this repository: https://github.com/joepickrell/1000-genomes-genetic-maps

import os 

chr = 20

infilename = "/media/kristiina/My Passport/Data/1KGdata/source/markers.txt"
mapfilename = "/path_to_HapMap2_data/HapmapII_GRCh37_RecombinationHotspots/genetic_map_GRCh37_chr"+ str(chr)+".txt"
outfilename = "5_snps_interpolated_HapMap2_map_20_tmp"

infile = open(infilename)
mapfile = open(mapfilename)


outfile = open(outfilename,"w")


posin = list()
rsin = list()
mappos = list()
mapgpos = list()


# for the given list of marker positions
line = infile.readline()
while line:
	line = line.strip().split()
	pos = int(line[0])
	rs = line[1]
	posin.append(pos)
	rsin.append(rs)
	line = infile.readline()



line = mapfile.readline()
line = mapfile.readline()

while line:
	line = line.strip().split()

	if "OMNI" in mapfilename:
		pos = int(line[0]) # 1000 Genomes OMNI map
	elif "HapmapII" in mapfilename:
		pos = int(line[1]) # HapMap2 map
	elif "plink" in mapfilename:
		pos = int(line[3]) # Beagle-supplied PLINK map

	if "OMNI" in mapfilename:
		gpos = float(line[2]) # 1000 Genomes OMNI map
	elif "HapmapII" in mapfilename:
		gpos = float(line[3]) # HapMap2 map
	elif "plink" in mapfilename:
		gpos = float(line[2]) # Beagle-supplied PLINK map




	mappos.append(pos)
	mapgpos.append(gpos)
	line = mapfile.readline()



print("Length posin : " + str(len(posin)))

index1 = 0
index2 = 0
while index1 < len(posin):
	pos = posin[index1]
	rs = rsin[index1]
	if pos == mappos[index2]:
		#the 1000 Genomes site was genotyped as part of the map
		print >> outfile, rs, pos, mapgpos[index2]
		# print rs, pos, mapgpos[index2]
		index1 = index1+1
	elif pos < mappos[index2]:
		#current position in interpolation before marker
		if index2 ==0:
			#before the first site in the map (genetic position = 0)
			print >> outfile, rs, pos, mapgpos[index2]
			index1 = index1+1
		else:
			#interpolate
			prevg = mapgpos[index2-1]
			prevpos = mappos[index2]
			frac = (float(pos)-float(mappos[index2-1]))/ (float(mappos[index2]) - float(mappos[index2-1]))
			tmpg = prevg + frac* (mapgpos[index2]-prevg)
			print >> outfile, rs, pos, tmpg
			# print rs, pos, tmpg
			index1 = index1+1
	elif pos > mappos[index2]:
		#current position in interpolation after marker
		if index2 == len(mappos)-1:
			#after the last site in the map (genetic position = maximum in map, note could try to extrapolate based on rate instead)
			print >> outfile, rs, pos, mapgpos[index2]
			# print rs, pos, mapgpos[index2]
			index1 = index1+1
		else:
			#increment the marker
			index2 = index2+1



### fix 0s at start of map

zero_poss = []

with open(outfilename) as infile:
	line = infile.readline()
	line = line.strip().split()
	dist = float(line[2])
	pos = int(line[1])

	while dist == 0.0:
		zero_poss.append(pos)
		line = infile.readline()
		line = line.strip().split()
		dist = float(line[2])
		pos = int(line[1])
	first_nonzero_pos = pos
	first_nonzero_dist = dist


cM_per_b = first_nonzero_dist/first_nonzero_pos


outfile = open(outfilename[:-4], 'w')

with open(outfilename) as infile:
	line = infile.readline()
	outfile.write(line)
	line = line.strip().split()
	dist = float(line[2])
	pos = int(line[1])
	assert dist == 0.0
	for line in infile:
		line_split = line.strip().split()
		dist = float(line_split[2])
		pos = int(line_split[1])
		if(dist == 0.0):
			line_split[2] = str(pos*cM_per_b)
			outfile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + "\n")
		else:
			outfile.write(line)


os.remove(outfilename)
