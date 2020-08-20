#!/usr/bin/env python3


__author__ = "Phuc HV"
__version__ = "1.0.0"



import sys, os
from shutil import copyfile
from subprocess import Popen
import datetime
import argparse
import re
from Bio import SeqIO




### STYLE

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def OKGREEN (msg):
	return bcolors.OKGREEN + msg + bcolors.ENDC
def OKBLUE (msg):
	return bcolors.OKBLUE + msg + bcolors.ENDC




### FUNCS


def runCMD(cmd, printcmd=True, runcmd = True, wait=True):
	cmd = " ".join(str(i) for i in cmd)
	procc = None
	if printcmd:
		t = '{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
		print ("[{}] {}".format(OKBLUE(t), OKGREEN(cmd)))
	if runcmd:
		procc = Popen(cmd, shell = True)
		if wait:
			procc.wait()
	return procc


def getGenomeSize(file):
	size = 0
	with open (file, "r") as f:
		for l in f:
			l =  l.strip()
			if l.startswith(">") or l.startswith(" "):
				continue
			size += len(l)
	return size


def extractSeq(seq,  location, method="share", minSeq=100):

	#location = [(startpos, endpos),]
	#method = "select" ; "ignore"
	result = []
	posCoor = []
	d = {0:0, len(seq)-1: len(seq)-1}

	for lo in location:
		startpos, endpos = lo

		if method == "share": 
			
			shareSeq = seq[startpos -1: endpos -1]
			if len(shareSeq) < minSeq:
				continue
			posCoor.append(lo) ## export coordinate  later
			result.append(shareSeq)

		if method == "none-share":
			idxStart, idxEnd = startpos -1, endpos -1
			d[idxStart] = idxEnd

	if method == "none-share":
		startpos = sorted(d.keys())


		for i in range(len(startpos)-1):
			k = startpos[i]
			remainStart = d[k]
			remainEnd = startpos[i+1]
			remainSeq = seq[remainStart:remainEnd +1]
			if len(remainSeq) < minSeq:
				continue
			posCoor.append((remainStart,remainEnd))
			result.append(remainSeq)


	return (result, posCoor)




def main():


	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-l', '--list-files', dest="listFiles", help='')
	parser.add_argument('-k', '--kmer-size', dest="kmerSize", default=31, help='')
	parser.add_argument('--method-B', dest="method", default="none-share", help='')
	parser.add_argument('-g', dest="nucmerG", default=10, help='')
	parser.add_argument('-b', dest="nucmerB", default=10, help='')
	parser.add_argument('-c', dest="nucmerC", default=10, help='')
	parser.add_argument('-i', dest="nucmerIdy", default=90, help='')
	parser.add_argument('--overlap-size', dest="overlap", default=100, help='')
	args = parser.parse_args()



	
	### PARAMS
	listFiles = args.listFiles
	kmerSize =  args.kmerSize
	nucmerG = args.nucmerG
	nucmerB = args.nucmerB
	nucmerC = args.nucmerC
	nucmerIdy = args.nucmerIdy
	overlap = args.overlap
	method = args.method

	### TOOLS

	nucmer = "nucmer"
	show_coords = "show-coords"

	###MAIN
	workDir = ""
	prefix = "tmp"
	DICT = {}
	#### copy file to tmp files  in oder to process
	maxGenomeSize = 0
	with open (listFiles, "r") as f:
		i = 0
		for l in f:
			l = l.strip()
			if len(l) > 1: ##  ignore empty lines
				tmp = "{}{}_{}.fna".format(workDir, prefix, i)
				copyfile(l, tmp)
				genomeSize = getGenomeSize(tmp)
				if genomeSize >= maxGenomeSize:
					maxGenomeSize = genomeSize
					largestGenome = tmp
				DICT[tmp] = genomeSize
				i +=1
	####<<


	files = list(DICT.keys())
	genomeA = "genomeA.fasta"
	genomeB = "genomeB.fasta"
	with open (genomeA, "w") as f1, open (genomeB, "w") as f2:
		
		for file in DICT:
			if file == largestGenome:
				copyfile(file, genomeA)
				continue
			with open(file, "r") as f:
				for l in f:
					l = l.strip()
					print (l, file=f2)


			

	#### align 2 genomes
	runCMD([nucmer,"--maxmatch", "-l", kmerSize, "-g", nucmerG, "-b", nucmerB, "-c", nucmerC, "-p", "nucmer_{}".format(prefix),genomeA, genomeB])	
	runCMD([show_coords, "nucmer_{}.delta".format(prefix), "> nucmer_{}.coords".format(prefix)])
	####<<


	#### merge 2 genomes


	commonRegion = {}


	with open ("nucmer_{}.coords".format(prefix), "r") as f:
		countLine = 0
		for l in f:
			countLine +=1

			if countLine >=6: ## above is header
				lobj = re.split(r" {2,}", l.strip())

				
				#[S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]

				S1, E1, S2, E2, L1, L2, IDY, TAG1, TAG2 = int(lobj[0]), int(lobj[1]), int(lobj[3]), int(lobj[4]), int(lobj[6]), int(lobj[7]), float(lobj[9]), lobj[10].split("\t")[0].replace("|", "").strip(), lobj[10].split("\t")[1]

				
				if IDY < nucmerIdy or  L1 < 3*overlap: ## filter regions consider as commons
					continue


				idx = len(commonRegion)

				commonRegion[idx] = [TAG1, S1 + overlap , E1 - overlap]

				if S2 < E2:

					commonRegion[idx] += [TAG2, S2 +overlap, E2 - overlap]
				else:
					commonRegion[idx] += [TAG2, E2 +overlap, S2 - overlap]



	## print GenomeA


	## print GenomeB 
	chromB = {}
	with open (genomeB, "r") as f:
		fasta_seqs  = SeqIO.parse(f, 'fasta')
		for fasta in fasta_seqs:
			chromB[fasta.id] = str(fasta.seq)


	## print none-share region from genomeB:
	shareSortById = {}
	for k in commonRegion:

		shareRegion = commonRegion[k]

		sharedId, start, end = shareRegion[3:]

		if sharedId not in shareSortById:
			shareSortById[sharedId] = [(start, end)]
		else:
			shareSortById[sharedId].append((start, end))


	noneShareBFile = open("genomeB_{}.fasta".format(method), "w")
	for seqId, pos in shareSortById.items():
		result, posCoor = extractSeq(chromB[seqId], pos, method=method)
		for i in range(len(posCoor)):
			header = "{}__{} @length={}@start={}@end={}@method={}".format(seqId,i,(posCoor[i][1]-posCoor[i][0]),posCoor[i][0],posCoor[i][1], method)
			seq = result[i]

			print (">{}\n{}".format(header, seq), file=noneShareBFile)

	##### write to files


	#### join to genomeA


	#with open (genomeA, "r") as f1, open()


if __name__ == "__main__":
	main()
				

				







