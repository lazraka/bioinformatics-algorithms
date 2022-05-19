import sys
import argparse
import numpy as np
import time
import zipfile


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads

    HINT: This might not work well if the number of reads is too large to handle in memory
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


"""
    TODO: Use this space to implement any additional functions you might need

"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    """
        TODO: Call functions to do the actual read alignment here

    """
import math
import os
#os.chdir('/Users/aminalazrak/Documents/UCLA/Spring 2020/Algorithms in Bioinformatics/CM122_starter_code-master/HP1/practice_W_1/')
os.chdir('/Users/aminalazrak/Documents/UCLA/Spring 2020/Algorithms in Bioinformatics/CM122_starter_code-master/HP2/hw2grad_M_1/')


with open('ref_hw2grad_M_2_chr_1.txt') as file:
    referencegenome=[line.rstrip() for line in file]
referencegenome=referencegenome[1:]
refgenome=''
for segs in referencegenome:
    refgenome=refgenome+segs

with open('reads_hw2grad_M_1_chr_1.txt') as file:
    readpairs=[line.rstrip() for line in file]

readpairs=readpairs[1:]

#choose kmer length
k=50
ksplit=int(k/5)

#split paired reads to individual reads
singleread=[]
for line in readpairs:
    comb=str.split(line,',')
    singleread.append(comb[0])
    singleread.append(comb[1])

#slide window along the genome, gathering up the reads/5 and storing them in a table
hashtable=[]
refgenome_visited=[]
a=0
while a+ksplit<=len(refgenome):
    if refgenome[a:a+ksplit] in refgenome_visited:
        hashtable[refgenome_visited.index(refgenome[a:a+ksplit])].append(a)
        #print(hashtable[refgenome_visited.index(refgenome[a:a+ksplit])])
    else:
        hashtable.append([refgenome[a:a+ksplit], a])
        refgenome_visited.append(refgenome[a:a+ksplit])
    a+=1

#Go through reads, dividing them into kmers/5, looking them up in table
#and store possible read mapping location
kmers=[]
kmerfrom=0 
for kmer in singleread:
    kmers.append([kmer[0:ksplit], 1, kmerfrom])
    kmers.append([kmer[ksplit:ksplit*2], 2, kmerfrom])
    kmers.append([kmer[ksplit*2:ksplit*3],3, kmerfrom])
    kmers.append([kmer[ksplit*3:ksplit*4],4, kmerfrom])
    kmers.append([kmer[ksplit*4:ksplit*5],5,kmerfrom])
    kmerfrom += 1

readmap=[]
kmerpos=[]

for c in kmers:
    if c[0] in refgenome_visited:
        readmap.append([hashtable[refgenome_visited.index(c[0])],c[1]])
        kmerpos.append([c[1], c[2]])
        
#at those locations, compare read to reference genome and find any difference,
#determine if error or SNP
potentialsnps=[]
readseen=[]
readbounds=[]
readcompares=[]
#use only 2900 bases because there is already 30X coverage so would only need 3000 reads
for rank in range(len(kmerpos[0:2900])):
    #if kmerpos[rank][0]==1
        if singleread[kmerpos[rank][1]] not in readseen:
            if readmap[rank][0][1]-((kmerpos[rank][0]-1)*10)>=0 and (readmap[rank][0][1]+(k-((kmerpos[rank][0]-1)*10)))<=len(refgenome):
                readcompare=singleread[kmerpos[rank][1]]
                readseen.append(readcompare)
                readbound=[j for j in range(readmap[rank][0][1]-((kmerpos[rank][0]-1)*10),(readmap[rank][0][1]+(k-((kmerpos[rank][0]-1)*10))))]
                readbounds.append(readbound)
                genomecompare=refgenome[readmap[rank][0][1]-((kmerpos[rank][0]-1)*10):(readmap[rank][0][1]+(k-((kmerpos[rank][0]-1)*10)))]
                if readcompare != genomecompare:
                    if sum(readcompare!=genomecompare for readcompare,genomecompare in zip(readcompare,genomecompare))<=2:
                        for base in range(len(readcompare)):
                            readbase=readcompare[base]
                            genomebase=genomecompare[base]
                            if readbase!=genomebase:
                                potentialsnps.append([genomebase,readbase,(readmap[rank][0][1])+base])
                                readcompares.append(readcompare)

#flatten readbounds to check if the potential snp is overlapped multiple times
boundflattened=[]
def removeNestings(l): 
    for i in l: 
        if type(i) == list: 
            removeNestings(i) 
        else: 
            boundflattened.append(i)
removeNestings(readbounds)

#find duplicates for each potential snp
def duplicates(lst,item):
    return [g for g, x in enumerate(lst) if x==item]

duplicatesnp=[]
snps=[]
for potentialsnp in range(len(potentialsnps)):
    #check all other mapped reads for overlapping regions,if no overlapping regions, throw it out (can't rely on it)
    #if overlapping regions
    if len(duplicates(boundflattened, potentialsnps[potentialsnp][2]))>1:
        duplicatesnp.append(duplicates(boundflattened, potentialsnps[potentialsnp][2]))
        snps.append(potentialsnps[potentialsnp])

#find all the bases from the overlapping reads at that region
consensus=[]
refbase=[]
for dupbase in duplicatesnp:
    snpcompare=[]
    for dups in range(len(dupbase)):
        snpbounds=readbounds[math.floor(dupbase[dups]/50)]
        snpcompare.append(readseen[math.floor(dupbase[dups]/50)][snps[duplicatesnp.index(dupbase)][2]-snpbounds[0]])
    consensus.append(snpcompare)
    refbase.append(snps[duplicatesnp.index(dupbase)][0])

#determine which base is seen the most often at that position with threshold
consensusbase=[]
for con in consensus:
    findcon=[]
    findcon.append(len(duplicates(con,'A'))/len(con))
    findcon.append(len(duplicates(con,'C'))/len(con))
    findcon.append(len(duplicates(con,'G'))/len(con))
    findcon.append(len(duplicates(con,'T'))/len(con))
    if max(findcon)>0.9:
        if findcon.index(max(findcon))==0:
            consensusbase.append('A')
        elif findcon.index(max(findcon))==1:
            consensusbase.append('C')
        elif findcon.index(max(findcon))==2:
            consensusbase.append('G')
        elif findcon.index(max(findcon))==3:
            consensusbase.append('T')
    else:
        consensusbase.append([])
    #else if overlaps same as potentialsnp = realsnp

#determine if snp matches the consensus or not
realsnps=[]
snpseen=[]
for snp1 in snps:
    if snp1[1]==consensusbase[snps.index(snp1)]:
        if snp1[2] not in snpseen:
            realsnps.append(snp1)
            snpseen.append(snp1[2])

#-----------------------------------------------------------------
    #snps = [['A', 'G', 3425]]
    insertions = [['ACGTA', 12434]]
    deletions = [['CACGG', 12]]

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in realsnps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
