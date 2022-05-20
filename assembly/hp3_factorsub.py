from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
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


"""
    TODO: Use this space to implement any additional functions you might need

"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    """
            TODO: Call functions to do the actual assembly here

    """
import os
import math
os.chdir('/Users/aminalazrak/Documents/UCLA/Spring 2020/Algorithms in Bioinformatics/CM122_starter_code-master/HP3/')
#os.chdir('/Users/aminalazrak/Downloads')

with open('reads_hw3all_A_3_chr_1.txt') as file:
    readpairs=[line.rstrip() for line in file]

readpairs=readpairs[1:]

#choose kmer length
k=30

#split paired reads to individual reads
singleread=[]
for line in readpairs:
    comb=str.split(line,',')
    singleread.append(comb[0])
    #singleread.append(comb[1])

#split individual reads into kmers, prekmers is the all the kmers before error correction
prekmers=[]
for read in singleread:
    nuc=0
    while nuc+k<=len(read):
        prekmers.append(read[nuc:nuc+k])
        nuc+=1

#function to find the number of times a prekmer appears in the reads
def duplicates(lst,item):
    return [g for g, x in enumerate(lst) if x==item]

#threshold for deciding if a prekmer is included
repeatmin=3

#determining number of times repeated and creating new array kmer of these correct kmers
kmersrepeat=[]
kmers=[]
count=0
for a in prekmers[0:6500]:
    count+=1
    if count % 500 == 0:
        print(count, " reads done")
    if len(duplicates(prekmers,a))>repeatmin:
        if a not in kmersrepeat:
            kmersrepeat.append([a, len(duplicates(prekmers,a))])
            kmers.append(a)

kmersyn=[]
kmersynrepeat=[]
factor=6
for b in kmersrepeat:
    kmersynrepeat.append([b[0],math.ceil(b[1]/factor)])
    for c in range(math.ceil(b[1]/factor)):
        kmersyn.append(b[0])
        
prefix=[]
store_pre=[]
edges=[]

for a in kmersyn:
    prefix.append(a[:-1])
    if a[:-1] not in store_pre:
        store_pre.append(a[:-1])

def duplicates(lst,item):
    return [g for g, x in enumerate(lst) if x==item]

for b in store_pre:
    if len(duplicates(prefix,b))==1:
        for c in kmersyn:
            if b==c[:-1]:
                #print(b+' -> '+c[1:])
                edges.append(b+' -> '+c[1:])
    else:
        dups=''
        dups_count=0
        for e in kmersyn:
            if b==e[:-1]:
                if dups_count==0:
                    dups=e[1:]
                    dups_count=dups_count+1
                else:
                    dups=dups+','+e[1:]
        #print(b+' -> '+dups)
        edges.append(b+' -> '+dups)
            
start_node=[]
conn_node=[]
index=[None]*len(edges)
num_edges=[]
for node in edges:
    start_node.append(str.split(node,' -> ')[0])
    conn_node.append(str.split(node,' -> ')[1])

for a in range(len(conn_node)):
    if ',' not in conn_node[a]:
        num_edges.append(1)
    else:
        conn_node[a]=str.split(conn_node[a],',')
        conn_node[a]=[h for h in conn_node[a]]
        num_edges.append(len(conn_node[a]))

flatten=[]
def removeNestings(l): 
    for i in l: 
        if type(i) == list: 
            removeNestings(i) 
        else: 
            flatten.append(i)
removeNestings(conn_node)

indegree=[]
outdegree=num_edges
for e in start_node:
    indegree.append(flatten.count(e))

for g in flatten:
    if g not in start_node:
        start_node.append(g)
        num_edges.append(0)
        indegree.append(1)
        conn_node.append([])
        #print(g)

cycle=[]
contig=[]
visited=[0]*len(start_node)

while sum(num_edges)>0:
    var=start_node[next((i for i, x in enumerate(num_edges) if x), None)]
    cycle=[var]

    while num_edges[start_node.index(var)]>0:           
        if num_edges[start_node.index(var)]==1:
            if type(conn_node[start_node.index(var)])==list:
                cycle.append(conn_node[start_node.index(var)][visited[start_node.index(var)]])
                conn_node[start_node.index(var)][visited[start_node.index(var)]]=[]
            else:
                cycle.append(conn_node[start_node.index(var)])
                conn_node[start_node.index(var)]=[]
            visited[start_node.index(var)]+=1
            num_edges[start_node.index(var)]-=1
        else:
            cycle.append(conn_node[start_node.index(var)][visited[start_node.index(var)]])
            num_edges[start_node.index(var)]-=1
            conn_node[start_node.index(var)][visited[start_node.index(var)]]=[]
            visited[start_node.index(var)]+=1
        var=cycle[-1]
    contig.append(cycle)

contigs=[]
contigslong=[]
for m in contig:
    contigstring=m[0]
    for n in range(1,len(m)):
        contigstring=contigstring+m[n][-1]
        contigs.append(contigstring)
    if len(contigstring)>70:
        if contigstring not in contigslong:
            contigslong.append(contigstring)

#------------------------------------

output_fn = args.output_file
zip_fn = output_fn + '.zip'
with open(output_fn, 'w') as output_file:
    output_file.write('>' + args.output_header + '\n')
    output_file.write('>ASSEMBLY\n')
    output_file.write('\n'.join(contigslong))
with zipfile.ZipFile(zip_fn, 'w') as myzip:
    myzip.write(output_fn)
