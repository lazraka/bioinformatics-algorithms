import os
import math
import numpy

os.chdir('/Users/aminalazrak/Documents/UCLA/Spring 2020/Algorithms in Bioinformatics/CM122_starter_code-master/HP4/hw4_r_4/')

with open('full_genome.txt') as file:
    refgenome=[line.rstrip() for line in file]
refgenome=refgenome[0]

with open('shuffled_reads.txt') as file:
    reads=[line.rstrip() for line in file]

ksplit= 50

hashtable=[]
refgenome_visited=[]
#a=0

exonbounds=[[90792,91412],[571929,572497],[1200222,1200856],[1353310,1354122],[41537,42880]]
for g in range(len(exonbounds)):
    a=exonbounds[g][0]
    while a+ksplit<=exonbounds[g][1]:
        if a % 10000 == 0:
            print(a, " kmers hashed")
        if refgenome[a:a+ksplit] in refgenome_visited:
            hashtable[refgenome_visited.index(refgenome[a:a+ksplit])].append(a)
        else:
            hashtable.append([refgenome[a:a+ksplit], a])
            refgenome_visited.append(refgenome[a:a+ksplit])
        a+=1

readmap=[]
readseen=[]
count=0
for c in reads:
    count+=1
    if count % 100000 == 0:
        print(count, " reads mapped")
    if c in refgenome_visited:
        readmap.append([hashtable[refgenome_visited.index(c)],reads.index(c)])
        readseen.append(c)

print("Finished mapping reads")

junction=[]
#junctions in gene 1
junction.append([refgenome[90792:90873]+refgenome[91139:91234],0,1,0]) #exons 0 and 1
junction.append([refgenome[91139:91234]+refgenome[91314:91412],1,2,0]) #exons 1 and 2

#junctions in gene 2
junction.append([refgenome[571929:571996]+refgenome[572247:572329],0,1,1]) #exons 0 and 1
junction.append([refgenome[572247:572329]+refgenome[572400:572498],1,2,1]) #exons 1 and 2

#junctions in gene 3
junction.append([refgenome[1200222:1200301]+refgenome[1200421:1200496],0,1,2]) #exons 0 and 1
junction.append([refgenome[1200421:1200496]+refgenome[1200779:1200857],1,2,2]) #exons 1 and 2

#junctions in gene 4
junction.append([refgenome[1353310:1353378]+refgenome[1353616:1353679],0,1,3]) #exons 0 and 1
junction.append([refgenome[1353310:1353378]+refgenome[1354069:1354123],0,3,3]) #exons 0 and 3
junction.append([refgenome[1353616:1353679]+refgenome[1353869:1353932],1,2,3]) #exons 1 and 2
junction.append([refgenome[1353869:1353932]+refgenome[1354069:1354123],2,3,3]) #exons 2 and 3

#junctions in gene 5
junction.append([refgenome[41803:41866]+refgenome[42157:42219],1,2,4]) #exons 1 and 2
junction.append([refgenome[42157:42219]+refgenome[42819:42881],2,4,4]) #exons 2 and 4

junctionsegments=[]
junctiontrans=[]
junctionratios=[]
junctionbounds=[80,94,66,81,78,74,67,67,62,62,62,61]

for k in range(len(junction)):
    m=junctionbounds[k]-50+1
    ratio=49
    while (m+ksplit)<=(junctionbounds[k]-1+50):
        junctionsegments.append([junction[k][0][m:(m+ksplit)],junction[k][1],ratio/50,junction[k][2],(50-ratio)/50,junction[k][3]])
        junctiontrans.append(junction[k][0][m:(m+ksplit)])
        m+=1
        ratio-=1
        

exons=[[90792,90872,91139,91233,91314,91412],[571929,571995,572247,572328,572400,572497],[1200222,1200300,1200421,1200495,1200779,1200856],[1353310,1353377,1353616,1353678,1353869,1353931,1354069,1354122],[41537,41597,41803,41865,42157,42218,42470,42539,42819,42880]]

readtoexon=[[0,0,0],[0,0,0],[0,0,0],[0,0,0,0],[0,0,0,0,0]]

for b in range(len(exons)):
    for d in range(int(len(exons[b])/2)):
        for e in range(len(readmap)):
            if readmap[e][0][1]>=exons[b][d*2] and readmap[e][0][1]<=exons[b][(d*2)+1]:
                #print("reached mapped exon")
                readtoexon[b][d]=readtoexon[b][d]+1
            
counter=0
for n in reads:
    counter+=1
    if counter % 100000 == 0:
        print(counter, "reads read")
    if n in junctiontrans:
        #break
        readtoexon[junctionsegments[junctiontrans.index(n)][5]][junctionsegments[junctiontrans.index(n)][1]]+=junctionsegments[junctiontrans.index(n)][2]
        readtoexon[junctionsegments[junctiontrans.index(n)][5]][junctionsegments[junctiontrans.index(n)][3]]+=junctionsegments[junctiontrans.index(n)][4]

##remainingreads=[]
##for p in reads:
##    if p not in readseen and p not in junctiontrans:
##        remainingreads.append(p)

##slide=49
##endjunction=[]
##endjunction.append([refgenome[90792:90873],1,0])
##endjunction.append([refgenome[91314:91412])
##refgenome[571929:571996]
##refgenome[572400:572498]
##refgenome[1200222:1200301]
##refgenome[1200779:1200857]
##refgenome[1353310:1353378]
##refgenome[1353616:1353679]
##refgenome[1354069:1354123]
##refgenome[1353616:1353679]
##refgenome[1354069:1354123]
##refgenome[41803:41866]
##refgenome[42819:42881]

exonlengths=[[80,94,98],[66,81,97],[78,74,77],[67,62,62,53],[60,62,61,69,61]]
isomatrix=[numpy.array([[80/50,0],[94/50,94/50],[98/50,98/50]]),numpy.array([[66/50],[81/50],[97/50]]),numpy.array([[78/50,78/50],[74/50,74/50],[0,77/50]]),numpy.array([[67/50,67/50,0],[62/50,0,62/50],[0,0,62/50],[0,53/50,53/50]]),numpy.array([[0],[62/50],[61/50],[0],[61/50]])]

betas = []
for f in range(len(readtoexon)):
    transiso=numpy.transpose(isomatrix[f])
    inverseiso=numpy.linalg.inv(transiso.dot(isomatrix[f]))
    xmult=inverseiso.dot(transiso)
    beta=xmult.dot(readtoexon[f])
    betas.append(beta.tolist())

abundance=[[0,0],[0],[0,0],[0,0,0],[0]]
for g in range(len(betas)):
    total=sum(betas[g])
    for h in range(len(betas[g])):
        abundance[g][h]=betas[g][h]/total

isoforms=[[refgenome[90792:90873]+refgenome[91139:91234]+refgenome[91314:91413],refgenome[91139:91234]+refgenome[91314:91413]],[refgenome[571929:571996]+refgenome[572247:572329]+refgenome[572400:572498]],[refgenome[1200222:1200301]+refgenome[1200421:1200496],refgenome[1200222:1200301]+refgenome[1200421:1200496]+refgenome[1200779:1200857]],[refgenome[1353310:1353378]+refgenome[1353616:1353679],refgenome[1353310:1353378]+refgenome[1354069:1354123],refgenome[1353616:1353679]+refgenome[1353869:1353932]+refgenome[1354069:1354123]],[refgenome[41803:41866]+refgenome[42157:42219]+refgenome[42819:42881]]]

coverage=[]
for i in range(len(abundance)):
    for j in range(len(abundance[i])):
        coverage.append((isoforms[i][j], abundance[i][j]))

print(coverage)
