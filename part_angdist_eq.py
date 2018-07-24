#!/usr/bin/env python

# filter particles remove over represented views
# keep the bast particles (highest MaxValueProbDistribution ) **** is this correct?
# TO DO:
#       thurough check that it is removing the right particles
#       test in dave's ribosomes - find more aniostropic datasets to test on

# FIXES/IMPROVEMENTS

# throws away maore parts than necessary: if bin is 1 part away from target and 10 subbins have the max number of particles it will throw away 1 from all 10
# rather than just throwing away 1. don't know if this s necessaraly a bad thing. 


vers = 0.1

import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

###---------function: read the star file get the header, labels, and data -------------#######
def read_starfile(f):
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    for i in alldata:
        if '#' in i:
            labelsdic[i.split('#')[0]] = int(i.split('#')[1])-1
        if len(i.split()) > 3:
            data.append(i.split())
        if len(i.split()) < 3:
            header.append(i.strip("\n"))
    return(labelsdic,header,data)
#---------------------------------------------------------------------------------------------#

if len(sys.argv) not in [4,5]:
    sys.exit('USAGE: part_angdist_eq.py <starfile> <# initial bins> <# std for filtering> --write_graphs')
write_graphs = False
if '--write_graphs' in sys.argv:
    write_graphs = True

(labels,header,data) = read_starfile(sys.argv[1])
nbins = int(sys.argv[2])
nstd = float(sys.argv[3])

def round_to_accrot(x, rotation_acc):
    '''round the number within the value'''    
    return(int(rotation_acc * round(float(x)/rotation_acc)))

def make_vecdic(bin):
    vectors_dic = {}        # vector:[particle,loglikylhood,(rot,tilt,psi),[entire line]]
    for i in data:
        vector = (round_to_accrot(float(i[labels['_rlnAngleRot ']]),bin),round_to_accrot(float(i[labels['_rlnAngleTilt ']]),bin),round_to_accrot(float(i[labels['_rlnAnglePsi ']]),bin))
        if vector not in vectors_dic:
            vectors_dic[vector] = [[i[labels['_rlnImageName ']],i[labels['_rlnMaxValueProbDistribution ']],(i[labels['_rlnAngleRot ']],i[labels['_rlnAngleTilt ']],i[labels['_rlnAnglePsi ']]),i]]
        else:
            vectors_dic[vector].append([i[labels['_rlnImageName ']],i[labels['_rlnMaxValueProbDistribution ']],(i[labels['_rlnAngleRot ']],i[labels['_rlnAngleTilt ']],i[labels['_rlnAnglePsi ']]),i])
    return(vectors_dic)
    
def initial_find_num_bins(nbins):
    '''find the right number of bins for all parts, with increment of 1'''
    print('\n*** getting sampling for {0} (rot,tilt,psi) bins ***'.format(int(nbins)))
    binsize = 1
    test = nbins+1
    while test > nbins:
        vecdic = make_vecdic(binsize)
        test = len(vecdic)
        bincount = []
        for i in vecdic:
            bincount.append(len(vecdic[i]))
        sys.stdout.write('.')
        sys.stdout.flush()
        #print('total/minbincount/maxbincount',test,min(bincount),max(bincount))
        binsize +=1
    print('using bin size of {0} degrees'.format(binsize))
    return(binsize,bincount,vecdic)

# get the inital rough bins
ibinsize,ibincounts,ivec_dic = initial_find_num_bins(int(sys.argv[2]))

# initial plot
if write_graphs == True:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    number = range(len(ibincounts))
    ax.hist(ibincounts)
    plt.ylabel('Particles', fontsize=10)
    plt.xlabel('Sub-bin #', fontsize=10)
    plt.savefig('plt.png')
    plt.close()

# mean # of parts in bin
meanperbin,std = (np.mean(ibincounts),np.std(ibincounts))   # mean parts/bin and standard dev
target_parts = int(meanperbin+(nstd*std))                   # target number of particles in bin

print('there are {0} bins with a mean of {1} particles/bin standard dev of {2}'.format(len(ibincounts),meanperbin,std))
print('the target number of particles/bin is {0} (mean + {1} std) '.format(target_parts,nstd))

#identify the bad bins
badbins = []
for i in ivec_dic:
    if len(ivec_dic[i]) > target_parts:
        badbins.append(i)
print(' {0} bins have more particles than the target'.format(len(badbins)))
print('make sure these numbers are reasonable...')
wait = raw_input('\npress enter to go!')



def make_subvecdic(vecdata,bin):
    vectors_dic = {}        # vector:[particle,loglikylhood,(rot,tilt,psi)]
    for i in vecdata:
        vector = (round_to_accrot(float(i[2][0]),bin),round_to_accrot(float(i[2][1]),bin),round_to_accrot(float(i[2][2]),bin))
        if vector not in vectors_dic:
            vectors_dic[vector] = [[i[0],i[1],(i[2][0],i[2][1],i[2][0]),i[3]]]
        else:
            vectors_dic[vector].append([i[0],i[1],(i[2][0],i[2][1],i[2][0]),i[3]])
    return(vectors_dic)

def find_num_subbins(subvecdata,nbins):
    '''find the right number of bins for all parts, with increment of 1'''
    if int(nbins) == 0:
        nbins=1
    print('*** getting sampling for {0} (rot,tilt,psi) bins ***'.format(int(nbins)))
    binsize = 1
    test = nbins
    while test >= nbins:
        vecdic = make_subvecdic(subvecdata,binsize)
        test = len(vecdic)
        bincount = []
        for i in vecdic:
            bincount.append(len(vecdic[i]))
        sys.stdout.write('.')
        sys.stdout.flush()
        #print('mincount/maxcount',min(bincount),max(bincount))
        binsize +=1
        if min(bincount) == max(bincount):
            break
    print('using bin size of {0} degrees - {1} bin(s)'.format(binsize,test))
    return(binsize,bincount,vecdic)


def gettotal(listolists):
    '''get the total number of parts in a bin'''
    running = []
    for i in listolists:
        running.append(len(listolists[i]))
    return(sum(running))
    
def get_to_target(datain,diccount,target):       # datain = vectors dic from bidbins
    '''get to the target number of particles by removing from most populated bins until targer number is reached'''
    maxno = max(diccount)
    total = gettotal(datain)
    zerocount = []                   #keep track of  how may groups have only 1 part 
    while total > target:
        print('total-target-maxno',total,target,maxno)
        for i in datain:
            if len(datain[i]) >= maxno:
                if len(datain[i]) > 1:
                    datain[i].sort(key=lambda x: x[1],reverse=True)
                    datain[i].pop()
                else:
                    if i not in zerocount:
                        zerocount.append(i)
            if len(zerocount) == len(datain):
                print('stopping because population of all sub-bins = 1')
                return(datain)
                
        maxno -=1
        total = gettotal(datain)
    return(datain)

# each bad bin divide it into n/10 subbins - reduce the number of particles until its total <= mean +n std dev
# set min number of bins to 2
# put these particles in a list of good parts
goodparts = []
bbcount = 1
for i in badbins:
    print('\nworking on bad bin #{0} - {1} particles'.format(bbcount,len(ivec_dic[i])))
    working_subbins = round(len(ivec_dic[i])/10,0)
    if working_subbins < 1:
        working_subbins = 1
    subsize,subcount,subvecdic = find_num_subbins(ivec_dic[i],working_subbins)
    

    number = range(len(subcount))
    
    print('cutting particles to get to the target')
    subvecdic_fixed = get_to_target(subvecdic,subcount,target_parts)
    print ('final total:',gettotal(subvecdic_fixed))
    newcount = []
    for i in subvecdic_fixed:
        newcount.append(len(subvecdic_fixed[i]))    
    number = range(len(newcount))
    if write_graphs == True:
        ax = fig.add_subplot(111)
        ax.bar(number,subcount,color='yellow')
        fig = plt.figure()
        plt.ylabel('Particles', fontsize=10)
        plt.xlabel('Bin #', fontsize=10)
        ax.bar(number,newcount,color='blue')
        plt.savefig('BB{0:03}-fixed.png'.format(bbcount))
        plt.close()
    
    print('**')
    print ('PARTICLES/SUB-BIN OLD',subcount)
    print ('PARTICLES/SUB-BIN NEW',newcount)
    print('**')
    bbcount +=1
    # then add them to the list of good particles:
    for i in subvecdic_fixed:
        for j in subvecdic_fixed[i]:
            goodparts.append(j)

# now put the original good particles in the particle list

for i in ivec_dic:
    if len(ivec_dic[i]) < target_parts:
        for j in ivec_dic[i]:
            goodparts.append(j)

# write the output
output = open('filtered_starfile.star','w')
output.write(header[0])
for i in header[1:]:
    output.write('\n{0}'.format(i))

print('** writing new star file - this could take a while if there is a large number of particles **')

for i in goodparts:
    output.write('\n{0}'.format('\t'.join(i[3])))

output.close()
    
print ('original file had {0} particles'.format(len(data)))
print ('filtered file has {0} particles'.format(len(goodparts)))
print('{0} particles were discarded'.format(len(data)-len(goodparts)))