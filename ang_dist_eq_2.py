#!/usr/bin/env python

import sys
import numpy as np

def read_3p1star(infile):
	header = []
	inDMCloop = False
	labels,data = {},[]
	dta = open(infile,'r').readlines()
	for i in dta:
		i=i.split()
		if 'data_particles' in i:
			inDMCloop = True
		if inDMCloop == True and len(i) > 0:
			if '_rln' in i[0]:
				labels[i[0]] = int(i[1].replace('#',''))-1
				header.append(' '.join(i))
			elif '@' in ''.join(i):		
				data.append(i)
			else:
				header.append(' '.join(i))
		else:
			header.append(' '.join(i))
	return(labels,header,data)


###---------function: read the star file get the header, labels, and data -------------#######
def read_starfile_new(f):
    inhead = True
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    count = 0
    labcount = 0
    for i in alldata:
        if '_rln' in i:
            labelsdic[i.split()[0]] = labcount
            labcount +=1
        if inhead == True:
            header.append(i.strip("\n"))
            if '_rln' in i and '#' in i and  '_rln' not in alldata[count+1] and '#' not in alldata[count+1]:
                inhead = False
        elif len(i.split())>=1:
            data.append(i.split())
        count +=1
    
    return(labelsdic,header,data)
#---------------------------------------------------------------------------------------------#
errmsg ='USAGE: ang-eq.py <selected from class2d starfile> <number of standard deviations>'

### get data
try:
	vers = open(sys.argv[1],'r').readlines()[:100]
	r3p1 = False
	for i in vers:
		if i.strip('\n') == '# version 30001':
			r3p1= True
			break
	if r3p1 == True:
		labels,header,data = read_3p1star(sys.argv[1])
	else:
		labels,header,data = read_starfile_new(sys.argv[1])
except:
	sys.exit('\nERROR: problem reading starfile\n{}'.format(errmsg))
	
try:
	numsd = float(sys.argv[2])
except:
	sys.exit('\nERROR: Bad number of standard deviations\n{}'.format(errmsg))

classes = {}
for i in data:
	try:
		classes[int(i[labels['_rlnClassNumber']])].append(i)
	except:
		classes[int(i[labels['_rlnClassNumber']])] = [i]
## totals
classtotals = []
for i in classes:
	classtotals.append(len(classes[i]))

## sort by MVPD
for i in classes:
		classes[i].sort(key=lambda x: x[labels['_rlnMaxValueProbDistribution']],reverse=True)

maxppc = int(round(np.mean(classtotals)+(numsd*np.std(classtotals)),0))

## trim the classes

print('''
          particles     mean MaxValProbDist     %
class   start   finish    original   final    culled''')

cl = list(classes)
cl.sort()
icount,fcount = [],[]
for i in cl:
	icount.append(len(classes[i]))
	if len(classes[i]) < maxppc:
		nnum = len(classes[i])
		fcount.append(len(classes[i]))
	else:
		nnum = maxppc
		fcount.append(maxppc)

	print('{0: 3d}     {1: 3d}     {2: 3d}       {3:0.3f}      {4:0.3f}    {5:0.2f}'.format(i,len(classes[i]),nnum,round(np.mean([float(x[labels['_rlnMaxValueProbDistribution']]) for x in classes[i]]),3),round(np.mean([float(x[labels['_rlnMaxValueProbDistribution']]) for x in classes[i][0:nnum+1]]),3),1-(float(nnum)/len(classes[i]))))
	classes[i] = classes[i][0:nnum+1]
## write output:
output = open('ang-eq_classes.star','w')
for i in header:
	output.write('{0}\n'.format(i))
for i in classes:
	for j in classes[i]:
		output.write('{0}\n'.format('  '.join(j)))
### screen barf
print('\nMinimum: {0}\tMaximum: {1}'.format(np.min(classtotals),np.max(classtotals)))
print('Mean   : {0:0.1f}\tStD:     {1:0.1f}'.format(np.mean(classtotals),np.std(classtotals)))
print('Maximum allowed particles/class: {0}'.format(maxppc) )
print('original set had {} particles'.format(sum(icount)))
print('culled set has {} particles'.format(sum(fcount)))
print('overall retention {0:0.3f}'.format(float(sum(fcount))/sum(icount)))
