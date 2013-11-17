#!/usr/bin/python

import random
from numpy import *
import ConfigParser

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--theta", dest="theta",
                  help="fundamental biodiversity parameter", metavar="FLOAT")

parser.add_option("-n", "--number", dest="N",
                  help="number of sites", metavar="FLOAT")

parser.add_option("-p", "--snumber", dest="P",
                  help="number of supersites", metavar="FLOAT")

parser.add_option("-m", "--migration", dest="migration", help="local migration rate m1,..,mN", metavar="FLOAT")

parser.add_option("-u", "--smigration", dest="u",
                  help="supersite migration rate n1,..,nP", metavar="FLOAT")

parser.add_option("-o", "--outfile", dest="outfile",
                  help="output file", metavar="FILE")

parser.add_option("-s", "--seed", dest="seed",
                  help="rng seed", metavar="integer")

parser.add_option("-x", "--toutfile", dest="toutfile",
                  help="toutput file", metavar="FILE")


(options, args) = parser.parse_args()

theta = float(options.theta)

nN = int(options.N)

nP = int(options.P);

dM = float(options.migration)

aadM = [[dM for i in range(nP)] for j in range(nN)]

dU = float(options.u)

adU = [dU]*nP

nJ = 1000

aanN = [[0 for i in range(nP)] for j in range(nN)]

aanT = [[0 for i in range(nP)] for j in range(nN)]
#number of tables at each site

random.seed(int(options.seed))

nS = 0

allIP = []
allTP = []

aanM = []
anU = [0]*p
anT = [0]*p

n = 0
print str(nN)

for p in range(nP):
	allI = []
	allT = []
	allM = []
	for i in range(nN):
		iList = []
		tList = []
		allI.append(iList)
		allT.append(tList)
	allIP.append(allI)
	allTP.append(allT)
	aanM.append(allM)

def migrate(nT,nM,anM,tList,iList,anT):
	adMSelect = [0]*nT
	adCMSelect = [0]*nT
				
	for k in range(nT):
		adMSelect[k] = double(anM[k])/double(nM)

		if(k > 0):
			adCMSelect[k] = adCMSelect[k - 1] 	

		adCMSelect[k] += adMSelect[k]

	drand3 = random.random()			
	l = 0
	while(drand3 > adCMSelect[l]):
		l=l+1

	tList.append(l)
	iList.append(1)
	anM[l] = anM[l] + 1
	nM = nM + 1
	anT[i] = anT[i] + 1


while (n < nJ):
	for p in range(nP):
		pIList = allIP[p]
		ptList = allTP[p]
		
		for i in range(nN):
			iList = pIList[i]
			tList = ptList[i]		
		
			dMigrate = aadM[p][i]/(double(aanN[p][i]) + aadM[p][i])

			drand = random.random()
		 
			if(drand < dMigrate):
				dMigrateP = adU[p]/(double(anU[p]) + adU[p])
			
				drand2 = random.random()
				if(drand2 < dMigrateP):
					anM.append(1)
					tList.append(nS)
					iList.append(1)
					nS = nS + 1
				else:
				
		else:
			iSelect = [0]*anT[i]
			iCSelect = [0]*anT[i]
			#print str(anT[i])
			for k in range(anT[i]):
				#print str(k) + " " + str(iList[k]) + " " + str(anN[i])
				iSelect[k] = double(iList[k])/double(anN[i])
				if(k > 0):
					iCSelect[k] = iCSelect[k - 1]
				iCSelect[k] += iSelect[k]
       			                 
			drand4 = random.random()
			l = 0
			while(drand4 > iCSelect[l]):
				l=l+1
			iList[l] = iList[l] + 1		
					
		anN[i] = anN[i] + 1
	n = n + 1
	if(n % 100 == 0):
		print str(n) + " " + str(nS) 
	
aanX = []
aanA = []
for i in range(nN):
	anX = [0]*nS
	anA = [0]*nS
	
	for j in range(anT[i]):
		anX[allT[i][j]] += allI[i][j]
		anA[allT[i][j]] += 1

	aanX.append(anX)
	aanA.append(anA)
	
f = open(options.outfile, 'w')
tList = range(1,nN +1)
sampleString = ',Sample'.join(map(str, tList))
f.write('OTUs,Sample'+sampleString+'\n')
u = 0
for i in range(nS):
	sampleString2=''
        for j in range(nN):
        	sampleString2 += ','+str(aanX[j][i])
	f.write('C'+repr(u)+sampleString2+'\n')
	u=u+1
f.closed

f = open(options.toutfile, 'w')
tList = range(1,nN +1)
sampleString = ',Sample'.join(map(str, tList))
f.write('OTUs,Sample'+sampleString+'\n')
u = 0
for i in range(nS):
	sampleString2=''
        for j in range(nN):
        	sampleString2 += ','+str(aanA[j][i])
	f.write('C'+repr(u)+sampleString2+'\n')
	u=u+1
f.closed
