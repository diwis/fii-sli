import sys
import subprocess
import os.path
import timeit


ontName=sys.argv[1]
threads=sys.argv[2]
sliFolder=sys.argv[3]
fiiFolder=sys.argv[4]
start=int(sys.argv[5])
stop=int(sys.argv[6])
synonyms=sys.argv[7]
interactions=sys.argv[8]

pthres='0.075'


sizes=[str(x) for x in range(start,stop+1) ]

iterations='1000000'

clean=fiiFolder + '/' + ontName + 'Clean.csv'
associations=fiiFolder + '/' + ontName + 'Associations.csv'
itemsets=fiiFolder + '/' + ontName + 'Itemsets.csv'


for size in sizes:
    if not os.path.exists(sliFolder):
        subprocess.call(['mkdir','-p',sliFolder])
    
    slifile=sliFolder + '/sli-' + size + '.txt'
    arguments=['./fii-create-sli',interactions,size,clean,itemsets,associations,iterations,threads,synonyms,slifile,pthres]

    subprocess.call(arguments)