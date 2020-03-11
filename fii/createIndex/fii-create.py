import createItemsets as ap
import os.path
import sys
import gc
import subprocess

ontology=sys.argv[1]
ontologyName=sys.argv[2]
indexFolder=sys.argv[3]

fileName=indexFolder + '/' + ontologyName

#read data
data={}
genes={}
names={}
f=open(ontology,'r')
for line in f:
	line=line.strip().split('|')
	if line[1] not in data:
		data[line[1]]=set()
	data[line[1]].add(line[0])
	names[line[1]]=line[2]
f.close()


for path in data:
	for gene in data[path]:
		if gene not in genes:
			genes[gene]=0
		genes[gene]+=1

#fix transactions and sort them
transactions=[]
trindices={}
paths=[]
i=0

for path in data:
	paths.append(path)

# transactions=sorted(transactions,key=lambda x: len(x))
### This is where the bug was fixed.
### The paths should first be sorted.
### Then the transaction order reflects the transaction id
paths=sorted(paths,key=lambda x: len(data[x]),reverse=False)

#add correct indices for each transaction

for path in paths:
	transactions.append(data[path])
	# paths.append(path)
	trindices[path]=i
	i+=1

print("Starting itemset discovery...")
res=ap.run(transactions)
print("Done.")

#contains itemset unique ids according to itemset
itemsets=res[0]
#the reverse of itemsets (itemset according to unique id)
ritemsets=res[1]
#contains unique ids according to transaction
relations=res[2]

#clear some memory
res=None
gc.collect()

#find all frequent genes
# freqlist=[]
freq=[]
for gene in genes:
	if genes[gene]>1:
		freq.append(frozenset({gene}))
		# freqlist.append(gene)


inames={}
i=0

#add unique identifiers
for gene in freq:
	name='set' + str(i)
	inames[gene]=name
	i+=1

for item in itemsets:
	name='set' + str(i)
	inames[item]=name
	i+=1


assoc={}
clean={}
existing=set()

print('Covering transactions...')
totali=len(data)
i=1.0
for path in paths:
	# print(path,i/totali)
	universe=set()
	assoc[path]=set()
	itemsets_cat=[]

	#get all frequent itemsets created by hackpriori related to the path
	if trindices[path] in relations:
		#add itemset according to their indices
		for num in relations[trindices[path]]:
			itemsets_cat.append(ritemsets[num])
		#sort itemsets from larget to smaller
		itemsets_cat=sorted(itemsets_cat,key=lambda x: len(x), reverse=True)

		#for every itemset, check that it is a subset of the transaction
		for item in itemsets_cat:
			if item.issubset(data[path]):
				#if it does not intersect with another itemset already selected
				if (len(universe&item)==0):
					universe|=item
					assoc[path].add(inames[item])
					existing.add(item)

	for gene in (data[path]-universe):
		geneset=frozenset((gene,))
		if geneset in inames:
			assoc[path].add(inames[geneset])
			universe.add(gene)
	clean[path]=data[path]-universe
	i+=1

itemsets_tmp=itemsets
itemsets=[]
for item in itemsets_tmp:
	if item in existing:
		itemsets.append(item)
#for item in itemsets:
#	print(item)
	
result=freq+itemsets
# result=sorted(itemsets, key=lambda x: len(x), reverse=False)

if not os.path.exists(indexFolder):
	subprocess.call(['mkdir', '-p', indexFolder])

print('Writing data to disk...')
if True:
	f=open(fileName + 'Itemsets.csv','w')
	for item in result:
		name=inames[item]
		for gene in item:
			f.write(gene + '|' + name + '\n')
	f.close


	g=open(fileName + 'Associations.csv','w')
	for path in paths:
		for name in assoc[path]:
			g.write(path + '|' + name + '\n')
	g.close()


	h=open(fileName + 'Clean.csv','w')
	for path in paths:
		for gene in clean[path]:
			h.write(gene + '|' + path + '|' + names[path] + '\n')
	h.close()
