import gc
import itertools


def run(transactions):
	itemsets={}
	#reverse itemsets
	ritemsets={}
	#itemset unique id counter
	iindex=0
	relations={}

	for i in range(len(transactions)):
		for j in range (i+1,len(transactions)):
			#get relevant transactions
			set1=transactions[i]
			set2=transactions[j]
			# print(i,j)
			#produce the itemset
			intersection=set1&set2
			length=len(intersection)
			
			if (length<2) or (length>255):
			# if (length<2):
				continue

			#add itemset to itemsets and get a unique id
			fitemset=frozenset(intersection)
			if fitemset not in itemsets:
				itemsets[fitemset]=iindex
				ritemsets[iindex]=fitemset
				iindex+=1 
			#add itemset counter to relations for i and j
			if i not in relations:
				relations[i]=set()
			if j not in relations:
				relations[j]=set()
			rindex=itemsets[fitemset]
			
			relations[i].add(rindex)
			relations[j].add(rindex)
	
	return [itemsets,ritemsets,relations]


