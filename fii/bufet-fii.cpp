/*
 * Copyright 2020 Konstantinos Zagganas and Thanasis Vergoulis
 * for the Information Management Systems Institute(IMSI) - "Athena" Research Center
 * 
 * This file is part of diwis/fii-sli.
 *
 * Foobar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * diwis/fii-sli is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 */

#define ARGUMENTS 10
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <random>
#include <vector>
#include <thread>
#include <cstdio>
#include <bitset>
#include <algorithm>
#include <string.h>

#define BSIZE 25000

using namespace std;

/*
 *
 *
 * Data types defined below
 *
 *
 */

/*
 * Interaction node used for temporary storage of interations
 */
struct interaction
{
	string gene;
	interaction * next;
};

/*
 * GO category node used for temporary storage of GO-gene associations
 */
struct goCatNode
{
	string category;
	string name;
	long int intersection=0;
	long int go_size=0;
	float overlap_proportion=0;
	double mean_overlap=0;
};

/*
 * Pvalue struct used for BH FDR correction
 */
struct pNode
{
	double pvalue;
	int index;
	string fdr="-";
};

struct counterType
{
	unsigned char * charPointer;
	int * intPointer;
};

/*
 * Hash table containing gene to internal id associations
 */
typedef unordered_map <string,int> genedata;
/*
 * Temporary hash table containing mirna-gene interactions as a linked list for each mirna
 */
typedef unordered_map <int, interaction * > geneInteractions;
/*
 * Temporary hash table containing GO-gene assosiations as a linked list for each category
 */
typedef unordered_map <string, interaction * > goTemp;
/*
 * Auxiliary pointer to a bitset
 */
typedef bitset<BSIZE> * bits;
/*
 * Auxiliary pointer to a vector
 */
typedef vector<int> * vec;
/*
 * Hash table containg interactions for each miRNA
 *
 * @key: internal miRNA id
 * @value: bit vector representing all internal gene ids in set
 */
typedef unordered_map <int, bits> final_interactions_type;
/*
 * Hash table containg GO-gene associations
 *
 * @key: go category
 * @value: bit vector representing all internal gene ids in set
 */
typedef unordered_map <string, vec> goGenes_type;
typedef unordered_map <int, vec> itemset_type;
/*
 * List of tokens
 */
typedef vector<string> token_list;
/*
 * Hash table to contain names for GO categories
 */
typedef unordered_map <string,string> go_names_type;

typedef vector<string> * string_list;

typedef vector<int>::iterator vit;

/*
 * Data structures
 */

genedata genes, miRNAs, frequent;
geneInteractions interactions;
goGenes_type goGenes,goFreq;

final_interactions_type finalInteractions;
bits * map_all;
vector<goCatNode> checkGO, noCheckGO;
vector<float> i_counts;
vector<long int> pvalues;
go_names_type goNames;
unordered_map<string,string_list> synonyms;
unordered_set<string> genesGo;
vector <pNode *> fdrList;
itemset_type freqGenes, freqLargeAssoc;
int uGenesAssoc[BSIZE];
vector<int>uGenes,nuGenes;
counterType * freqCounters;
int totalu=0, totalnu=0;
/*
 * Global integer variables
 */
int gene_count=1, miRNA_count=1,group_size=0,GOcount=0,freq_count=0, thread_count,total_input;
long int iterations;
string species;

/*
 * Function definitions(see explanations above each function)
 */
void getInteractions(string);
void getClean(string);
void getRandom(int,int,int);
void writeOutput(string);
void findIntersections(int,int);
void fixInteractions();
void getMirnas(string);
void getFrequent(string);
void getAssociations(string);
void frequentIntersectionsSingle(int,int);
void frequentIntersectionsMulti(int,int);
void calculateCounts();
string trim(string mystr);
string trim_chars_left(string,string);
void getSynonyms(string);
void clearMemory();
void benjaminiHochberg();
void createCounters();
bool pcomparison(pNode * n1, pNode * n2);
bool icomparison(pNode * n1, pNode * n2);



int main(int argc, char* argv[])
{
	clock_t start,stop;
	if (argc<ARGUMENTS)
	{
		cout << "Not enough arguments\n";
		exit(1);
	}
	
	thread_count=atoi(argv[6]);
	thread *t= new thread[thread_count];
	species="human";
	iterations=atoi(argv[5]);
    iterations-=iterations % thread_count;
	map_all=new bits[iterations];
	cout << "Reading GO category data" << endl;
	getClean(argv[4]);
	cout << "Reading frequent itemsets" << endl;
	getFrequent(argv[8]);
	cout << "Reading associations" << endl;
	getAssociations(argv[9]);
	cout << "Reading synonym data" << endl;
	getSynonyms(argv[7]);
	cout << "Reading interaction data" << endl;
	getInteractions(argv[1]);

	
	cout << "Synonym matching" << endl;
	fixInteractions();
	cout<< "Calculating query GO overlap" << endl;
	start=clock();
	getMirnas(argv[3]);
	cout << "Getting Random miRNA groups" << endl;
	/*
     * Spawn multiple threads to calculate unions
     */
	if (thread_count>1)
	{
		if ((group_size==1) && (iterations< finalInteractions.size()))
		{
			getRandom(group_size,0,1);
			// getRandomSame(group_size,0,1,argv[15]);
		}
		else
		{
			for (int u=0; u<thread_count; u++)
				t[u]=thread(getRandom,group_size,u,thread_count);
				// t[u]=thread(getRandomSame,group_size,u,thread_count,argv[15]);
		
			for (int u=0; u<thread_count; u++)
				t[u].join();
		}
	}
	else
		getRandom(group_size,0,1);
	calculateCounts();\

	createCounters();
	cout << "Getting GO overlap for " << iterations << " random miRNA groups" << endl;
	/*
     * Spawn multiple threads to calculate intersections
     */


	/*
	 * Find intersection for single-gene itemsets
	 */
    if (thread_count>1)
	{
		for (int u=0; u<thread_count; u++)
			t[u]=thread(frequentIntersectionsSingle,u,thread_count);
		
		for (int u=0; u<thread_count; u++)
			t[u].join();
	}
	else
		frequentIntersectionsSingle(0,1);

	/*
	 * Find intersection for multi-gene itemsets
	 */
	if (thread_count>1)
	{
		for (int u=0; u<thread_count; u++)
			t[u]=thread(frequentIntersectionsMulti,u,thread_count);
		
		for (int u=0; u<thread_count; u++)
			t[u].join();
	}
	else
		frequentIntersectionsMulti(0,1);

	/*
	 * Use index and calculate final intersections
	 */
	if (thread_count>1)
	{
		for (int u=0; u<thread_count; u++)
			t[u]=thread(findIntersections,u,thread_count);
		
		for (int u=0; u<thread_count; u++)
			t[u].join();
	}
	else
		findIntersections(0,1);

	cout << "Making Benjamini-Hochberg corrections" << endl;
	benjaminiHochberg();
	cout << "Writing final output" << endl;
	writeOutput(argv[2]);

	return 0;
}

/*
 * The following function reads the differentially expressed miRNAs 
 * in the file provided by the user
 *
 * @param filename: the input file specified by the user
 */
void getMirnas(string filename)
{
	
	string line,miRNA;
	int notFound=0,found=0,intersection;
	float target_genes;
	bitset<BSIZE> gene_map;
	/*
     * Read file
     */
	ifstream inFile;
	inFile.open(filename);

	if (inFile.is_open())
	{

		while(getline(inFile,line))
		{
			miRNA=trim(line);
			/*
             * Calculate genes targeted by the given miRNAs
             * 
             * If they do not exist in our hash table
             * increase the counter to be printed as soon as the
             * process ends
             */
            
			if (miRNA=="") continue;
			if (miRNAs.find(miRNA)!=miRNAs.end())
			{
				found++;
				gene_map |= (*finalInteractions[miRNAs[miRNA]]);

			}
			else
			{
				notFound++;
			}
		}

		inFile.close();
		cout << "Found " << found << " differentially expressed miRNAs" << endl;
		if (notFound)
			cout << notFound << " miRNAs were not found in the set of interactions" << endl;
		group_size=found;
		total_input=found+notFound;
		target_genes=gene_map.count();
		/*
         * Calculate GO categories associated with given miRNAs.
         * Those are the candidates for which to calculate an empirical p-value.
         * The rest have a pvalue=1.
         *
         * Add the results to a special vector
         */
		for (goGenes_type::iterator git=goGenes.begin(); git!=goGenes.end(); git++)
		{	
			/*
             * If the miRNA set has any common genes with a go category, then 
             * add it to the list as a candidate or else add it to a special list
             */
            int total=git->second->size();
            intersection=0;
			string category=git->first;
			int goFreqSize=0;

			if (goFreq.find(category)!=goFreq.end())
			{
				vector <int> * fmap=goFreq[category];
				int totalf=fmap->size();
				for (int i=0; i<totalf; i++)
				{
					vector <int> * gfmap=freqGenes[(*fmap)[i]];
					int totalgf=gfmap->size();
					for (int j=0; j<totalgf; j++)
					{
						if (gene_map[(*(gfmap))[j]]==1)
							intersection++;
					}
				}
				goFreqSize=totalf;
			}
			
			/*
			 * Not executed for empty vectors
			 * (categories consisting only of itemsets)
			 */
			for (int k=0; k<total; k++)
			{
				if (gene_map[(*(git->second))[k]]==1)
					intersection++;
			}

			if (intersection>0)
			{
				goCatNode newGO;
				
				newGO.category=git->first;
				newGO.go_size= git->second->size() + goFreqSize;
				newGO.intersection=intersection;
				newGO.overlap_proportion= intersection/target_genes;
				newGO.name=goNames[git->first];
				checkGO.push_back(newGO);
				pvalues.push_back(0);
				GOcount++;

			}
			else
			{
				goCatNode newGO;
				newGO.category=git->first;
				newGO.go_size= git->second->size() + goFreqSize;
				newGO.intersection=0;
				newGO.overlap_proportion=0.0;
				newGO.name=goNames[git->first];
				noCheckGO.push_back(newGO);
			}

		}
	}

}

/*
 * This function reads interactions from a file provided by the user in the alternative form
 *
 * It saves the interactions in a temporary hash table of linked lists
 * @param filename : the file name as provided by the user
 */
void getInteractions(string filename)
{
	
	string line,gline;
	ifstream inFile;
	int index;
	
	inFile.open(filename);
	
	if (inFile.is_open())
	{
		while (getline(inFile,line))
		{
			string miRNA, gene;
			token_list tokens,gtokens;

			line=trim(line);
			
			index=line.find_first_of("|");
			miRNA=line.substr(0,index);
			gene=line.substr(index+1);
			
			/* 
			 * If miRNA does not exist in the respective hash tables
			 * assign it an internal id and add it
			 */			
			
			if (miRNAs.find(miRNA) == miRNAs.end())
			{
				miRNAs[miRNA]=miRNA_count++;
			}
			if (interactions.find(miRNAs[miRNA])==interactions.end())
			{
				interactions[miRNAs[miRNA]]=nullptr;
			}
			
			/*
			 * Add interaction to temporary list
			 */
			
			interaction * newInteraction= new interaction();

			newInteraction->gene=gene;
			newInteraction->next=interactions[miRNAs[miRNA]];
			interactions[miRNAs[miRNA]]=newInteraction;
		}
		inFile.close();
		
	}
}

/*
 * This function reads GO-gene associations as given by the user in an alternative form
 * Genes that do not exist in the interactions provided by the user are 
 * assigned an internal id and added to the gene hash table 
 *
 * @param filename: filename provided by the user
 */
void getClean(string filename)
{
	
	string line;
	ifstream inFile;
	goTemp go_tmp;
	inFile.open(filename);
	unordered_map<string,unordered_set<int> *> tempGO;
	
	/*
     * Temporary saving of data in a hash table containing linked lists
     * like we did for the interactions
     */
	if (inFile.is_open())
	{
		string category, gene, name,domain;
		int index;

		while (getline(inFile,line))
		{
			index=0;
			
			line=trim(line);
			
			index=line.find_first_of("|");
			gene=line.substr(0,index);
			line=line.substr(index+1);
			index=line.find_first_of("|");
			category=line.substr(0,index);
			line=line.substr(index+1);
			name=line;
			
			

			if (genes.find(gene) == genes.end())
			{
				genes[gene]=gene_count++;
			}
			
			/*
			 * Add name to GO names
			 */
			goNames[category]=name;
			
			/*
			 * If category does not exist add it
			 */ 
			if (tempGO.find(category)==tempGO.end())
			{
				tempGO[category]= new unordered_set<int>;
			}
			
			tempGO[category]->insert(genes[gene]);
		}
		inFile.close();
		
		/*
		 * Create a vector for each category's genes
		 */
		for (unordered_map<string,unordered_set<int> *>::iterator it=tempGO.begin(); it!=tempGO.end(); it++)
		{
			if (goGenes.find(it->first)==goGenes.end())
			{
					goGenes[it->first]= new vector <int>;
			}
			for (unordered_set<int>::iterator git=(it->second)->begin(); git!=(it->second)->end();git++)
			{	
			
				goGenes[it->first]->push_back((*git));
			}
			delete it->second;
		}
	}
}

void getFrequent(string filename)
{
	
	string line;
	ifstream inFile;
	inFile.open(filename);
	unordered_map<int,unordered_set<int> *> tempSet;
	
	/*
     * Temporary saving of data in a hash table containing linked lists
     * like we did for the interactions
     */
	if (inFile.is_open())
	{
		string set, gene, name,domain;
		int index;

		while (getline(inFile,line))
		{
			index=0;
			
			line=trim(line);
			
			index=line.find_first_of("|");
			gene=line.substr(0,index);
			set=line.substr(index+1);
			
			//cout << set << " " << gene << endl;
			if (genes.find(gene) == genes.end())
			{
				genes[gene]=gene_count++;
			}
			if (frequent.find(set) == frequent.end())
			{
				frequent[set]=freq_count++;
			}
			int id=frequent[set];
			/*
			 * If category does not exist add it
			 */ 
			if (tempSet.find(id)==tempSet.end())
			{
				tempSet[id]= new unordered_set<int>;
			}
			
			tempSet[id]->insert(genes[gene]);
		}
		inFile.close();
		
		/*
		 * Create a vector for each category's genes
		 */
		for (unordered_map<int,unordered_set<int> *>::iterator it=tempSet.begin(); it!=tempSet.end(); it++)
		{
			if (freqGenes.find(it->first)==freqGenes.end())
			{
				freqGenes[it->first]= new vector <int>;
			}
			for (unordered_set<int>::iterator git=(it->second)->begin(); git!=(it->second)->end();git++)
			{	
				freqGenes[it->first]->push_back((*git));
			}
			delete it->second;
		}

		for (itemset_type::iterator it=freqGenes.begin(); it!=freqGenes.end(); it++)
		{
			if (it->second->size()==1)
			{	
				//gene->freqset id
				uGenesAssoc[(*it->second)[0]]=it->first;
				if ((it->first)>totalu)
				{
					totalu=it->first;
				}
			}
			else
			{
				if ((it->first)>totalnu)
				{
					totalnu=it->first;
				}
			}
		}
		//increasing gives total num of elements
		totalu++;
		totalnu++;
		//fix associations
		for (int i=totalu; i<totalnu; i++)
		{
			vector<int> * fmap= freqGenes[i];
			freqLargeAssoc[i]=new vector <int>;
			int totalgenes=fmap->size();
			for (int k=0; k<totalgenes; k++)
			{
				freqLargeAssoc[i]->push_back(uGenesAssoc[(*fmap)[k]]);
			}
		}
	}
}

void getAssociations(string filename)
{
	
	string line;
	ifstream inFile;
	inFile.open(filename);
	unordered_map<string,unordered_set<int> *> tempGO;
	
	/*
     * Temporary saving of data in a hash table containing linked lists
     * like we did for the interactions
     */
	if (inFile.is_open())
	{
		string path,freq;
		int index;

		while (getline(inFile,line))
		{
			index=0;
			
			line=trim(line);
			
			index=line.find_first_of("|");
			path=line.substr(0,index);
			freq=line.substr(index+1);

			/* 
			 * if category only consists of itemsets
			 * add a vector for this category
			 */

			if (goGenes.find(path)==goGenes.end())
			{
				goGenes[path]=new vector <int>;
			}	
			
			if (frequent.find(freq)==frequent.end())
			{
				cout << "There was a problem with the itemset data. Aborting" << endl;
				exit(2);
			}
			if (tempGO.find(path)==tempGO.end())
			{
				tempGO[path]= new unordered_set<int>;
			}
			
			tempGO[path]->insert(frequent[freq]);

		}
		inFile.close();

		for (unordered_map<string,unordered_set<int> *>::iterator it=tempGO.begin(); it!=tempGO.end(); it++)
		{
			if (goFreq.find(it->first)==goFreq.end())
			{
					goFreq[it->first]= new vector <int>;
			}
			for (unordered_set<int>::iterator git=(it->second)->begin(); git!=(it->second)->end();git++)
			{	
				goFreq[it->first]->push_back((*git));
			}
			delete it->second;
		}

		
	}
}

/*
 * This function creates a list of synonyms for each gene name
 * 
 * @param filename : filename provided by the user
 */
void getSynonyms(string filename)
{
	ifstream inFile;
	string line,gene, taxid, taxonomy, synline,synonym;
	int index=0;

	taxid=(species=="human") ? "9606" : "10090";

	inFile.open(filename);
	if (inFile.is_open())
	{
		getline(inFile,line);

		while(getline(inFile,line))
		{
			line=trim(line);
			index=line.find_first_of("\t");
			taxonomy=line.substr(0,index);
			if (taxonomy!=taxid)
				continue;
			line=line.substr(index+1);
			index=line.find_first_of("\t");
			line=line.substr(index+1);
			index=line.find_first_of("\t");
			gene=line.substr(0,index);
			line=line.substr(index+1);
			index=line.find_first_of("\t");
			line=line.substr(index+1);
			index=line.find_first_of("\t");
			synline=line.substr(0,index);
			
			
			if (synline!="-")
			{
				string_list synonyms_tmp=new vector<string>;

				synonyms_tmp->push_back(trim(gene));
				while((index=synline.find_first_of("|"))!=-1)
				{
					synonym=synline.substr(0,index);
					synline=synline.substr(index+1);
					synonyms_tmp->push_back(trim(synonym));
				}
				/*
				 * Do not forget to add the last synonym
				 */
				synonym=synline;


				synonyms_tmp->push_back(trim(synonym));
			

				for (int k=0; k<synonyms_tmp->size(); k++)
				{
					if (synonyms.count((*synonyms_tmp)[k])==0)
						synonyms[(*synonyms_tmp)[k]]=synonyms_tmp;
			
				}
			}

		}
	}

}

/*
 * This function creates bitsets for each miRNA in the set of interactions and does the basic gene matching
 * of the original empiricalGO.py script
 */
void fixInteractions()
{
	for (geneInteractions::iterator it=interactions.begin(); it!=interactions.end(); it++)
	{

		for (interaction * oldInt=it->second;oldInt!=nullptr; oldInt=oldInt->next)
		{
			string gene=oldInt->gene;
			if (genes.find(gene)==genes.end())
			{
				vector<string> * alternatives;
				if (synonyms.find(gene)!=synonyms.end())
				{
					alternatives=synonyms[gene];
				}
				else
				{
					alternatives=nullptr;
				}
				if (alternatives!=nullptr)
				{
					for (int k=0; k<alternatives->size(); k++)
					{
						if (genes.find((*alternatives)[k])!=genes.end())
						{
							oldInt->gene=(*alternatives)[k];
							break;
						}
					}
				}
			}
		}
	}
	for (geneInteractions::iterator it=interactions.begin(); it!=interactions.end(); it++)
	{
		
		finalInteractions[it->first]= new bitset<BSIZE>;

		for (interaction * oldInt=it->second;oldInt!=nullptr; oldInt=oldInt->next)
		{
			string gene=oldInt->gene;
			if (genes.find(gene)==genes.end())
			{
				genes[gene]=gene_count++;
			}
			(*finalInteractions[it->first])[genes[gene]]=1;
		}
	}

}

/*
 * Get a number of random miRNA sets of size m and calculate their interactions 
 * using bitwise operations.
 * The number of random groups is specified by the "iterations" variable
 * The result  will be saved in a vector of bitsets
 * Internal IDs are used
 *
 * @param size: size of the random miRNA sets
 * @param t_num: the thread number
 * @param inc : increment step (the number of threads to be used)
 */

void getRandom(int size,int t_num, int inc)
{
	/*
     * initialize random number generator
     */
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> randMirna(1,miRNA_count-1);

	long int total_interactions=finalInteractions.size();

	if ((size==1) && (total_interactions < iterations))
	{
		int j=0;
		cout << "You have selected " << iterations << " iterations, but your query contains only 1 miRNA. " 
									 << total_interactions << " iterations will be performed." << endl;
		iterations=total_interactions;
		for (final_interactions_type::iterator fit=finalInteractions.begin(); fit!=finalInteractions.end(); fit++, j++)
		{
			map_all[j]=fit->second;
		}

	}
	else
	{
		for (int i=t_num; i<iterations; i+=inc)
		{
		
			/*
         	 * Get n random internal IDs where n=size
         	 * 
         	 * Bitwise OR for to calculate targeted genes
         	 */
        	bits gene_map=new bitset<BSIZE>;
		
			for (int j=0; j<size; j++)
			{
				int id;
				id=randMirna(gen);
				(*gene_map) |=(*finalInteractions[id]);
			}
			map_all[i]=gene_map;
		}
	}

}


void frequentIntersectionsSingle(int t_num,int inc)
{
	for (int git=t_num; git<totalu; git+=inc)
	{
		freqCounters[git].charPointer=new unsigned char[iterations];

		int geneid=(*freqGenes[git])[0];
		for (int i=0; i<iterations; i++)
		{
			int intersection=0;
			bits gene_map=map_all[i];

			if ((*gene_map)[geneid]==1)
			{
				intersection++;
			}

			freqCounters[git].charPointer[i]=intersection;
		}
	}
	//cout << totalu << " " << totalnu << endl;
	
}

void frequentIntersectionsMulti(int t_num,int inc)
{
	unsigned char * intermediateChar=new unsigned char[iterations];
	int * intermediateInt = new int[iterations];

	for (int fit=totalu+t_num; fit<totalnu; fit+=inc)
	{
		
		//cout << fit << " " << freqGenes[fit]->size() << endl;
		vector <int> * fmap=freqLargeAssoc[fit];
		int totalgenes=fmap->size();

		if (totalgenes<=255)
		{
			memset(intermediateChar,0,iterations*sizeof(unsigned char));

			for (int k=0; k<totalgenes; k++)
			{
				// int intersection=0;
				// bits gene_map=map_all[i];
				int geneid=(*fmap)[k];

				for (int i=0; i<iterations; i++)
				{
				
					intermediateChar[i]+=freqCounters[geneid].charPointer[i];

				}

			}

			// for (int i=0; i<iterations; i++)
			// {
			// 	freqCounters[fit][i]=intermediate[i];
			// }
			freqCounters[fit].charPointer=new unsigned char[iterations];
			freqCounters[fit].intPointer=nullptr;
			memcpy(freqCounters[fit].charPointer, intermediateChar, iterations*sizeof(unsigned char));
		}
		else
		{
			memset(intermediateInt,0,iterations*sizeof(int) );

			for (int k=0; k<totalgenes; k++)
			{
				// int intersection=0;
				// bits gene_map=map_all[i];
				int geneid=(*fmap)[k];

				for (int i=0; i<iterations; i++)
				{
				
					intermediateInt[i]+=freqCounters[geneid].charPointer[i];

				}

			}

			// for (int i=0; i<iterations; i++)
			// {
			// 	freqCounters[fit][i]=intermediate[i];
			// }
			freqCounters[fit].intPointer=new int[iterations];
			freqCounters[fit].charPointer=nullptr;
			memcpy(freqCounters[fit].intPointer, intermediateInt, iterations*sizeof(int));
		

		}

		
		
	}
}

/*
 * This function calculates the intsections for all random miRNA sets
 * for the candidate GO categories and calculate the sets with greater overlap
 * than the queried one
 *
 * @param t_num: the thread number
 * @param inc : increment step (the number of threads to be used)
 */
void findIntersections(int t_num,int inc)
{
	bits gene_map;
	vec go_map;

	int parameter=3;
	/*
	 * Depending on the thread select specific GO categories to check
	 * without overlap between cores
	 */
	long int total=checkGO.size();
	int * intersections=new int[iterations]();
	for (int git=t_num; git<total; git+=inc)
	{	
		vector <int> * freq_map;
		int totalf; 
		int total_k;
		string category;
		int i,k,j;
		double overlap,baseline_overlap;

		category=checkGO[git].category;
		baseline_overlap=checkGO[git].overlap_proportion;
		
		if (goFreq.find(category)!=goFreq.end())
		{
			freq_map = goFreq[category];
			totalf=goFreq[category]->size();
		}
		else
		{
			totalf=0;
		}

		go_map=goGenes[category];
		total_k=go_map->size();

		memset(intersections,0,iterations*sizeof(int));
		
		// cout << "totalf=" << totalf << endl;
		// cout << "totalf=" << totalf << " " << category << endl; 
		for (k=0; k<totalf; k++)
		{
			

			if (freqCounters[(*freq_map)[k]].charPointer!=nullptr)
			{
				unsigned char * counterMap = freqCounters[(*freq_map)[k]].charPointer;
				for (i=0; i<iterations; i++)
				{	
					
					intersections[i]+=counterMap[i];
				}
			}
			else
			{
				int * counterMap = freqCounters[(*freq_map)[k]].intPointer;
				for (i=0; i<iterations; i++)
				{	
					
					intersections[i]+=counterMap[i];
				}
			}
			
		}

		for (i=0; i<iterations; i++)
		{	
			int intersection=intersections[i];

			gene_map=map_all[i];

			for (k=0; k<total_k; k++)
			{ 
				if ((*gene_map)[(*go_map)[k]]==1)
					intersection++;
			} 

			overlap=intersection/i_counts[i];
			//checkGO[git].mean_overlap+=overlap;

			if (overlap >= baseline_overlap)
			{	
				pvalues[git]++;
			}
		}
	}
}

/*
 * Write final output to a file
 *
 * @param filename: filename provided by the user
 */
void writeOutput(string filename)
{
	ofstream outFile;
	int total_check=checkGO.size(), total_n_check=noCheckGO.size();

	outFile.open(filename);

	if (outFile.is_open())
	{
		/*
		 * File header
		 */
		outFile << "GO-term-ID\tGO-term-size\tObserved-Target-Gene-Overlap-Proportio\t";
		outFile << "Mean-Random-Simulated-MicroRNA-Target-Overlap-Proportion\tOne-sided-empirical-p-value\tBenjamini-Hochberg-0.05-FDR" << endl;
		for (int i=0; i< total_check ; i++)
		{
			double pval=fdrList[i]->pvalue;
			outFile << checkGO[i].category << "~" << checkGO[i].name << "\t";
			outFile << checkGO[i].go_size << "\t" << checkGO[i].overlap_proportion << "\t";
			outFile << checkGO[i].mean_overlap/iterations << "\t";
			//if (pval< 0.001) outFile << scientific; 
			outFile << pval  << "\t" << fdrList[i]->fdr << endl;
		}
		
		for (int i=0; i < total_n_check ; i++)
		{
			double pval=fdrList[i+total_check]->pvalue;
			outFile << noCheckGO[i].category << "~" << noCheckGO[i].name << "\t";
			outFile << noCheckGO[i].go_size << "\t" << noCheckGO[i].overlap_proportion << "\t";
			outFile << noCheckGO[i].mean_overlap/iterations << "\t";
			outFile << pval  << "\t" << fdrList[i+total_check]->fdr << endl;
		}
		outFile.close();
	}
}

/*
 * This function calculates the number of genes targeted by each random miRNA group
 */
void calculateCounts()
{
	for (int i=0; i< iterations; i++)
	{
		i_counts.push_back((float)(map_all[i])->count());
	}
}


/*
 * This function calculates the BH FDR correction
 */
void benjaminiHochberg()
{
	int total_p=checkGO.size(), total_no_p=noCheckGO.size();
	int k,i,maxi=0;
	double star=0.05, doubleStar=0.01;
	
	for (k=0; k<total_p; k++)
	{
		pNode * tmp=new pNode;
		
		tmp->pvalue=(double)pvalues[k]/iterations;
		tmp->index=k;
		fdrList.push_back(tmp);

	}
	
	for (i=0; i<total_no_p; i++)
	{
		pNode * tmp=new pNode;
		
		tmp->pvalue=1.0;
		tmp->index=i+k;
		fdrList.push_back(tmp);
	}

	sort(fdrList.begin(), fdrList.end(), pcomparison);

	int total_fdr=fdrList.size();

	for (int j=0; j<total_fdr; j++)
	{
		double starcheck;
		i=j+1;
		starcheck=(i*star)/total_fdr;
		if (fdrList[j]->pvalue < starcheck)
		{
			maxi=j;
		}
	}
	if (maxi>total_fdr)
	{
		maxi=total_fdr;
	}
	for (int j=0; j<maxi; j++)
	{
		fdrList[j]->fdr="*";
	}

	
	for (int j=0; j<total_fdr; j++)
	{
		double starcheck;
		i=j+1;
		starcheck=(i*doubleStar)/total_fdr;
		if (fdrList[j]->pvalue < starcheck)
		{
			maxi=j;
		}
	}
	if (maxi>total_fdr)
	{
		maxi=total_fdr;
	}
	for (int j=0; j<maxi; j++)
	{
		fdrList[j]->fdr="**";
	}

	sort(fdrList.begin(), fdrList.end(), icomparison);

}

/*
 * Auxiliary functions
 */
string trim(string mystr)
{
	int start,stop;

	start=mystr.find_first_not_of("\n\t ");
	stop=mystr.find_last_not_of("\n\t ");

	return mystr.substr(start,stop-start+1);

}

string trim_chars_left(string mystr,string chars)
{
	int start;

	start=mystr.find_first_not_of(chars);

	return mystr.substr(start);
}

bool pcomparison(pNode * n1, pNode * n2)
{
	return (n1->pvalue < n2->pvalue);
}

bool icomparison(pNode * n1, pNode * n2)
{
	return (n1->index < n2->index);
}
void createCounters()
{
	int total=freqGenes.size();
	freqCounters=new counterType [total];
	// for (int i=0; i<total; i++)
	// {
		// freqCounters[i]=new counterType [iterations];
		// for (int j=0; j<iterations; j++)
		// {
		// 	freqCounters[i][j]=0;
		// }
	// }
}