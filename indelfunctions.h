
// indel data structures and functions for split-read indel alignment and detection  //
#ifndef INC_indelfunctions_H
#define INC_indelfunctions_H

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "bamread.h"

extern int MRLENGTH;
extern int LOWQV;
extern int COMPLEX_INDELS;
extern int MIN_READS_FOR_INDEL;
extern int TARGETED; // targeted sequencing 
// "SPINDEL" conflicts with SPINDEL of variant calling 

// data structure for storing information about indels
typedef struct SPINDEL
{
	char code;  char partial; // gapped alignment, split-read alignment, de novo assembly 
	int chrom; int position; int length; char strand; int ambiguity; int firstbase,lastbase;
	int l1; 
	uint32_t cigar[2]; // 4D 20I | 10I 
	short readsf; short readsr; short reads; short newreads; int cluster; // for storing information about the cluster of indels
	char* insertedseq;  
	//int L,R; int match; int mismatch; int ambiguity;
	//char* readid; char* sequence; //char* cigar; 
	//int matechrom; int mateposition; char matestrand; 
} SPINDEL; 

// also determine the rightmost start position for each indel in addition to the leftmost start position
int leftjustifyindel(SPINDEL* indel,char* sequence,int l,REFLIST* reflist);
int compare_indel(const void* a,const void* b);
int parse_gapped_read(struct alignedread* read,REFLIST* reflist,int current,SPINDEL* indelreadlist,int* indelreads,char code);
void cluster_indelreads(SPINDEL* indelreadlist,int indelreads,REFLIST* reflist,int partial);
void print_clusters(SPINDEL* indelreadlist,int indelreads,REFLIST* reflist,int partial,FILE* outfp);

#endif


