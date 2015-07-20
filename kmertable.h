#ifndef INC_kmertable_H
#define INC_kmertable_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "readfasta.h"

extern int MAXMATCHES;
extern int MAX_KMER_HITS;
extern int MAX_KMER_COUNT;

extern unsigned int BTI[];

// kmertable similar to SSAHA, kmer represented as 32bit integer, pointers also represented as integers...
//index is pointer to first occurence of k-mer in locarray, locarray is three-tuple (chrom,position,next)
//change int to int32_t 10/30/13
typedef struct 
{ 
	int* index; short* counts; int gap; int** locarray; int kmer; 
} KMERTABLE;

// information for hit (aligned) between a query (read) and the reference sequence 
typedef struct 
{ 
	// pos is diagonal coordinate for alignment used for clumping matches
	int ch; int pos; char strand; short offset;
	short first, last; // first and last are coordinates in read that are spanned by this hit without any gaps, last position is included
	int reffirst; int reflast; 
	short el,er; // bases extended to left and bases extended to right 
	char valid; 
	int mb; int mm; 
	//int indelpos; char type; char* bases; // indel information to which the read is aligned to
	uint32_t* cigarlist; int32_t cigs; int max_len; // for realloc
	float bestscore; int previous; int DL,IL;  // for keeping track of longest path in directed acyclic graph
	float clippenalty; // clipping penalty
} MATCH; 

unsigned int kmertoint(char* seq,int l);

int compare_hits(const void *pa, const void *pb);
int compare_hits_diagonalpos(const void *pa, const void *pb);
int compare_hits_querypos(const void *pa, const void *pb);

void build_kmertable_notworking(REFLIST* reflist,int kmer,KMERTABLE* kmertable,int gap);

int build_kmertable(REFLIST* reflist,int kmer,KMERTABLE* kmertable,int gap, int chromosome,int start, int end);

int init_kmertable(REFLIST* reflist,int kmer,KMERTABLE* kmertable,int genomelength);

// build kmertable from fasta file
int build_kmertable_full(REFLIST* reflist,int kmer,KMERTABLE* kmertable,int gap);

int build_hashtables(REFLIST* reflist,int kmer,KMERTABLE* kmertable,KMERTABLE* kmertable1,int ntables);

#endif
