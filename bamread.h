#ifndef INC_bamread_H
#define INC_bamread_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "../readfasta.h"
#include "sam.h"

extern int QVoffset;

extern char INT_CIGAROP[];

struct alignedread
{
	int mtid; // matetid
	
	char* readid; short flag; char strand; char matestrand; 
	char* chrom; int position; char* matechrom; int mateposition;
	int mquality; int IS; 
	int tid;
	char* sequence; char* quality; int readlength;
	uint32_t* cigarlist; int32_t cigs; 
	uint32_t* fcigarlist; int32_t fcigs;

	int l1,l2,mismatches; //int indels; // no of mismatches and no of insertions/deletions
	int alignedbases; int clipped, gaps, span;
	int cflag; 
	int XC; // XC value for clipped reads
	//int max_lengths[4]; // max length of sequence/quality, cigarlist, fcigarlist, readid
};

// data structure to hold information about realignment of read using known indel
struct REALIGNMENT
{
	int varid; // index to variantlist, variant that is included in new cigar
	uint32_t* cigarlist; int32_t cigs; int newpos;
	int added; int mismatches; int delta; 
};

int fetch_func(const bam1_t *b, void* data,struct alignedread* read);

void free_readmemory(struct alignedread* read);

int update_bam_record(bam1_t* b,int newpos,uint32_t* cigarlist,int32_t cigs,int flag);

void RC(char* sequence,char* revsequence,int l);

void read_stats(struct alignedread* read1,int* lq1,int* missing1);

int parse_cigar(struct alignedread* read,REFLIST* reflist,uint32_t* fcigarlist);

int combine_indelcigars(struct alignedread* read,REFLIST* reflist);

int modify_insertion_cigar(struct alignedread* read);

int validate_bam_header(char* bamfile,REFLIST* reflist);

#endif
