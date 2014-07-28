// functions for storing, clustering indels identified from parsing a bam/sam file as well as split-read alignment 
// this are used in split-read alignment code
#include "indelfunctions.h"
//int MRLENGTH = 50; int LOWQV =2;
int COMPLEX_INDELS = 0;

#include "printindels.c"

int compare_indel(const void* a,const void* b) // for sorting list of indels...
{
	const SPINDEL *ia = (const  SPINDEL*)a; const  SPINDEL *ib = (const  SPINDEL*)b;
	if (ia->chrom == ib->chrom) 
	{
		if (ia->position == ib->position) 
		{
			if (ia->length == ib->length) return (int)ia->code -(int)ib->code;
			else return ia->length -ib->length;
		}
		else return ia->position -ib->position;
	}
	else return ia->chrom - ib->chrom;
}

// also add multiple indels on same read into separate record for realignment...

// insert indels into priority queue sorted in reverse order by position, new indels on top of heap 
int parse_gapped_read(struct alignedread* read,REFLIST* reflist,int current,SPINDEL* IndelList,int* indelreads,char code)
{
	int i=0,j=0; int ir = *indelreads; int l1=0,l2=0,op=0,ol=0,ol1=0;
	if (read->gaps > 1 && COMPLEX_INDELS >=2) // check if read has multiple DI cigar values 
	{
		combine_indelcigars(read,reflist); 
	}
	//fprintf(stderr," %c %d %d %s,%c,%s,%d %s %s\n",read->strand,l1,l2,read->readid,read->strand,read->matechrom,read->mateposition,read->sequence,read->cigar);
	// for DI consecutive cigars -> create indel haplotypes  TODO jan7 2013
	for (i=0;i<read->cigs;i++)
	{
                op = read->cigarlist[i]&0xf; ol =read->cigarlist[i]>>4;
		if (i < read->cigs-1 && op == BAM_CDEL && (read->cigarlist[i+1]&0xf) == BAM_CINS && COMPLEX_INDELS >= 0)  // block substitution 4D 20I
		{
			ol1 = read->cigarlist[i+1]>>4;
			IndelList[ir].strand = read->strand; IndelList[ir].chrom = current; 
			IndelList[ir].position = read->position+l2; IndelList[ir].length =  ol1-ol; 
			IndelList[ir].cigar[0] = read->cigarlist[i]; IndelList[ir].cigar[1] = read->cigarlist[i+1]; IndelList[ir].l1 = l1;
			IndelList[ir].code = code; IndelList[ir].partial = '0'; 
			IndelList[ir].insertedseq = calloc(sizeof(char),ol1+1); 
			for (j=0;j<ol1;j++) IndelList[ir].insertedseq[j] = read->sequence[l1+j];
			ir++; (*indelreads)++;
			l2 += ol; l1 += ol1;
			i++; // extra cigar 
		}
		else if (op == BAM_CMATCH || op == 8) { l1 += ol; l2 += ol; }  // 'M' or 'X'
		else if (op == BAM_CSOFT_CLIP) l1 += ol;
		else if (op == BAM_CDEL)
		{
			if (read->quality[l1] < LOWQV+QVoffset) { l2 += ol; continue; } 
			IndelList[ir].strand = read->strand; IndelList[ir].chrom = current; IndelList[ir].position = read->position+l2; IndelList[ir].length = -1*ol; 
			IndelList[ir].code = code; IndelList[ir].partial = '0'; 
			IndelList[ir].cigar[0] = read->cigarlist[i]; IndelList[ir].cigar[1] = 0; IndelList[ir].l1 = l1;
			IndelList[ir].firstbase = read->position+l2; IndelList[ir].lastbase = read->position+l2+ol; 
			leftjustifyindel(&IndelList[ir],read->sequence,read->readlength,reflist);
			ir++; (*indelreads)++; l2 += ol;
		}
		else if (op == BAM_CINS) 
		{
			// added extra filters for novoalign indels at very end of reads..
			if (i ==0 || i == read->cigs-1 || read->quality[l1] < LOWQV+QVoffset || (i ==1 && (read->cigarlist[0]&0xf)==BAM_CMATCH && (read->cigarlist[0]>>4) < 4) || (i==read->cigs-2 && (read->cigarlist[read->cigs-1]&0xf)==BAM_CMATCH && (read->cigarlist[read->cigs-1]>>4) < 4 )) 
			{ 
				l1 += ol; continue; 
			} 
			IndelList[ir].strand = read->strand; IndelList[ir].chrom = current; IndelList[ir].position = read->position+l2;  IndelList[ir].length = ol; 
			IndelList[ir].cigar[0] = read->cigarlist[i]; IndelList[ir].cigar[1] = 0; IndelList[ir].l1 = l1;
			IndelList[ir].code = code; IndelList[ir].partial = '0';
			IndelList[ir].firstbase = l1; IndelList[ir].lastbase = l1+ol;
			//IndelList[ir].matestrand =read->strand; IndelList[ir].matechrom = current;IndelList[ir].mateposition = read->mateposition;
			leftjustifyindel(&IndelList[ir],read->sequence,read->readlength,reflist);
			IndelList[ir].insertedseq = calloc(sizeof(char),ol+1); 
			for (j=0;j<ol;j++) IndelList[ir].insertedseq[j] = read->sequence[IndelList[ir].firstbase+j];

			ir++; (*indelreads)++;	l1 += ol;
		}
		else if (op == BAM_CREF_SKIP) l2 += ol;
		else if (op == BAM_CHARD_CLIP) l2 +=0; 
		else { read->cflag =1; break; }
	}
	return 1;
}
// also determine the rightmost start position for each indel in addition to the leftmost start position
int leftjustifyindel(SPINDEL* indel,char* sequence,int l,REFLIST* reflist)
{
	int t1,t2,t3,offset,lmax;
	if (indel->length < 0) // deletion
	{
		t1 = indel->lastbase-1; t2 = indel->firstbase-1; t3=0;
		while (t1 > 0 && t2 > 0)
		{
			if (reflist->sequences[indel->chrom][t1] == reflist->sequences[indel->chrom][t2]) { t1--; t2--; t3++;}
			else break;
		} offset = t3;

		t1 = indel->lastbase; t2 = indel->firstbase; t3=0; lmax = reflist->lengths[indel->chrom];
		while (t1 < lmax && t2 < lmax)
		{
			if (reflist->sequences[indel->chrom][t1] == reflist->sequences[indel->chrom][t2]) { t1++; t2++; t3++;}
			else break;
		} 
		indel->ambiguity = offset + t3; // distance between the leftmost start position and the rightmost start position
		indel->firstbase -= offset; indel->lastbase -=offset; indel->position -= offset;
	}
	else if (indel->length > 0)
	{
		t1 = indel->lastbase-1; t2 = indel->position-1; t3=0;
		//t1 = indel->lastbase-1; t2 = indel->firstbase-1; t3=0;
		while (t1-indel->length > 0 && t2 > 0)
		{
			if (sequence[t1] == reflist->sequences[indel->chrom][t2]) { t1--; t2--; t3++;}
			//if (sequence[t1] == sequence[t2]) { t1--; t2--; t3++;} 
			else break;
		}offset = t3;
		indel->firstbase -= offset; indel->lastbase -=offset; indel->position -= offset;

		t1 = indel->lastbase; t2 = indel->firstbase; t3=0; lmax =  reflist->lengths[indel->chrom];
		while (t1 < l && t2 < l)
		{
			if (sequence[t1] == sequence[t2]) { t1++; t2++; t3++;}
			else break;
		} 
		indel->ambiguity =  t3; // distance between the leftmost start position and the rightmost start position
	}
	//fprintf(stderr,"indel %d %d offset %d t3 %d \n",indel->position,indel->length,offset,t3);
	return 1;
}

