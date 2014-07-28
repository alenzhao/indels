#include "kmertable.h"

// functions and data structures for kmertable (index of 14-mers in the reference genome)
// this is used for split-read mapping as well as realignment of indel spanning reads 

int MAXMATCHES = 1000; // max size of kmer match list -> should depend on read length but constant for now
int MAX_KMER_COUNT = 64; // k-mers more frequent than this are not used...
int MAX_KMER_HITS = 5; 

// bug detected may 17 2012, BTI['S'] should be 0 and not '4'
unsigned int BTI[] = {
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  3, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  3, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};

unsigned int BTIoriginal[] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 0, 4, 4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 0, 4, 4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};



unsigned int kmertoint(char* seq,int l)
{
	int i=0; unsigned int hash =BTI[(int)seq[0]];
	for (i=1;i<l;i++) { hash <<=2; hash += BTI[(int)seq[i]]; } 
	return hash;
}

// diagonal start position
int compare_hits_diagonalpos(const void *pa, const void *pb)
{
        const MATCH *a = (MATCH*) pa; const MATCH *b = (MATCH*) pb;
	if (a->pos == b->pos) return a->first -b->first;
	else return a->pos -b->pos; 
}

int compare_hits_querypos(const void *pa, const void *pb)
{
        const MATCH *a = (MATCH*) pa; const MATCH *b = (MATCH*) pb;
	if (a->first == b->first) 
	{
		return a->last -b->last;
	}
	else return a->first -b->first; 
}

int compare_hits(const void *pa, const void *pb)
{
	const MATCH *a = (MATCH*) pa; const MATCH *b = (MATCH*) pb;
	if (a->mm == b->mm)
	{
		if (a->ch == b->ch) 
		{
			if (a->pos == b->pos) return a->strand-b->strand;
			else return a->pos-b->pos;
		}
		else return a->ch-b->ch;
	}
	else return a->mm-b->mm;
}


void build_kmertable_notworking(REFLIST* reflist,int kmer,KMERTABLE* kmertable,int gap)
{
	int s=0;
	for (s=0;s<reflist->ns;s++)
	{
		if (s < 5) fprintf(stderr,"building table for  %d %s %d \n",s,reflist->names[s],reflist->lengths[s]);
		build_kmertable(reflist,kmer,kmertable,gap,s,0,reflist->lengths[s]);
	}
}


int build_kmertable(REFLIST* reflist,int kmer,KMERTABLE* kmertable,int gap, int chromosome,int start, int end)
{
	int s=chromosome,current =0,i=0,j=0; unsigned int seq=0,seq1=0,seq2=0;  int initialize =1;
	int kmin = kmer; if (kmer%2==1) kmin--;
	int kplus = kmer; if (kmer%2==1) kplus++;
	// if kmer=13, kmer/2 = 6 so first half is 6bp second half is 7bp 
	for(i=start;i< end-kmer-gap;i++)
	{
		if (initialize ==1)
		{
			seq1=0; seq2=0;
			for (j=0;j<kmer/2;j++)
			{
				seq1 <<=2; seq1 += BTI[(int)reflist->sequences[s][i+j]]; 
			}
			for  (j=kmer/2;j<kmer;j++) { seq2 <<=2; seq2 += BTI[(int)reflist->sequences[s][i+j+gap]]; }
			initialize = 0;
		}
		else
		{
			seq1 -= (BTI[reflist->sequences[s][i-1]]<<(kmin-2)); seq1 <<=2; seq1 += BTI[reflist->sequences[s][i+kmer/2-1]];
			seq2 -= (BTI[reflist->sequences[s][i+kmer/2+gap-1]]<<(kplus-2)); seq2 <<=2;  seq2 += BTI[reflist->sequences[s][i+kmer+gap-1]];
		}

		seq = (seq1<<(kplus)) + seq2;
		//printf("%d %d %d %d\n",seq1,seq2,i,seq);
		if (kmertable->counts[seq] < MAXMATCHES)  // only add for counts below max value of short
		{
			kmertable->locarray[current][2] = kmertable->index[seq];
			kmertable->locarray[current][0] = s; kmertable->locarray[current][1] = i;
			kmertable->index[seq] = current; kmertable->counts[seq]++;
			current++;
		}
	}
	return 1;
}



int init_kmertable(REFLIST* reflist,int kmer,KMERTABLE* kmertable,int genomelength)
{
	int s=0; int tablesize = 1<<(2*kmer);
	unsigned long memallocated = 0;
	kmertable->index = (int*)malloc(sizeof(int)*tablesize); memallocated += 4*(unsigned long)tablesize;
	kmertable->counts = (short*)malloc(sizeof(short)*tablesize); memallocated += 2*(unsigned long)tablesize;
	kmertable->locarray = (int**)malloc(sizeof(int*)*genomelength); memallocated += 16*genomelength;
	for (s=0;s<genomelength;s++) kmertable->locarray[s] = (int*)malloc(sizeof(int)*3);
	fprintf(stderr,"table size %d maxkmers %d MEMory allocated %ld %d\n",tablesize,genomelength,memallocated,MAX_KMER_HITS);
	return 1;

}



// build kmertable from fasta file
int build_kmertable_full(REFLIST* reflist,int kmer,KMERTABLE* kmertable,int gap)
{
	int s=0,current =0,i=0,j=0; unsigned int seq=0,seq1=0,seq2=0;  int kmersadded=0;
	for (s=0;s<reflist->ns;s++)
	{
		if (s < 3) fprintf(stderr,"ref %d %s %d \n",s,reflist->names[s],reflist->lengths[s]);
		i=0; seq1=0; seq2=0;
		for (j=0;j<kmer/2;j++)
		{
			seq1 <<=2; seq1 += BTI[reflist->sequences[s][i+j]];

		}
		for (j=0;j<kmer/2;j++) { seq2 <<=2; seq2 += BTI[reflist->sequences[s][i+kmer/2+gap+j]]; }

		while (i < reflist->lengths[s]-kmer-gap)
		{
			if (i > 0)
			{
				seq1 -= (BTI[reflist->sequences[s][i-1]]<<(kmer-2)); seq1 <<=2; seq1 += BTI[reflist->sequences[s][i+kmer/2-1]];
				seq2 -= (BTI[reflist->sequences[s][i+kmer/2+gap-1]]<<(kmer-2)); seq2 <<=2; seq2 += BTI[reflist->sequences[s][i+kmer+gap-1]];
			}
			seq = (seq1<<kmer) + seq2;
			//printf("%d %d %d %d\n",seq1,seq2,i,seq);
			if (kmertable->counts[seq] < 10000)  // only add for counts below max value of short
			{
				kmertable->locarray[current][2] = kmertable->index[seq];
				kmertable->locarray[current][0] = s;
				kmertable->locarray[current][1] = i;
				kmertable->index[seq] = current; kmertable->counts[seq]++;
				current++;
			}
			i++;
			kmersadded++;
		}
	}
	fprintf(stderr,"table build with %d kmers added \n",kmersadded);
	return 1; 
}


int build_hashtables(REFLIST* reflist,int kmer,KMERTABLE* kmertable,KMERTABLE* kmertable1,int ntables)
{
	int genomelength =0,s=0,tablesize = 1<<(2*kmer);
	fprintf(stderr,"table size %d \n",tablesize);
	for (s=0;s<reflist->ns;s++) genomelength += reflist->lengths[s]-kmer;

	kmertable->index = (int*)malloc(sizeof(int)*tablesize);
	kmertable->counts = (short*)malloc(sizeof(short)*tablesize);
	kmertable->locarray = (int**)malloc(sizeof(int*)*genomelength);
	for (s=0;s<genomelength;s++) kmertable->locarray[s] = (int*)malloc(sizeof(int)*3);
	build_kmertable_full(reflist,kmer,kmertable,0);

	if (ntables > 1)
	{
		kmertable1->index = (int*)malloc(sizeof(int)*tablesize);
		kmertable1->counts = (short*)malloc(sizeof(short)*tablesize);
		kmertable1->locarray = (int**)malloc(sizeof(int*)*genomelength);
		for (s=0;s<genomelength;s++) kmertable1->locarray[s] = (int*)malloc(sizeof(int)*3);
		build_kmertable_full(reflist,kmer,kmertable1,kmer/2);
	}
	return 1;
}


