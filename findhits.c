#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>

int find_hits_forseed(char* sequence,int l,REFLIST* reflist,unsigned int hash, KMERTABLE* kmertable,MATCH* mlist,int* nmatches,char matestrand,int matepos,int shift,int maxcoord);
int read_firstpartmapped(char* sequence,int l,REFLIST* reflist,MATCH mlist[],int* clusters);
int read_secondpartmapped(char* sequence,int l,REFLIST* reflist,MATCH mlist[],int* clusters);
//int cluster_hits(char* sequence,int l,MATCH* mlist,int nmatches,REFLIST* reflist);


// extend these seeds until the number that match drops...
// AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA highly repetitive reads, low-complexity -> matches are in same region
// if the number of hits in the target region is too many, extend right single  base by base until number of hits drops
int extend_repetitive_hits(char* sequence,int l,REFLIST* reflist,MATCH* mlist,int first,int last,int* extend)
{
	int i=0,hits = last-first; int j=0, offset = mlist[first].last+1;  //int kmer = mlist[first].last-mlist[first].first+1;
	// offset is same for all the hits in the last = shift+kmer and is the first position after the original seed hit
	while (hits >= MAX_KMER_HITS && (j+ offset < l) )
	{ 	
		for (i=first;i<last;i++)
		{
			if (mlist[i].valid == '0') continue; 
			else if (mlist[i].reflast+1 < reflist->lengths[mlist[i].ch] && sequence[mlist[i].last+1] != 'N' && sequence[mlist[i].last+1] != reflist->sequences[mlist[i].ch][mlist[i].reflast+1]) 
			{
				mlist[i].valid = '0'; hits--; continue;
			}
			mlist[i].last++; mlist[i].reflast++;
		}j++;
		//fprintf(stdout,"j %d %d | ",j,hits);
	}
	*extend = j; // number of bases extended..

	/*
	// if still high and reached end of read, extend to left...
	j=mlist[first].first-1;	while (hits >= MAX_KMER_HITS && j >=0) 
	{
		for (i=first;i<last;i++)
		{
			if (mlist[i].valid == '0') continue; 
			else if (sequence[mlist[i].first-1] != 'N' && sequence[mlist[i].first-1] != reflist->sequences[mlist[i].ch][mlist[i].reffirst-1]) 
			{
				mlist[i].valid = '0'; hits--; continue;
			}
			mlist[i].first--; mlist[i].reffirst--;
		}j--;
	}*/
	if (hits >= MAX_KMER_HITS || hits ==0) return hits;
	fprintf(stdout,"offset %d trying to reduce number of hits %d by extension j:%d -> new_hits %d \n",mlist[first].first,last-first,j,hits); 
	// now shorten the list by removing the valid = 0 hits.. and extend their right position to last + j
	j = first; for (i=first;i<last;i++)
	{
		if (mlist[i].valid == '0') continue;
		mlist[j].ch = mlist[i].ch; mlist[j].pos = mlist[i].pos; mlist[j].offset = mlist[i].offset; 
		mlist[j].first = mlist[i].first; mlist[j].last = mlist[i].last; 
		mlist[j].reffirst = mlist[i].reffirst; mlist[j].reflast = mlist[i].reflast; 
		mlist[j].valid = '1'; j++;
	}
	return hits;
}

int find_hits_forseed(char* sequence,int l,REFLIST* reflist,unsigned int hash, KMERTABLE* kmertable,MATCH* mlist,int* nmatches,char matestrand,int matepos,int shift,int maxcoord)
{
	int first = kmertable->index[hash]; 
	if (kmertable->counts[hash] >= MAX_KMER_COUNT) return -2;  // too many hits for seed 
	else if (first < 0) return 0;
	int delta = -1*l; // related to insert size and if two ends can overlap

	int nm = *nmatches; // detect if kmer maps to low-complexity region with multiple matches that are shifted by 1-3 bp 
	int flag=0; int ch; int i=0,j=0;
	int reduced_hits=0;

	while (first > 0 && *nmatches < MAXMATCHES)
	{
		//fprintf(stdout,"%s %d %d %d %d \n",query,kmertable->locarray[first][1],matepos+offset,maxdistance,LR);
		// chromosome filter removed since it is local hashtable
		if ( (matestrand == '-' && matepos+shift -kmertable->locarray[first][1]< maxcoord && matepos+shift -kmertable->locarray[first][1] > delta) || (matestrand =='+' && kmertable->locarray[first][1]-matepos-shift > delta && kmertable->locarray[first][1]-matepos-shift < maxcoord)) flag = 1;
		else flag =0;

		if (flag ==1 && kmertable->gap ==0) 
		{
			mlist[(*nmatches)].ch = kmertable->locarray[first][0]; mlist[(*nmatches)].pos = kmertable->locarray[first][1]-shift;
			mlist[(*nmatches)].offset = shift;
			mlist[(*nmatches)].reffirst = kmertable->locarray[first][1]; mlist[(*nmatches)].reflast = kmertable->locarray[first][1]+kmertable->kmer+kmertable->gap-1;
			mlist[*nmatches].first = shift; mlist[*nmatches].last = shift + kmertable->kmer-1+kmertable->gap;  // last position is inclusive
			mlist[*nmatches].valid = '1';
			(*nmatches)++;
		}
		else if (flag ==1 && kmertable->gap > 0)
		{
			ch = kmertable->locarray[first][0]; 
			i=kmertable->kmer/2; 
			//if (kmertable->locarray[first][1] < 50 || kmertable->locarray[first][1] >= reflist->lengths[ch]-50) 
			//fprintf(stdout,"pos %d %d ch%d shift %d %d\n",kmertable->locarray[first][1],reflist->lengths[ch],ch,shift,l);
			while (sequence[shift+i] == reflist->sequences[ch][kmertable->locarray[first][1]+i] && i < kmertable->kmer/2+kmertable->gap) i++;
			j= kmertable->kmer/2+kmertable->gap; 
			while (sequence[shift+j-1] == reflist->sequences[ch][kmertable->locarray[first][1]+j-1] && j > kmertable->kmer/2) j--;
			if (j-i ==1) 
			{
				//fprintf(stdout,"gap match shift %d ref %d i,j %d %d %d\n",shift,kmertable->locarray[first][1],i,j,kmertable->kmer);
				mlist[(*nmatches)].ch = kmertable->locarray[first][0]; 
				mlist[(*nmatches)].pos = kmertable->locarray[first][1]-shift;
				mlist[(*nmatches)].offset = shift;
				mlist[(*nmatches)].reffirst = kmertable->locarray[first][1]; 
				mlist[(*nmatches)].reflast = kmertable->locarray[first][1]+i-1;
				mlist[*nmatches].first = shift; mlist[*nmatches].last = shift + i-1;
				mlist[*nmatches].valid = '1'; (*nmatches)++;

				mlist[(*nmatches)].ch = kmertable->locarray[first][0]; mlist[(*nmatches)].pos = kmertable->locarray[first][1]-shift;
				mlist[(*nmatches)].offset = shift+j;
				mlist[(*nmatches)].reffirst = kmertable->locarray[first][1]+j; mlist[(*nmatches)].reflast = kmertable->locarray[first][1]+kmertable->kmer+kmertable->gap-1;
				mlist[*nmatches].first = shift+j; mlist[*nmatches].last = shift + kmertable->kmer-1+kmertable->gap;
				mlist[*nmatches].valid = '1';
				(*nmatches)++;

				
			}
			//if (j-i ==1) fprintf(stdout,"GOOD "); else fprintf(stdout,"FAIL ");
		}		
		first = kmertable->locarray[first][2]; 
	}
	int extend =0;
	if ((*nmatches)-nm >= MAX_KMER_HITS && kmertable->gap ==0)
	{
		//fprintf(stdout,"reduce %d:%d | ",*nmatches,nm);
		reduced_hits =  extend_repetitive_hits(sequence,l,reflist,mlist,nm,*nmatches,&extend);
		*nmatches = nm+reduced_hits;
	}
	if ((*nmatches)-nm >= MAX_KMER_HITS || ((*nmatches)-nm >= 2*MAX_KMER_HITS && kmertable->gap > 0))
	{
		*nmatches = nm; 
		fprintf(stdout,"repetitve %d:%d %d | ",shift,kmertable->counts[hash],kmertable->gap); return -1; 
	}
	return kmertable->counts[hash];

}

// a read whose first half is mapped (70M30S cigar) : check if it spans a deletion/insertion at end
int read_firstpartmapped(char* sequence,int l,REFLIST* reflist,MATCH mlist[],int* clusters)
{
	char* bigstring = (char*)malloc(1024); char* smallstring = (char*)malloc(1024); char* pch;
	int j=0,t3=0,delta=0,offset=0,wts=40; // changed this value
	int s1 = mlist[0].pos; int chrom1 = mlist[0].ch; 
	int t1; int minl=6; int match=0,mism=0;

	// bigstring is unmapped portion of read, smallstring is 6bp string that is immediately next in refsequence 
	// if there is an insertion, smallstring will map with offset > 0 in bigstring 
	t1 = mlist[0].last+1; s1 = mlist[0].reflast-mlist[0].last; // these two match against each other
	for (j=t1;j<l;j++) bigstring[j-t1] = sequence[j]; bigstring[j-t1] = '\0';
	for (j=s1+t1;j<s1+t1+minl;j++) smallstring[j-s1-t1] = reflist->sequences[chrom1][j]; smallstring[j-s1-t1] = '\0';
	//fprintf(stdout,"\ncheck for insertion fl1 > 0 fl2=0 %d %s %s |",t1,smallstring,bigstring);

	if ((pch= strstr(bigstring,smallstring)) != NULL)
	{
		match = 0; mism =0; delta = (int)(pch-bigstring); t3= t1; 
		while (t3+delta  < l && s1+t3 < reflist->lengths[chrom1])
		{
			if (reflist->sequences[chrom1][s1+t3] == sequence[t3+delta]) match++;else mism++; t3++;
		}
		if ( (match >= 8 && mism <= 1) || (match >= 7 && mism ==0) )
		{
			mlist[*clusters].reffirst = s1+t1; mlist[*clusters].reflast = s1+t3-1;
			mlist[*clusters].first = t1+(int)(pch-bigstring); mlist[*clusters].last = t3+(int)(pch-bigstring)-1; 
			mlist[*clusters].el = 0; mlist[*clusters].er = 0; 
			mlist[*clusters].pos = mlist[*clusters].reffirst-mlist[*clusters].first;
			(*clusters)++;
			fprintf(stdout,"end_indel_1 %d:%d delta %d %s %s | ",match,mism,delta,smallstring,bigstring);
			free(smallstring);free(bigstring); return 1;
		}
	}

	// bigstring is reference sequence, smallstring is string 'minl' from end of read 
	for (j=l-minl-1;j<l;j++) smallstring[j-(l-minl-1)] = sequence[j]; smallstring[j-(l-minl-1)] = '\0';
	for (j=t1+s1;j<s1+t1+wts+20;j++) bigstring[j-t1-s1] = reflist->sequences[chrom1][j]; bigstring[j-t1-s1] = '\0';
	pch = strstr(bigstring,smallstring); 
	offset = l-minl-1-t1;
	if (pch != NULL && offset > 0)
	{
		delta = (int)(pch-bigstring-offset); t3 = l-1; match =0; mism =0;
		while (t3 >=t1 && mism < 2) // extending till two mismatches encountered 
		{
			if (reflist->sequences[chrom1][s1+t3+delta] == sequence[t3]) match++;   else mism++; t3--;
		}
		if (reflist->sequences[chrom1][s1+delta+t3] != sequence[t3]) { t3++; mism--; } 
		if ( (match >= 9 && mism <= 1) || (match >= 8 && mism ==0) )
		{
			mlist[*clusters].reffirst = s1+t3+delta; mlist[*clusters].reflast = s1+l-1+delta;
			mlist[*clusters].first = t3; mlist[*clusters].last = l-1; 
			mlist[*clusters].el = 0; mlist[*clusters].er = 0; 
			mlist[*clusters].pos = mlist[*clusters].reffirst-mlist[*clusters].first;
			(*clusters)++;
			fprintf(stdout,"end_indel_2 %d:%d delta %d %s %s ",match,mism,delta,smallstring,bigstring);
			//fprintf(stdout,"end deletion-L %d:%d delta %d %s %s",match,mism,delta,smallstring,bigstring);
			free(smallstring);free(bigstring); return 1;
		}
	}
	free(smallstring);free(bigstring); return 0;
}

// a read whose second part mapped (40S60M cigar) : check for deletion/insertion at start of read 
int read_secondpartmapped(char* sequence,int l,REFLIST* reflist,MATCH mlist[],int* clusters)
{
	char* bigstring = (char*)malloc(1024); char* smallstring = (char*)malloc(1024); char* pch;
	int j=0,t3=0,t4=0,delta=0,wts=40; int chrom = mlist[0].ch;
	int mapped=0;   int minl=6; int match=0,mism=0;
	int t2 = mlist[0].first-1; //first position on read before match 
	int s2 = mlist[0].reffirst-mlist[0].first; // first position on reference that is matched to t2

	for (j=0;j<=t2;j++) bigstring[j] = sequence[j]; bigstring[j] = '\0';
        for (j=s2+t2+1-minl;j<=s2+t2;j++) smallstring[j-(s2+t2+1-minl)] = reflist->sequences[chrom][j];
        smallstring[j-(s2+t2+1-minl)] = '\0';

        if ((pch = strstr(bigstring,smallstring)) != NULL)
        {
                //fprintf(stdout,"check start t2:%d s2:%d %s %s ",t2,s2,smallstring,bigstring);
                match = 0; mism =0; delta = (int)(pch-bigstring); t3= t2; t4 = delta;
                // extend the match to the left 
                while (t4+minl-1 >= 0 )
                {
                        if (reflist->sequences[chrom][s2+t3] == sequence[t4+minl-1]) match++; else mism++;
                        t3--; t4--;
                }
                if ( (match >= 8 && mism <= 1) || (match >= 7 && mism ==0) )
                {
                        fprintf(stdout,"start_indel_1 %d:%d %s:%s ",match,mism,smallstring,bigstring);
                        mlist[*clusters].reffirst = s2+t3+1; mlist[*clusters].reflast = s2+t2;
                        mlist[*clusters].first = t4+minl; mlist[*clusters].last = delta+minl-1;
                        mlist[*clusters].el = 0; mlist[*clusters].er = 0;
                        mlist[*clusters].pos = mlist[*clusters].reffirst-mlist[*clusters].first;
                        (*clusters)++;
                        mapped =1; free(smallstring);free(bigstring); return mapped;
                }
        }

	// smallstring is reference sequence 0-6 bp, bigstring is reference sequence uptil the clipped
	for (j=s2+t2-wts;j<s2+t2+1;j++) bigstring[j-(s2+t2-wts)] = reflist->sequences[chrom][j]; bigstring[j-(s2+t2-wts)] = '\0';
	for (j=0;j<minl;j++) smallstring[j] = sequence[j]; smallstring[j] = '\0';
	pch = strstr(bigstring,smallstring);
	if (pch != NULL)
	{
		match = 0; t3 = t2; mism =0; delta = (int)(pch-bigstring-wts+t2); t3 = minl; match = minl;
		while (t3 <=t2 && mism < 2) 
		{              
			if (reflist->sequences[chrom][s2+delta+t3] == sequence[t3]) match++; else mism++;
			t3++; 
		}
		if (reflist->sequences[chrom][s2+delta+t3] != sequence[t3]) { t3--; mism--; } 
		if ( (match >= 9 && mism <= 1) || (match >= 8 && mism ==0) )
		{
			mlist[*clusters].reffirst = s2+delta; mlist[*clusters].reflast = s2+delta+t3;
			mlist[*clusters].first = 0; mlist[*clusters].last = t3; 
			mlist[*clusters].el = 0; mlist[*clusters].er = 0; 
			mlist[*clusters].pos = mlist[*clusters].reffirst-mlist[*clusters].first;
			(*clusters)++;
			fprintf(stdout,"start_indel_2 %d:%d delta %d %s:%s ",match,mism,delta,smallstring,bigstring);
			mapped =1; free(smallstring);free(bigstring); return mapped;
		}
	}
	free(smallstring);free(bigstring); return mapped;
}

/*
	// special case when first hit that has same diagonal as previous hit...
	if ((*nmatches) > 0 && mlist[(*nmatches)-1].pos == (kmertable->locarray[first][1]-shift) && (mlist[(*nmatches)-1].last-mlist[(*nmatches)-1].first) == kmertable->kmer && kmertable->gap ==0) 
	{
		i = shift-mlist[(*nmatches)-1].first;
		if (i > 0 && i <= kmertable->kmer) 
		{ 
			mlist[(*nmatches)-1].last = shift + kmertable->kmer-1; mlist[(*nmatches)-1].reflast = kmertable->locarray[first][1] + kmertable->kmer-1;
			first = kmertable->locarray[first][2]; 
			fprintf(stdout,"fast hit \n");
		} 
		else if (i < 0 && i >= -1*kmertable->kmer) 
		{ 
			mlist[(*nmatches)-1].first = shift; mlist[(*nmatches)-1].reffirst = kmertable->locarray[first][1];
			first = kmertable->locarray[first][2]; 
			fprintf(stdout,"fast hit \n");
		} 
	}
*/
