
int generate_cigar_alignment(struct alignedread* read,MATCH* mlist,int clusters,REFLIST* reflist,int bestnode);

// pairwise clumping of hits, need to check if three clusters need to be merged merge hits that are on the same diagonal or small distance prior to this 
// the hits are sorted by query position, DAG -> heaviest path in the graph -> find heaviest path ending at each node/cluster
int clump_hits(struct alignedread* read, MATCH* mlist,int clusters,REFLIST* reflist)
{
	if (clusters >= 30) return -1;
	float MT = 1, MP = 3, GO = 4.5, GE = 1.5; // log-scale penalty
	int i=0,j=0,k=0; int delta = 0; int joinflag = 1; int b1 =0,b2=0,b12=0; // no of bases of read covered by a pair of hits
	float score = 0; int bnode=0; int overlap =0;
	int Dlength = 0, Ilength =0; 
	int MIN_SCORE = 30;
	int clipped_j=0,clipped_bnode=0;

	// add clipping penalty for hits that don't cover entire read...
	// extend the first or last hit in list to cover the beginning or end of read
	for (i=0;i<clusters;i++) 
	{ 
		mlist[i].bestscore = (mlist[i].last-mlist[i].first+1)*MT; 
		if (mlist[i].first > 0) mlist[i].clippenalty = 0;
		mlist[i].previous = -1; 
	} 
	//  the only case not covered is when there is a substitution close to insertion/deletion -> treat as DI or ??

	// complexity is O(n^2) where n is no of clusters, could be reduced by not considering all nodes prior to it but only a subset
	for (j=1;j< clusters;j++)
	{
		for (i=j-1;i>=0;i--) // for all nodes prior to it starting from the previous node
		{
			delta = mlist[i].pos - mlist[j].pos;
			// allow 10 bp extension of insertion
			if (mlist[i].reflast >= mlist[j].reflast+15 || mlist[i].last > mlist[j].last) continue; // impossible to join hits
			b1 = mlist[i].last-mlist[i].first+1; b2 = mlist[j].last-mlist[j].first+1; b12 = b1; 
			if (mlist[j].first > mlist[i].last) b12 += b2; 
			else if (mlist[j].last > mlist[i].last) b12 += mlist[j].last-mlist[i].last; 
			overlap = mlist[i].last - mlist[j].first+1;  // if overlap is < 0 => no overlap
			//if (b12 < b1 + 4 || b12 < b2 + 4) joinflag = 0; // filter removed 07/24/13 check this
			//fprintf(stdout,"%d:%d flag %d \n",i,j,joinflag); if (joinflag ==0) continue; 

			score = 0;
			// if delta ==0 -> consider mismatch score as alternate
			if (delta ==0)
			{
				score = mlist[i].bestscore + (mlist[j].last-mlist[j].first+1)*MT;
				if (overlap > 0) score -= MT*overlap; 
				k=0;
				while (mlist[i].last+k+1 < mlist[j].first) 
				{
					if (mlist[i].reflast+1+k < reflist->lengths[mlist[i].ch] && read->sequence[mlist[i].last+1+k] == reflist->sequences[mlist[i].ch][mlist[i].reflast+1+k]) score += MT;
					else score -= MP; 
					k++;
				}
			}
			else if (delta < 0)  // deletion of bases 
			{
				if (overlap >= 0)  // deletion with overlap >= 0 as expected 
				{
					Dlength = -1*delta; Ilength =0;
					score = mlist[i].bestscore -MT*overlap + (mlist[j].last-mlist[j].first+1)*MT;
					score -= GO; score -= log(-1*delta)*GE; 
					if (score >= MIN_SCORE) fprintf(stdout,"D-clump:%d:%d:%d:%d:score:%0.1f:%dD ",i,j,delta,b12,score,Dlength);
				}
				else // deletion with some insertion or mismatches  D > I
				{
					Dlength = -1*(overlap+delta); Ilength = -1*overlap;
					score = mlist[i].bestscore + (mlist[j].last-mlist[j].first+1)*MT; 
					score -= GO*2; score -= log(Ilength)*GE; score -= log(Dlength)*GE;
					if (score >= MIN_SCORE) fprintf(stdout,"DI-clump:%d:%d:%d:%d:score:%0.1f:%dD_%dI ",i,j,delta,b12,score,Dlength,Ilength);
				}
			}
			else if (delta > 0)  // deletion of bases
			{
				overlap *= -1; 
				if (overlap <= delta) // size of insertion is less than the gap between hits -> simple event 
				{
					Ilength = delta; Dlength =0;
					score = mlist[i].bestscore + (mlist[j].last-mlist[j].first+1)*MT - (delta-overlap)*MT; 
					score -= GO; score -= log(Ilength)*GE;
					if (score >= MIN_SCORE) fprintf(stdout,"I-clump:%d:%d:%d:%d:score:%0.1f:%dI ",i,j,delta,b12,score,Ilength);
				}
				else
				{
					Ilength = overlap; Dlength = overlap-delta;
					score = mlist[i].bestscore + (mlist[j].last-mlist[j].first+1)*MT; 
					score -= GO*2; score -= log(Ilength)*GE; score -= log(Dlength)*GE;
					if (score >= MIN_SCORE) fprintf(stdout,"ID-clump:%d:%d:%d:%d:score:%0.1f:%dD_%dI ",i,j,delta,b12,score,Dlength,Ilength);
					
				}
			}
			if (score > mlist[j].bestscore && delta ==0) 
			{ 
				mlist[j].bestscore = score; mlist[j].previous = i; mlist[j].DL =0; mlist[j].IL = 0;
				//mlist[i].last = mlist[j].first-1; mlist[i].reflast = mlist[j].reffirst-1;
				fprintf(stdout,"NOINDEL:%d:%d:%d:%d:score:%0.1f:%d ",i,j,delta,b12,score,k);
			} 
			else if (score > mlist[j].bestscore && delta != 0) 
			{ 
				mlist[j].bestscore = score; mlist[j].previous = i; mlist[j].DL = Dlength; mlist[j].IL = Ilength; 
				// length of insertions and deleted bases from previous hit
			} 
			// store the delta in score plus the cigar implied by the joining of the two fragments 
		}
	}

	// bestscore of alignment ending at cluster 'j'
	bnode = 0; 
	for (j=1;j< clusters;j++) 
	{
		// consider clipping penalty for end of read 
		clipped_bnode = read->readlength-1-mlist[bnode].last; clipped_j = read->readlength-1-mlist[bnode].last;
			
		if (mlist[j].bestscore > mlist[bnode].bestscore) bnode = j; 
		// find best path in graph ending at cluster 'j'
	}
	if (mlist[bnode].bestscore < MIN_SCORE) return 0;
	return generate_cigar_alignment(read,mlist,clusters,reflist,bnode);
}

// char INT_CIGAROP[] = {'M','I','D','N','S','H','P','E','X'};
// generate cigar of best alignment using backtrace in node list // bnode is last node in best path in dynamic programming
// return value is 1 if new alignment replaces original one
int generate_cigar_alignment(struct alignedread* read,MATCH* mlist,int clusters,REFLIST* reflist,int bnode) 
{
	// generate cigar for best path... indels are pushed leftward 
        static uint32_t cigarlist[1024]; int cigs=0; int l = read->readlength;

	int node = bnode; int gapped=0; int i=0,mm_tract=0, overlap=0;
	fprintf(stdout," | bestpath %.1f %d",mlist[node].bestscore,node);
	while (mlist[node].previous >= 0) 
	{ 
		if (mlist[node].DL > 0 || mlist[node].IL > 0) gapped++;
		fprintf(stdout,"->%d",mlist[node].previous); 
		node =mlist[node].previous; 
	} 
	node = bnode; //fprintf(stdout," | gapped %d | ",gapped);
	
	if (mlist[node].last < l-1) cigarlist[cigs++] = ((l-mlist[node].last-1)<<4) + 4;
	while (mlist[node].previous >= 0) 
	{
		if ((mlist[node].last-mlist[node].first+1-overlap) > 0) cigarlist[cigs++] = (mlist[node].last-mlist[node].first+1-overlap)<<4; 
		if (mlist[node].DL > 0 && mlist[node].IL > 0) 
		{ 
			cigarlist[cigs++] = ((mlist[node].IL)<<4) + 1;  cigarlist[cigs++] = ((mlist[node].DL)<<4) + 2;  
			overlap =0;
		} 
		else if (mlist[node].DL > 0) 
		{
			cigarlist[cigs++] = ((mlist[node].DL)<<4) + 2;
			overlap = mlist[mlist[node].previous].last-mlist[node].first+1; 
			if (overlap < 0) overlap =0; 
		}
		else if (mlist[node].IL > 0) 
		{
			cigarlist[cigs++] = ((mlist[node].IL)<<4) + 1;
			overlap = mlist[mlist[node].previous].last-(mlist[node].first-mlist[node].IL)+1; 
			if (overlap < 0) overlap =0;
		}
		else 
		{
			// same diagonal so some mismatches 
			mm_tract = mlist[node].first - mlist[mlist[node].previous].last-1;
			//if (mm_tract > 0) fprintf(stdout,"mm stretch %d \n",mm_tract);
			if (mm_tract > 0)  cigarlist[cigs++] = (mm_tract<<4)+8; 
			overlap =0;
		}
		node =mlist[node].previous; 
	} 
	//fprintf(stdout,"%dM",mlist[node].last-mlist[node].first+1);
	cigarlist[cigs++] = (mlist[node].last-mlist[node].first+1-overlap)<<4;
	if (mlist[node].first > 0) cigarlist[cigs++] = ((mlist[node].first)<<4) + 4;
	fprintf(stdout,"\n%s %d ",read->chrom,mlist[node].reffirst+1);
	for (i=cigs-1;i>=0;i--) fprintf(stdout,"%d%c:",cigarlist[i]>>4,INT_CIGAROP[cigarlist[i]&0xf]); fprintf(stdout,"\t");
	for (i=0;i<read->fcigs;i++) fprintf(stdout,"%d%c:",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]); fprintf(stdout,"\t");

	int clipped_new =0; 
	if ((cigarlist[0]&0xf) == BAM_CSOFT_CLIP) clipped_new += cigarlist[0]>>4; 
	if ((cigarlist[cigs-1]&0xf) == BAM_CSOFT_CLIP) clipped_new += cigarlist[cigs-1]>>4; 

	if (gapped > 0 && read->gaps ==0) 
	{
		fprintf(stdout," GAPALIGN-NEW clip:%d | ",clipped_new); 
		// reallocate read cigar and copy cigarlist to it 
		if (gapped <= 2 && clipped_new <= 10 && clipped_new < read->clipped)
		{ 
			if (read->cigs < cigs) read->cigarlist = (uint32_t*)realloc(read->cigarlist,sizeof(uint32_t)*cigs);
			for (i=cigs-1;i>=0;i--) read->cigarlist[cigs-1-i] = cigarlist[i]; read->cigs = cigs; read->position = mlist[node].reffirst; // 0 offset not 1
			read->gaps = gapped;
			return 1;
		}
		else fprintf(stdout," unused ");
	}
	else if (gapped > 0 && read->gaps > 0)  // trust new mapping ?? useful for complex indels
	{
		int cigar_match =1;
		for (i=cigs-1;i>=0;i--) 
		{
			if (cigs-1-i >= read->cigs || read->cigarlist[cigs-1-i] != cigarlist[i]) cigar_match = 0;
		}
		fprintf(stdout,"GAPALIGN-match clip:%d CM:%d | ",clipped_new,cigar_match);
		// reallocate read cigar and copy cigarlist to it
		if (gapped <= read->gaps && gapped <= 2 && clipped_new <= read->clipped && cigs < read->cigs) 
		{ 
			if (read->cigs < cigs) read->cigarlist = (uint32_t*)realloc(read->cigarlist,sizeof(uint32_t)*cigs);
			for (i=cigs-1;i>=0;i--) read->cigarlist[cigs-1-i] = cigarlist[i]; read->cigs = cigs; read->position = mlist[node].reffirst;	
			read->gaps = gapped;
			return 1;
		}
		else if (cigar_match ==0) fprintf(stdout," unused ");
	}
	else if (read->gaps > 0) fprintf(stdout,"GAPALIGN-fail | ");
	return 0;
}

