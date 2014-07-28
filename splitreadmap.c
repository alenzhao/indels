

//for (i=0;i<nmatches;i++) { if (mlist[i].valid != '0') fprintf(stdout,"%d:%d:%d:%d:%d  ",i,mlist[i].pos,mlist[i].pos+mlist[i].offset,mlist[i].first,mlist[i].last);}
// clusters initial hits and returns number of non-overlapping merged hits
// use next pointer in mlist to reduce iterating over all matches
int cluster_hits(char* sequence,int l,MATCH* mlist,int nmatches,REFLIST* reflist)
{
	int i=0,j=0,clusters=0; 
	if (nmatches > 1) qsort(mlist,nmatches,sizeof(MATCH),compare_hits_diagonalpos); 
	clusters =nmatches; 
	while (i < nmatches-1) // merge neighboring matches that are on same diagonal and define adjacent hits
	{
		j = i+1;  // pos is diagonal coordinate, first is read start pos, hits sorted by pos and then by first 
		while (mlist[i].pos == mlist[j].pos && (mlist[j].first - mlist[i].last) <= 1 && j < nmatches) 
		{
			if (mlist[j].last > mlist[i].last) { mlist[i].last = mlist[j].last; mlist[i].reflast = mlist[j].reflast; } 
			mlist[j].valid = '0';
			j++; clusters--;
		}
		i= j; 
	}
	// extend each unique match left and right
	for (i=0;i<nmatches;i++)
	{
		if (mlist[i].valid =='0') continue;
		while (mlist[i].last+1 < l && (sequence[mlist[i].last+1] == reflist->sequences[mlist[i].ch][mlist[i].reflast+1]))
		{
			mlist[i].last++; mlist[i].reflast++;
		}
		while (mlist[i].first-1 >=0 && sequence[mlist[i].first-1] == reflist->sequences[mlist[i].ch][mlist[i].reffirst-1])
		{
			mlist[i].first--; mlist[i].reffirst--; 
		}
	}

	// this step is unnecessary 
	for (i=1;i<nmatches;i++)
	{
		if (mlist[i].valid =='0') continue; 
		j= i+1; while (j < nmatches && mlist[j].valid == '0') j++;
		if (j < nmatches && mlist[i].pos == mlist[j].pos && mlist[i].first== mlist[j].first && mlist[i].last == mlist[j].last) 
		{ 
			mlist[i].valid = '0'; clusters--; 
		}       
	}

	// remove non valid hits from list to make smaller list
	j=0; for (i=0;i<nmatches;i++) 
	{
		if (mlist[i].valid != '0' && j <i) 
		{
			mlist[j].pos = mlist[i].pos; mlist[j].ch = mlist[i].ch; mlist[j].offset = mlist[i].offset; 
			mlist[j].first= mlist[i].first; mlist[j].last = mlist[i].last; 
			mlist[j].reffirst= mlist[i].reffirst; mlist[j].reflast = mlist[i].reflast; 
		}
		if (mlist[i].valid != '0') j++; 
	}
	nmatches = j;
	return nmatches;
}

void print_hits(struct alignedread* read,REFLIST* reflist,MATCH* mlist,int nmatches,int clusters)
{
	int  i=0; char matestrand = '+'; if ((read->flag & 32) == 32) matestrand = '-';  // set the strand = strand of mate pair

	if (nmatches < 2 && read->gaps > 0)
	{
		 fprintf(stdout,"read %s %s \t matepair %d:%d strand %c rl:%d | hits %d:%d | ",read->sequence,read->readid,read->mtid,read->mateposition,matestrand,read->readlength,clusters,nmatches);
		if (read->cigs ==0) fprintf(stdout," UNMAPPED ");
		else
		{
			for (i=0;i<read->fcigs;i++) fprintf(stdout,"%d%c:",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]);
			fprintf(stdout," %d | ",read->mismatches);
		}
		fprintf(stdout,"NOIT %d IS %d\n",nmatches,read->IS);

	}
	else if (nmatches >=2 && clusters ==1) 
	{
		// determine if start/end of read has no hits (find coordinates start:end) and do mapping with smaller k-mer 
		if (mlist[0].first < 4 && mlist[0].last-mlist[0].first >= 30)
		{
			//fprintf(stdout,"indel-near-end ");  read_firstpartmapped(read->sequence,read->readlength,reflist,mlist,&clusters);
		}
		else if (mlist[0].last -mlist[0].first >= 30 && read->readlength-mlist[0].last < 4)
		{
			//fprintf(stdout,"indel-near-start "); read_secondpartmapped(read->sequence,read->readlength,reflist,mlist,&clusters);
		}
	}
	else if (nmatches >=2 && clusters >=2) 
	{
		fprintf(stdout,"read %s %s \t matepair %d:%d strand %c rl:%d | hits %d:%d | ",read->sequence,read->readid,read->mtid,read->mateposition,matestrand,read->readlength,clusters,nmatches);
		if (read->cigs ==0) fprintf(stdout," UNMAPPED "); 
		else 
		{
			for (i=0;i<read->fcigs;i++) fprintf(stdout,"%d%c:",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]); 
			fprintf(stdout," %d | ",read->mismatches);
		}
		qsort(mlist,clusters,sizeof(MATCH),compare_hits_querypos);  // sort by position on query -> DAG 
		for (i=0;i<clusters;i++) fprintf(stdout,"%d:%d:%d-%d:%d-%d:%0.1f:%d ",i,mlist[i].pos,mlist[i].reffirst,mlist[i].reflast,mlist[i].first,mlist[i].last,mlist[i].bestscore,mlist[i].previous);
		if (read->gaps >=2) fprintf(stdout,"MGAP "); //fprintf(stdout,"\n");
	}
	//fprintf(stdout,"NM<2 read %s %s \t matepair %d:%d strand %c rl:%d | hits %d:%d %s\n",read->sequence,read->readid,read->matech,read->mateposition,matestrand,read->readlength,clusters,nmatches,read->cigar);
}

// if mate is mapped on + strand, then second read should be mapped on - strand for paired-end Illumina reads 
//  5'-------------------read1_s--------read1-e---------------------------------------------------------------------------->3'
//  3'------------------------------------------------------------read2_e-----------------read2-s--------------------------<5'

// bug in this function when readlength (due to hard clipping) is even smaller than kmer !! 10/30/2013 
// index is build only for positive strand, so if read is supposed to mapped on negative strand, we reverse complement it 
int splitreadmapping(struct alignedread* read,KMERTABLE* kmertable,KMERTABLE* kmertable1,REFLIST* reflist,unsigned char* shortseq,MATCH* mlist)
{
	if (read->readlength < 20) return 0; // minimum readlength of 20 required
	int matepos =read->mateposition;

	if (read->mateposition >= reflist->lengths[read->tid]) fprintf(stderr,"error %d %d %d \n",read->tid,read->mateposition,reflist->lengths[read->tid]);
	// seg fault due to reads that map close to ends of chromosome, condition below fixes it temporarily
	//if (read->mateposition >= reflist->lengths[read->tid]-5000 || read->mateposition <= 5000) return 0;

	char matestrand = '+'; if ((read->flag & 32) == 32) matestrand = '-';  // set the strand = strand of mate pair
	int i; int shift=0; int h=0; unsigned int hash =0; int nmatches=0; int flag =0; int clusters=0;
	//char* kmerhit = calloc(sizeof(char),read->readlength+1); for (i=0;i<read->readlength;i++) kmerhit[i] = '0'; kmerhit[i] ='\0';
	int realigned =0;

	if (matestrand == '-')   // read should be mapped on positive strand...
	{
		// read could be mapped or unmapped, if unmapped, read->flag & 16 should equal 0  
		if ((read->flag & 16) == 16) RC(read->sequence,read->sequence,read->readlength);  matepos += read->readlength; 
	} 
	else // read should be mapped on negative strand, reverse complement if original sequence is on positive strand
	{ 
		if ((read->flag & 16) == 0) RC(read->sequence,read->sequence,read->readlength);  
	}

	int last =0; int pmatches=0; int ps = 0; int s=0;
	while (flag ==0)
	{
		if (shift+kmertable->kmer > read->readlength) 
		{ 
			shift = read->readlength-kmertable->kmer; flag = 1; 
			if (shift <0) break; // added 10/30/13 to avoid bugs for very short reads
		} 
		for (i=0;i<kmertable->kmer;i++) shortseq[i] = read->sequence[i+shift]; shortseq[i] = '\0';
		hash = BTI[shortseq[0]]; for (i=1;i<kmertable->kmer;i++) { hash <<=2; hash += BTI[shortseq[i]]; } 
		h = find_hits_forseed(read->sequence,read->readlength,reflist,hash,kmertable,mlist,&nmatches,matestrand,matepos,shift,MAX_IS+MAXDEL);
		if (nmatches-pmatches > 0) 
		{
			//for (i=last;i<shift+kmertable->kmer;i++) kmerhit[i] = '1'; last = shift + kmertable->kmer;
		}
		else last = shift; 
		if (h ==0 && shift <= read->readlength-kmertable1->kmer-kmertable->gap)
		{
			for (i=0;i<kmertable1->kmer/2;i++) shortseq[i] = read->sequence[i+shift];
			for (i=kmertable1->kmer/2;i<kmertable1->kmer;i++) shortseq[i] = read->sequence[i+shift+kmertable1->gap]; shortseq[i] = '\0';
			hash = BTI[shortseq[0]]; 
			for (i=1;i<kmertable1->kmer;i++) { hash <<=2; hash += BTI[shortseq[i]]; } 
			h = find_hits_forseed(read->sequence,read->readlength,reflist,hash,kmertable1,mlist,&nmatches,matestrand,matepos,shift,MAX_IS+MAXDEL);
		}
		if (shift == read->readlength-kmertable->kmer) break; 
		//if (nmatches -pmatches <= 0 && ps == kmertable->kmer/2 ) s = -ps; 
		if (nmatches -pmatches <= 0 ) s = 1; // too many hits for this k-mer so small shift for next
		else s = kmertable->kmer/2; 
		shift += s;
		//fprintf(stdout,"%s:%d:%d:%d ",shortseq,shift,h,nmatches); //if (h >0) nmatches = h;
		pmatches = nmatches; ps =s;
	}
	// search the second hash table for part of the read with no covering hit..
	if (nmatches >=2)  clusters = cluster_hits(read->sequence,read->readlength,mlist,nmatches,reflist); // cluster short hits to form ungapped hits
	else clusters = nmatches;

	if (clusters >= 1) print_hits(read,reflist,mlist,nmatches,clusters);
	if (nmatches >=2 && clusters >=2) 
	{
		realigned = clump_hits(read,mlist,clusters,reflist); 
		fprintf(stdout,"RL:%d\n\n",realigned);
	}
	return realigned;
	//free(kmerhit);
}

