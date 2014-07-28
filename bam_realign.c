/* functions for comparing an aligned sequence read to the set of variants to identify alleles and haplotype-informative reads */
//CL is cigar offset that is being evaluated 4S 60M 10I 5M  then CL = 1 for 60M, cigarlist is copy of bam->cigar
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//int MIN_BQ = 10;
//char INT_CIGAROP[] = {'M','I','D','N','S','H','P','E','X'};
#include "realignment.c"

int bamread_realignment(struct alignedread* read,VARIANT* varlist,REFLIST* reflist,struct REALIGNMENT realignments[],int* ovlist,int ov);

// identify list of indels that overlap the read 
int extract_variants_read(struct alignedread* read,CHROMVARS* chromvars,VARIANT* varlist,int chrom,REFLIST* reflist,struct REALIGNMENT realignments[])
{
	int start = read->position; int end = start + read->span; int ss=0,firstvar=0,j=0,ov=0,i=0; // int lastvar=0;
	int l1=0,l2=0; // l1 is advance on read, l2 is advance on reference genome
	int op=0,ol=0,delta=0,indelmatch=0;

	// better strategy is to search for XC flag from BWA and use it to determine if we should ignore clipped portion
	// only extend into unmapped region if it is of good quality using base quality values 
	// assumption that soft clip is either the first or last cigar, clip 'H' from cigarlist for this to be always true 
	op = read->cigarlist[0]&0xf;  ol = read->cigarlist[0]>>4; if (op == BAM_CSOFT_CLIP) start -= ol; 
	op = read->cigarlist[read->cigs-1]&0xf;  ol = read->cigarlist[read->cigs-1]>>4; if (op == BAM_CSOFT_CLIP) end += ol;
	// only extend into unmapped region if it is of good quality 
	
	j = (int)(start/BSIZE); if (j >= chromvars[chrom].blocks) return 0; // another BUG april29 2011 found here 
	ss = chromvars[chrom].intervalmap[j];  if (ss < 0 || ss >= VARIANTS) return 0;
	// check if ss is less than first variant for the chromosome 'chrom', if so assign it to the first variant 
	if (ss < chromvars[chrom].first) ss = chromvars[chrom].first; 

	//static int ovlist[1024]; 
	int* ovlist = calloc(100,sizeof(int)); 
	if (varlist[ss].position <= end)
	{
		while (varlist[ss].position+delta > start && ss >0 ) ss--; 
		if (varlist[ss].type < 0) delta = -1*varlist[ss].type; else delta = 0;
		while(varlist[ss].position+delta < start && ss < VARIANTS && ss <= chromvars[chrom].last) ss++; firstvar = ss;
		//fprintf(stdout,"variant %s:%d:%s:%s ",varlist[ss].chrom,varlist[ss].position,varlist[ss].RA,varlist[ss].AA);
		while (varlist[ss].position <= end && ss < VARIANTS && ss <= chromvars[chrom].last && ov < 100)
		{
			ovlist[ov] = ss; 
			if (varlist[ss].rightshift < 0) calculate_rightshift(varlist,ss,reflist,read->tid);  // calculate rightshift
			ov++; ss++;
		}
		//lastvar = ss;
	}
	//fprintf(stderr,"chrom %d variants %d ss %d first %d-%d span %d-%d\n",chrom,VARIANTS,ss,chromvars[chrom].first,chromvars[chrom].last,start,end);
	if (ov < 1)  { 	free(ovlist); return 0; }
	ss = firstvar; // use variable firstvar to store first variant that overlaps this read
	//fprintf(stdout,"%s %s %s %s:%d ",read->sequence,read->quality,read->readid,read->chrom,read->position);
	//for (i=0;i<read->fcigs;++i) fprintf(stdout,"%d%c:",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]);
	//fprintf(stdout,"ss %d %d %d %d-%d ov %d\n",firstvar,varlist[ss].position,read->position,start,end,ov);

	j=0;
	for (i=0;i<read->cigs;i++)
	{
		//fprintf(stdout,"%c %d \t",(char)read->cigarlist[i+1],read->cigarlist[i]); 
		if (varlist[ss].position > end) break; 
		op = read->cigarlist[i]&0xf; ol = read->cigarlist[i]>>4;
		if (op == BAM_CMATCH)
		{
			if (varlist[ss].type < 0) delta = -1*varlist[ss].type; else delta= 0;
			while (varlist[ss].position+delta >= read->position+l2 && varlist[ss].position+delta < read->position +l2 + ol)
                        {
                                ss++; //ovlist[ov].var = ss; ovlist[ov].cig = i; ovlist[ov].position = l2; ss++; ov++;  
                        }
                        while (varlist[ss].position >= read->position + l2 && varlist[ss].position < read->position + l2 + ol && varlist[ss].position + delta >= end)
                        {
                                ss++; //ovlist[ov].var = ss; ovlist[ov].cig = i; ovlist[ov].position = l2; ss++; ov++;  
                        }
			l1 += ol; l2 +=ol;
		}
		else if (op == BAM_CINS) 
		{ 
			if (varlist[ss].position == read->position+l2 && varlist[ss].type == ol) 
			{
				//fprintf(stdout,"insertion in read %d %s %d %dI %d %d\n",ss,varlist[ss].chrom,varlist[ss].position,ol,read->position+l2,varlist[ss].type);
				while (ovlist[j] < ss) j++; if (ovlist[j] ==ss) ovlist[j] = -1; 
				indelmatch = 1; 
				ss++; 
			}
			else if (varlist[ss].position >= read->position+l2 && varlist[ss].position < read->position +l2 +ol) ss++;
			l1 += ol;
		}
		else if (op == BAM_CDEL) 
		{ 
			if (varlist[ss].position == read->position+l2 && -1*varlist[ss].type == ol) 
			{
				//fprintf(stdout,"deletion in read %d %s %d %dD \n",ss,varlist[ss].chrom,varlist[ss].position,ol);
				while (ovlist[j] < ss) j++; if (ovlist[j] ==ss) ovlist[j] = -1; 
				indelmatch = 1; 
				ss++;
			}
			else if (varlist[ss].position >= read->position+l2 && varlist[ss].position < read->position +l2 + ol) ss++; 
			l2 += ol;
		}
		else if (op == BAM_CREF_SKIP) l2 += ol; 
		else if (op == BAM_CSOFT_CLIP && i ==0)
		{
			if (varlist[ss].type < 0) delta = varlist[ss].length; else delta = 0;
                        while (varlist[ss].position+delta >= read->position+l2-ol && varlist[ss].position < read->position +l2) ss++;
			l1 += ol; 
		}
		else if (op == BAM_CSOFT_CLIP) 
		{
			while (varlist[ss].position >= read->position+l2 && varlist[ss].position < read->position +l2+ol) ss++;
			l1 += ol; 
		}
		else if (op == BAM_CHARD_CLIP) l2 += 0; 
	}
	if (indelmatch ==1 && ov ==1) 
	{
		// if indelmatch ==1, change the ovlist[i] = -1
		free(ovlist); return 0; 
	}
	int best = bamread_realignment(read,varlist,reflist,realignments,ovlist,ov);
	free(ovlist); //free(ncigarlist); 
	return best;
}

// for each indel variant that overlaps the read (allowing for clipped portion as well), evaluate if modifying the alignment using this indel improves the alignment (extends into clipped or reduces mismatches at ends of read)
int bamread_realignment(struct alignedread* read,VARIANT* varlist,REFLIST* reflist,struct REALIGNMENT realignments[],int* ovlist,int ov)
{
	int i=0,j=0; 
	fprintf(stdout,"%s %s %s %s:%d ",read->sequence,read->quality,read->readid,read->chrom,read->position);
	for (i=0;i<read->fcigs;++i) fprintf(stdout,"%d%c:",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]);
	fprintf(stdout," XM:%d XC:%d\n",read->mismatches,read->XC);

	int alignments =0; int rlflag =0; int best = -1; int bestscore[2] = {-100,100};  int ss=0;
	//uint32_t* ncigarlist = calloc(100,sizeof(uint32_t)); uint32_t ncigs =0;
	for (i=0;i<ov;i++)
	{
		if (ovlist[i] < 0) continue;  ss = ovlist[i];
		fprintf(stdout,"variant %s:%d:%s:%s AB:%d ",varlist[ss].chrom,varlist[ss].position,varlist[ss].RA,varlist[ss].AA,varlist[ss].rightshift);
		print_indelhaplotypes(varlist,ss,reflist,read->tid);
		realignments[0].cigs = 0; rlflag = 0;
		rlflag = realignread_LR(read,reflist,varlist,ss,&realignments[0]); 
		if (rlflag ==1)
		{
			realignments[alignments+1].cigs = realignments[0].cigs; 
			for (j=0;j<realignments[0].cigs;j++) realignments[alignments+1].cigarlist[j] = realignments[0].cigarlist[j]; 
			realignments[alignments+1].added = realignments[0].added;  
			realignments[alignments+1].newpos = realignments[0].newpos;  
			realignments[alignments+1].delta = realignments[0].delta; realignments[alignments+1].mismatches = realignments[0].mismatches; 
			alignments++;
		}

		fprintf(stdout," === ");
		realignments[0].cigs = 0; rlflag = 0;
		rlflag = realignread_RL(read,reflist,varlist,ss,&realignments[0]); 
		if (rlflag ==1)
		{
			realignments[alignments+1].cigs = realignments[0].cigs; 
			for (j=0;j<realignments[0].cigs;j++) realignments[alignments+1].cigarlist[j] = realignments[0].cigarlist[j]; 
			realignments[alignments+1].added = realignments[0].added; 
			realignments[alignments+1].newpos = realignments[0].newpos;  
			realignments[alignments+1].delta = realignments[0].delta; realignments[alignments+1].mismatches = realignments[0].mismatches; 

			// check if cigars are identical or not...to avoid duplicate LR and RL alignments
			if (realignments[0].cigs != realignments[alignments].cigs || alignments ==0) alignments++; 
			else
			{
				for (j=0;j<realignments[0].cigs;j++) 
				{
					if (realignments[alignments].cigarlist[j] != realignments[0].cigarlist[j]) { alignments++; break; } 
				}
			}
		}
		fprintf(stdout,"\n");
		// realignread function call
	}
	if (alignments ==0) 
	{
		fprintf(stdout,"\n");
		return -1;
	}

	// ERR024173.17348618
	// find best realignment if there are multiple ones, if equally good return -1, else return index of best realignment...
	// having two identical alignments (RL and LR) is possible -> treat as one...
	for (i=0;i<alignments;i++) 
	{
		if (realignments[i+1].delta > bestscore[0] && realignments[i+1].mismatches < bestscore[1]) 
		{
			bestscore[0] = realignments[i+1].delta; bestscore[1] = realignments[i+1].mismatches; best = i+1;
		}
		else if (realignments[i+1].delta == bestscore[0]) 
		{ 
			//fprintf(stdout," %d:%d:%d:%d ",realignments[i+1].delta,realignments[i+1].mismatches,bestscore[0],i); 
			best = -1; break; 
		} 
	}
	//if (alignments > 1 && best >= 0 && realignments[best].delta < 2) best = -1; // for multiple realignment case, require delta >=2

	if (best >= 0) fprintf(stdout,"best %d %d %d \t",best,realignments[best].delta,realignments[best].mismatches);
	else fprintf(stdout,"best -1 ");
	fprintf(stdout,"realignments %d\n\n",alignments);

	return best; 
}


