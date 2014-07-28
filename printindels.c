// functions for storing, clustering indels identified from parsing a bam/sam file as well as split-read alignment 
// this are used in split-read alignment code
int MIN_READS_FOR_INDEL = 3;

// for each indel cluster -> require one good aligned read to avoid aggressive indels... 
// mark duplicates among indel reads (91M2I7M, 91M2I7M -> count as one) 
// for indels identified exclusively from split-read mapping, require at least one indel on both fwd and reverse strands

// sort indel list by starting position and  cluster the indels
void cluster_indelreads(SPINDEL* indelreadlist,int indelreads,REFLIST* reflist,int partial)
{
	int i=0,j=0,cluster=0;//,t1=0,t2=0,match=0,mism=0;
	fprintf(stderr,"no of reads that support indels %d sorting list of indelreads by coordinates \n",indelreads);
	fprintf(stdout,"no of reads that support indels %d sorting list of indelreads by coordinates \n",indelreads);
	qsort(indelreadlist,indelreads,sizeof(SPINDEL),compare_indel);
	cluster =0;
	for (i=0;i<indelreads;i++)
	{
		indelreadlist[i].readsf =0; indelreadlist[i].readsr =0; indelreadlist[i].newreads =0; indelreadlist[i].reads = 0;
		if (indelreadlist[i].strand == '+') indelreadlist[i].readsf =1; else indelreadlist[i].readsr =1;
		if (indelreadlist[i].code == 'S') indelreadlist[i].newreads =1;  // reads found from split-read mapping
		else indelreadlist[i].reads = 1;
		indelreadlist[i].cluster = i; // clusterid
	}
	for (i=0;i<indelreads-1;i++)
	{
		if (indelreadlist[i].cluster !=i) continue;
		for (j=i+1;j<indelreads && indelreadlist[j].position <= indelreadlist[i].position;j++)
		{
			if (indelreadlist[i].position == indelreadlist[j].position && indelreadlist[i].length == indelreadlist[j].length && indelreadlist[i].partial != '1' && indelreadlist[j].partial != '1')
			{
				cluster = indelreadlist[i].cluster;
				indelreadlist[cluster].readsf += indelreadlist[j].readsf; indelreadlist[j].readsf = 0;
				indelreadlist[cluster].readsr += indelreadlist[j].readsr; indelreadlist[j].readsr = 0;
				indelreadlist[cluster].newreads += indelreadlist[j].newreads; indelreadlist[j].newreads = 0;
				indelreadlist[cluster].reads += indelreadlist[j].reads; indelreadlist[j].reads = 0;
				indelreadlist[j].cluster = cluster;
			}
			//if (indelreadlist[i].position < indelreadlist[j].position) break;
		}
	}
}

void print_clusters(SPINDEL* indelreadlist,int indelreads,REFLIST* reflist,int print_partial_indels,FILE* outfp)
{
	int FLANKING = 100;
	static int printvcf = 0; char code; int k; int near_interval= 1; int k1=0;
	printvcf++;
	if (printvcf ==1) fprintf(outfp,"##fileformat=VCFv4.0\n##INDELS-from-BAMfile\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

	if (TARGETED ==1) reflist->cinterval = reflist->first_interval_chrom[indelreadlist[0].chrom]; 

	int i=0,t1=0;
	for (i=0;i<indelreads;i++)
	{
		if (indelreadlist[i].cluster !=i || indelreadlist[i].readsf + indelreadlist[i].readsr < 3) continue;
		// filter on small indels, require at least three reads for indels less than length 3
		if (abs(indelreadlist[i].length) < 5 && indelreadlist[i].readsf + indelreadlist[i].readsr < MIN_READS_FOR_INDEL) continue;
		if ((indelreadlist[i].readsf + indelreadlist[i].readsr == indelreadlist[i].newreads) && (indelreadlist[i].readsf < 1 || indelreadlist[i].readsr < 1)) continue;
		if (indelreadlist[i].partial == '1') continue; 

		if (TARGETED ==1)
		{
			k = indelreadlist[i].position; near_interval = 0;
			while (reflist->intervallist[reflist->cinterval].end < k && reflist->cinterval < reflist->intervals-1) reflist->cinterval++; 
                        if (k >= reflist->intervallist[reflist->cinterval].start && k <= reflist->intervallist[reflist->cinterval].end )near_interval = 2;
			else if (reflist->intervallist[reflist->cinterval].start-k > 0 && reflist->intervallist[reflist->cinterval].start-k < FLANKING ) near_interval = 1;
			else if (reflist->cinterval-1 >=0 && k-reflist->intervallist[reflist->cinterval-1].end > 0 &&  k-reflist->intervallist[reflist->cinterval-1].end < FLANKING) near_interval = 1;
			if ( (indelreadlist[i].cigar[0]&0xf) == BAM_CDEL) 
			{
				k += indelreadlist[i].cigar[0]>>4;
				k1 = reflist->cinterval; while (reflist->intervallist[k1].end < k && k1 < reflist->intervals-1) k1++;
                        	if (k >= reflist->intervallist[k1].start && k <= reflist->intervallist[k1].end )near_interval = 2;
				else if (reflist->intervallist[k1].start-k > 0 && reflist->intervallist[k1].start-k < FLANKING ) near_interval = 1;
				else if (k1-1 >=0 && k-reflist->intervallist[k1-1].end > 0 &&  k-reflist->intervallist[k1-1].end < FLANKING) near_interval = 1;
			}
			if (near_interval ==0) continue;
			//if (reflist->cinterval >= 1) fprintf(outfp,"near interval %d %d \t",reflist->intervallist[reflist->cinterval-1].start,reflist->intervallist[reflist->cinterval-1].end);
			//fprintf(outfp,"near interval %d %d val %d\n",reflist->intervallist[reflist->cinterval].start,reflist->intervallist[reflist->cinterval].end,near_interval);
		}
 
		//fprintf(stdout,"VCFINDEL:%s\t%d\t%s\t",reflist->names[indelreadlist[i].chrom],indelreadlist[i].position,indelreadlist[i].code);
		if ( (indelreadlist[i].cigar[0]&0xf) == BAM_CDEL && indelreadlist[i].cigar[1] ==0) 
		{
			fprintf(outfp,"%s\t%d\t.\t",reflist->names[indelreadlist[i].chrom],indelreadlist[i].position);
			for(t1=indelreadlist[i].position-1;t1<indelreadlist[i].position+(indelreadlist[i].cigar[0]>>4);t1++) fprintf(outfp,"%c",reflist->sequences[indelreadlist[i].chrom][t1]);
			fprintf(outfp,"\t");
			for(t1=indelreadlist[i].position-1;t1<indelreadlist[i].position;t1++) fprintf(outfp,"%c",reflist->sequences[indelreadlist[i].chrom][t1]);
		}
		else if ( (indelreadlist[i].cigar[0]&0xf) == BAM_CINS && indelreadlist[i].cigar[1] ==0) 
		{
			fprintf(outfp,"%s\t%d\t.\t",reflist->names[indelreadlist[i].chrom],indelreadlist[i].position);
			for(t1=indelreadlist[i].position-1;t1<indelreadlist[i].position;t1++) fprintf(outfp,"%c",reflist->sequences[indelreadlist[i].chrom][t1]);
			fprintf(outfp,"\t%c",reflist->sequences[indelreadlist[i].chrom][indelreadlist[i].position-1]);
			for(t1=0;t1<(indelreadlist[i].cigar[0]>>4);t1++) fprintf(outfp,"%c",indelreadlist[i].insertedseq[t1]);
			//for (t1=indelreadlist[i].firstbase-1;t1<indelreadlist[i].lastbase;t1++) fprintf(outfp,"%c",
			
		}
		else if ( (indelreadlist[i].cigar[0]&0xf) == BAM_CDEL && (indelreadlist[i].cigar[1]&0xf) ==BAM_CINS) 
		{
			fprintf(outfp,"%s\t%d\t.\t",reflist->names[indelreadlist[i].chrom],indelreadlist[i].position);
			for(t1=indelreadlist[i].position-1;t1<indelreadlist[i].position+(indelreadlist[i].cigar[0]>>4);t1++) fprintf(outfp,"%c",reflist->sequences[indelreadlist[i].chrom][t1]);

			fprintf(outfp,"\t%c",reflist->sequences[indelreadlist[i].chrom][indelreadlist[i].position-1]);
			for(t1=0;t1<(indelreadlist[i].cigar[1]>>4);t1++) fprintf(outfp,"%c",indelreadlist[i].insertedseq[t1]);
			
		}
		//fprintf(outfp,"%s\t%s",indelreadlist[i].reference,indelreadlist[i].alternate);
		fprintf(outfp,"\t.\tPASS\tSREADS=%d,%d,%d,%d;",indelreadlist[i].reads,indelreadlist[i].readsf,indelreadlist[i].readsr,indelreadlist[i].newreads); 
		if (indelreadlist[i].reads ==0) code = 'S'; else code = 'G';
		fprintf(outfp,"HOMLEN=%d;SVLEN=%d;TYPE=%c",indelreadlist[i].ambiguity,indelreadlist[i].length,code);
		if (indelreadlist[i].ambiguity >0) fprintf(outfp,";HOMSEQ=");
		for (t1=indelreadlist[i].position;t1<indelreadlist[i].position+indelreadlist[i].ambiguity;t1++) fprintf(outfp,"%c",reflist->sequences[indelreadlist[i].chrom][t1]); 
		if (TARGETED ==1 && near_interval ==2)
		{
			fprintf(outfp,";INTERVAL=%d:%d",reflist->intervallist[reflist->cinterval].start,reflist->intervallist[reflist->cinterval].end);
		}
		fprintf(outfp,"\n");
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////// CODE THAT was used in previous version //////////

/*
// print indel haplotype consensus sequences
if (indelreadlist[i].length < 0)
{
	fprintf(stdout,">%s,%d,D,",reflist->names[indelreadlist[i].chrom],indelreadlist[i].position+1);
	if (-1*indelreadlist[i].length < 25) { for (t1=indelreadlist[i].firstbase;t1<indelreadlist[i].lastbase;t1++) fprintf(stdout,"%c",reflist->sequences[indelreadlist[i].chrom][t1]); }
	else
	{
		for (t1=0;t1<5;t1++) fprintf(stdout,"%c",reflist->sequences[indelreadlist[i].chrom][t1+indelreadlist[i].firstbase]);
		fprintf(stdout,"<%d>",-1*indelreadlist[i].length-10);
		for (t1=5;t1>0;t1--) fprintf(stdout,"%c",reflist->sequences[indelreadlist[i].chrom][indelreadlist[i].lastbase-t1]);
	}
	leftlen = MRLENGTH-indelreadlist[i].ambiguity-1; rightlen = MRLENGTH-1;
	fprintf(stdout,",%d,%d,%dM%dD%dM,%d\t",indelreadlist[i].position-leftlen+1,indelreadlist[i].position+rightlen+1,leftlen,-1*indelreadlist[i].length,rightlen,indelreadlist[i].ambiguity);

	for (t1=indelreadlist[i].position-leftlen;t1<indelreadlist[i].position;t1++) fprintf(stdout,"%c",reflist->sequences[indelreadlist[i].chrom][t1]);
	for (t1=indelreadlist[i].position-indelreadlist[i].length;t1<indelreadlist[i].position-indelreadlist[i].length+MRLENGTH;t1++) fprintf(stdout,"%c",reflist->sequences[indelreadlist[i].chrom][t1]);
	fprintf(stdout,"\n\n");
}
if (indelreadlist[i].length > 0)
{
	fprintf(stdout,">%s,%d,I,",reflist->names[indelreadlist[i].chrom],indelreadlist[i].position+1);
	for (t1=indelreadlist[i].firstbase;t1<indelreadlist[i].lastbase;t1++) fprintf(stdout,"%c",indelreadlist[i].sequence[t1]);
	leftlen = MRLENGTH-indelreadlist[i].ambiguity-1; rightlen = MRLENGTH-1;
	fprintf(stdout,",%d,%d,%dM%dI%dM,%d\t",indelreadlist[i].position-leftlen+1,indelreadlist[i].position+rightlen+1,leftlen,indelreadlist[i].length,rightlen,indelreadlist[i].ambiguity);

	for (t1=indelreadlist[i].position-leftlen;t1<indelreadlist[i].position;t1++) fprintf(stdout,"%c",reflist->sequences[indelreadlist[i].chrom][t1]);
	for (t1=indelreadlist[i].firstbase;t1<indelreadlist[i].lastbase;t1++) fprintf(stdout,"%c",indelreadlist[i].sequence[t1]);
	for (t1=indelreadlist[i].position;t1<indelreadlist[i].position-indelreadlist[i].length+MRLENGTH;t1++) fprintf(stdout,"%c",reflist->sequences[indelreadlist[i].chrom][t1]);
	fprintf(stdout,"\n\n");
}
		if (indelreadlist[i].partial == '1' && partial ==1)
		{
			fprintf(stdout,"%s %s %d %d %d %s ",indelreadlist[i].code,reflist->names[indelreadlist[i].chrom],indelreadlist[i].position+1,indelreadlist[i].L,indelreadlist[i].match,indelreadlist[i].readid);
			for (j=0;j<indelreadlist[i].L;j++) fprintf(stdout,"%c",indelreadlist[i].sequence[j]);
			fprintf(stdout,"_"); for (j=indelreadlist[i].L;j<indelreadlist[i].rl;j++) fprintf(stdout,"%c",indelreadlist[i].sequence[j]);
			fprintf(stdout,"\n");
		}
		else if (indelreadlist[i].partial != '1')
		{
			fprintf(stdout,"%s %s %d %d %d ",indelreadlist[i].code,reflist->names[indelreadlist[i].chrom],indelreadlist[i].position+1,indelreadlist[i].position+indelreadlist[i].ambiguity+1,indelreadlist[i].length);
			if (indelreadlist[i].length < 0) { for (t1=indelreadlist[i].firstbase;t1<indelreadlist[i].lastbase && t1 < indelreadlist[i].firstbase+20;t1++) fprintf(stdout,"%c",reflist->sequences[indelreadlist[i].chrom][t1]); }
			else{ for (t1=indelreadlist[i].firstbase;t1<indelreadlist[i].lastbase;t1++) fprintf(stdout,"%c",indelreadlist[i].sequence[t1]); }
			fprintf(stdout," %c %d %d %s,%c,%s,%d ",indelreadlist[i].strand,indelreadlist[i].L,indelreadlist[i].R,indelreadlist[i].readid,indelreadlist[i].matestrand,reflist->names[indelreadlist[i].matechrom],indelreadlist[i].mateposition);
			fprintf(stdout,"%s\n",indelreadlist[i].sequence);
		}
			// partially mapped reads 
			if (indelreadlist[i].position == indelreadlist[j].position && strcmp(indelreadlist[i].code,indelreadlist[j].code)==0 && strcmp(indelreadlist[i].code,"partial_R") ==0)
			{
				match =0; mism =0; t1 = indelreadlist[i].L-1; t2 = indelreadlist[j].L-1; 
				while (t1 > 0 && t2 > 0)
				{
					if (indelreadlist[i].sequence[t1] == indelreadlist[j].sequence[t2]) match++; else mism++; 
					t1--; t2--;
				}
				if ( (match-mism >= 6 && mism <= 1)  || (match-mism >= 12 && mism <=2))
				{
					cluster = indelreadlist[i].cluster;
					indelreadlist[cluster].readsf += indelreadlist[j].readsf; indelreadlist[j].readsf = 0;
					indelreadlist[cluster].readsr += indelreadlist[j].readsr; indelreadlist[j].readsr = 0;
					indelreadlist[cluster].newreads += indelreadlist[j].newreads; indelreadlist[j].newreads = 0;
					indelreadlist[j].cluster = cluster;
					//printf("matching partial_R mapped reads %d %d pos %d \n",match,mism,indelreadlist[i].position);
				}
			}
			if (indelreadlist[i].position == indelreadlist[j].position && strcmp(indelreadlist[i].code,indelreadlist[j].code)==0 && strcmp(indelreadlist[i].code,"partial_L") ==0)
			{
				match =0; mism =0; t1 = indelreadlist[i].L; t2 = indelreadlist[j].L; 
				while (t1 < indelreadlist[i].rl && t2 < indelreadlist[j].rl)
				{
					if (indelreadlist[i].sequence[t1] == indelreadlist[j].sequence[t2]) match++; else mism++; 
					t1++; t2++;
				}
				if ( (match-mism >= 6 && mism <= 1)  || (match-mism >= 12 && mism <=2))
				{
					cluster = indelreadlist[i].cluster;
					indelreadlist[cluster].readsf += indelreadlist[j].readsf; indelreadlist[j].readsf = 0;
					indelreadlist[cluster].readsr += indelreadlist[j].readsr; indelreadlist[j].readsr = 0;
					indelreadlist[cluster].newreads += indelreadlist[j].newreads; indelreadlist[j].newreads = 0;
					indelreadlist[j].cluster = cluster;
					//printf("matching partial_L mapped reads %d %d pos %d \n",match,mism,indelreadlist[i].position);
				}
			}
		// only print clusters corresponding to reads mapped partially, if partial is set to 1
		if (indelreadlist[i].partial == '1' && indelreadlist[i].readsf + indelreadlist[i].readsr >=2 && partial ==1)
		{
			fprintf(stdout,"PARINDEL:%s %s %d %d,%d ",indelreadlist[i].code,reflist->names[indelreadlist[i].chrom],indelreadlist[i].position+1,indelreadlist[i].readsf,indelreadlist[i].readsr);
			fprintf(stdout,"%d %d %s ",indelreadlist[i].L,indelreadlist[i].match,indelreadlist[i].readid);
			for (j=0;j<indelreadlist[i].L;j++) fprintf(stdout,"%c",indelreadlist[i].sequence[j]);
			fprintf(stdout,":"); for (j=indelreadlist[i].L;j<indelreadlist[i].rl;j++) fprintf(stdout,"%c",indelreadlist[i].sequence[j]);
			fprintf(stdout,"\n");
			continue;
		}



		if (indelreadlist[i].length < 0)  // print deletion 
		{
			//if (indelreadlist[i].length >= -25) 
			{
				for (t1=indelreadlist[i].firstbase;t1<indelreadlist[i].lastbase;t1++) fprintf(outfp,"%c",reflist->sequences[indelreadlist[i].chrom][t1]); 
			}
		//	else // for long deletions, we were printing compact deletion -> disable for now
			{
		//		for (t1=indelreadlist[i].firstbase;t1<indelreadlist[i].firstbase+5;t1++) fprintf(outfp,"%c",reflist->sequences[indelreadlist[i].chrom][t1]); 
		//		fprintf(outfp,"<%d>",-1*indelreadlist[i].length-10);
		//		for (t1=indelreadlist[i].lastbase-5;t1<indelreadlist[i].lastbase;t1++) fprintf(outfp,"%c",reflist->sequences[indelreadlist[i].chrom][t1]); 
			}
		}
		else
		{ 
			for (t1=indelreadlist[i].firstbase;t1<indelreadlist[i].lastbase;t1++) fprintf(outfp,"%c",indelreadlist[i].sequence[t1]); 
		}
		fprintf(outfp," %c %d %d %s,%c,%s,%d ",indelreadlist[i].strand,indelreadlist[i].L,indelreadlist[i].R,indelreadlist[i].readid,indelreadlist[i].matestrand,reflist->names[indelreadlist[i].matechrom],indelreadlist[i].mateposition);
		fprintf(outfp,"%s ",indelreadlist[i].sequence);
		//                fprintf(stdout,"%d ",indelreadlist[i].ambiguity);
		//if (indelreadlist[i].code[0] == 'G') fprintf(stdout,"%s",indelreadlist[i].cigar); 
		fprintf(outfp,"\n");
 */






