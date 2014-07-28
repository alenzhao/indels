
//bam_mate.c from samtools has useful functions for doing this
void fix_insertsizes_buffer(struct BAM_record** readbuffer,int breads,int maxISIZE)
{
	// inferred insert size would be the difference between the 5' positions of the two reads
	// code to find mate for each realigned read and fix insert size for each read in pair

	// fix insert sizes should be done using original read positions, only change insert size and mate positions of reads
	int i=0,j=0,delta=0;
	for (i=0;i<breads;i++)
	{
		delta = readbuffer[i]->b->core.pos - readbuffer[i]->newpos;

		if (readbuffer[i]->realigned != 1 || (readbuffer[i]->b->core.flag & 1) ==0) continue; 
		if (readbuffer[i]->b->core.tid != readbuffer[i]->b->core.mtid && readbuffer[i]->b->core.mtid >=0) continue; // different chromosomes

		if (delta != 0 && readbuffer[i]->b->core.isize > 0 && readbuffer[i]->b->core.isize <= maxISIZE) 
		{
			j= i+1; while (j < breads && readbuffer[j]->b->core.pos <= readbuffer[i]->b->core.mpos) 
			{
				if (readbuffer[j]->b->core.pos == readbuffer[i]->b->core.mpos && strcmp(bam1_qname(readbuffer[i]->b),bam1_qname(readbuffer[j]->b)) == 0)
				{
					readbuffer[i]->b->core.isize += delta; 	readbuffer[j]->b->core.isize -= delta;
					readbuffer[j]->b->core.mpos = readbuffer[i]->newpos;
				}
				j++; 
			}
			//first read in pair is realigned at start -> start position--, IS++
		}
		else if (readbuffer[i]->b->core.isize < 0 && readbuffer[i]->b->core.isize >= -1*maxISIZE) 
		{
			// second read in pair is realigned at end -> IS++ |  go backwards in list to find mate
			j = i-1; while (j >=0 && readbuffer[j]->b->core.pos >= readbuffer[i]->b->core.mpos)
			{
				if (readbuffer[j]->b->core.pos == readbuffer[i]->b->core.mpos && strcmp(bam1_qname(readbuffer[i]->b),bam1_qname(readbuffer[j]->b)) == 0) 
                                {	
					// end_read = (old.start + span) - (old.start - newstart)
                                        int32_t end_read = bam_calend(&readbuffer[i]->b->core,bam1_cigar(readbuffer[i]->b)) - delta; 

					// what if the position of mate also changed ?? is this still correct
                                        readbuffer[i]->b->core.isize = readbuffer[j]->newpos - end_read; 
                                        readbuffer[j]->b->core.isize = end_read-readbuffer[j]->newpos;
					readbuffer[j]->b->core.mpos = readbuffer[i]->newpos;
                                }
                                j--;
			}
			
		}
		else if (delta != 0 && (readbuffer[i]->b->core.flag & 8) == 8) // when we change read position -> need to do same for mate if it is unmapped 
		{
			j = i+1; if (j < breads && strcmp(bam1_qname(readbuffer[i]->b),bam1_qname(readbuffer[j]->b)) == 0)
			{
				readbuffer[j]->b->core.pos = readbuffer[i]->newpos; readbuffer[j]->b->core.mpos = readbuffer[i]->newpos;
				//uint32_t* cigar = bam1_cigar(readbuffer[i]->b); 
				//for (k=0;k<readbuffer[i]->b->core.n_cigar;k++) fprintf(stdout,"%d%c:",cigar[k]>>4,INT_CIGAROP[cigar[k]&0xf]);
				//fprintf(stdout,"fixed insert size %s %d %d %d\n",bam1_qname(readbuffer[i]->b),readbuffer[i]->b->core.pos+1,readbuffer[i]->b->core.isize,delta);
			}
			j = i-1; if (j >= 0 && strcmp(bam1_qname(readbuffer[i]->b),bam1_qname(readbuffer[j]->b)) == 0)
			{
				readbuffer[j]->b->core.pos = readbuffer[i]->newpos; readbuffer[j]->b->core.mpos = readbuffer[i]->newpos;
			}
		}
	}
}

void sort_reads_buffer(struct BAM_record** readbuffer,int breads,int* readstofix, int* moves)
{
	struct BAM_record* temp;
	int i=0,j=0,k=0;
	*readstofix = 0; *moves=0;
	for (i=0;i<breads;i++) readbuffer[i]->b->core.pos = readbuffer[i]->newpos; // update position to newpos before sorting

	for (i=1;i<breads;i++)
	{
		if (readbuffer[i]->b->core.pos < readbuffer[i-1]->b->core.pos)
		{
			j=i-2; while (j>= 0 && readbuffer[i]->b->core.pos < readbuffer[j]->b->core.pos) j--; 
			temp = readbuffer[i]; for (k=i-1;k>j;k--) readbuffer[k+1] = readbuffer[k]; readbuffer[j+1] = temp; 
			*moves += i-j-1; (*readstofix)++;
		}
	}
	for (i=breads-2;i>=0;i--) // check if adjacent read pairs are out of order (r[i] > r[i+1]) 
	{
		if (readbuffer[i]->b->core.pos > readbuffer[i+1]->b->core.pos) 
		{
			j=i+2; while (j < breads && readbuffer[i]->b->core.pos > readbuffer[j]->b->core.pos) j++; 
			temp = readbuffer[i]; for (k=i;k<j-1;k++) readbuffer[k] = readbuffer[k+1]; readbuffer[j-1] = temp; 
			*moves += j-i-1; (*readstofix)++;
		}
	}
}
int clean_bam_buffer(struct BAM_record** readbuffer,int breads,samfile_t* fout,int newchrom) 
{
	int i=0,j=0,k=0; 
	int lastreadpos =0;
	j = breads-1;  while (j >= 0 && (readbuffer[j]->b->core.flag & 4) ==4) j--;  
	if (j >=0) lastreadpos = readbuffer[j]->b->core.pos; else return 0;

	int moves =0,readstofix =0;
	struct BAM_record* temp; 

	fix_insertsizes_buffer(readbuffer,breads,5000);
	sort_reads_buffer(readbuffer,breads,&readstofix,&moves);

	// find the first read for which the mate is not in list, yet to arrive -> cannot print to file yet
	// maximum leftshift for a read (cannot print reads that start after that position) 
	for (i=0;i<breads;i++)
	{	
		if (i > 0 && readbuffer[i]->b->core.pos < readbuffer[i-1]->b->core.pos) fprintf(stderr,"sort error %d %d %d:%d \n",i-1,i,readbuffer[i-1]->b->core.pos,readbuffer[i]->b->core.pos);
		if ( newchrom == 0 && ( (readbuffer[i]->b)->core.flag & 4) ==0 && (readbuffer[i]->b)->core.pos >= lastreadpos-3000) 
		{
			fprintf(stderr,"i %d readpos %d last pos %d \n",i,(readbuffer[i]->b)->core.pos,lastreadpos);
			break;
		}
		if (OUTPUT_BAM==1 && fout != NULL) samwrite(fout,readbuffer[i]->b); 	
		else if (OUTPUT_BAM==2 && readbuffer[i]->realigned ==1 && fout != NULL) samwrite(fout,readbuffer[i]->b); 	
		else if (OUTPUT_BAM ==3) fprintf(stdout,"READ %d %d %d %d mate:%d RL:%d\n",(readbuffer[i]->b)->core.flag,readbuffer[i]->b->core.qual,readbuffer[i]->oldpos,readbuffer[i]->b->core.pos,readbuffer[i]->b->core.mpos,readbuffer[i]->realigned);
	}
	fprintf(stderr,"printed %d reads to bam file fixed %d/%d |  reads in buffer: %d last position %d\n\n",i,readstofix,moves,breads-i,lastreadpos);

	// move the reads from (k-breads) to top of buffer 
	for (k=0,j=i;j<breads;j++) 
	{
		temp = readbuffer[k]; readbuffer[k]= readbuffer[j]; readbuffer[j] = temp; k++;
	}
	return k;
}


