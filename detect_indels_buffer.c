
int SINGLE_READS = 0;
int MAXDEL = 5000;
int LOWQV = 2; float LOWQVfrac = 0.1;
int HASHTABLES = 2; // 1 table for now
int SPLITMAP = 0;
int chunksize = 0;

#include "findhits.c" // functions for finding hits and clustering them
#include "clumphits.c"
#include "splitreadmap.c" // functions to do split-read mapping of single read with mate mapped
#include "bam_realign.c"

// global variables that point to kmertables used for split-read mapping
KMERTABLE* kmertable; KMERTABLE* kmertable1;

// data structure for indel: uint32_t position, uint16_t type, uint16_t counts[2], list of pointers to reads with indel...

struct BAM_record
{
	bam1_t* b;  
	uint32_t newpos; 
	int8_t realigned; char tag; uint8_t matequal; 
	//int32_t mate; // index for mate 
	uint32_t oldpos;
	//struct realignment* RL; list of alternate alignments... newpos, newcigar, cigs, next alignment;;;
	uint32_t* newcigar; uint32_t cigs; 
	//uint32_t newIS; // b->core.pos/tid/mpos/isize has all the information we need
	// keep info for additional mappings from split-read alignments and scores for each | scores change depending on variant list 
	// EM algorithm to maximize likelihood of variants + read-alignments....
};

void read_stats(struct alignedread* read1,int* lq1,int* missing1)
{
        *lq1 = 0; *missing1 = 0; int i=0;
        for (i=0;i<read1->readlength;i++)
        {
                if (read1->sequence[i] == 'N') (*missing1)++; if (read1->quality[i] <= LOWQV+QVoffset) (*lq1)++;
        }
}

// mark reads that are candidates for realignment using split-read mapping code 
void extract_candidatereads_buffer(struct BAM_record** readbuffer,int breads,REFLIST* reflist,samfile_t *fp)
{
	fprintf(stderr,"marking candidate reads for split-read mapping\n");
	int i=0,lq=0,missing=0;	
	struct alignedread* read = calloc(1,sizeof(struct alignedread));       
        bam1_t* b;
	int spreads[3] = {0,0,0}; 

	for (i=0;i<breads;i++)
	{
		if (readbuffer[i]->realigned >=0) continue; // read in buffer from previous round, do not re-analyze 
		readbuffer[i]->tag = '-'; 
		if ( (readbuffer[i]->b->core.flag & 1) ==0 && SINGLE_READS ==0) continue; // discard single-end reads 
		if ((readbuffer[i]->b->core.flag & 12) == 12) continue; // both ends are not mapped
                b = readbuffer[i]->b;
                fetch_func(b, fp,read); 

		// should we include extra case where read is mapped but with low mapping quality but its mate is mapped with high MQ ??
		// insert size is not very good...
		
		if ( (b->core.flag & 4) == 0 && b->core.qual >= MIN_MQ) // mapped with good mapping quality 
		{
			read_stats(read,&lq,&missing);
                        if ((read->clipped-read->XC >= 5 || read->gaps >1) && read->IS != 0 && abs(read->IS) < MAX_IS + MAXDEL+read->readlength)
			{
				if ((float)missing/read->readlength <=0.05) { readbuffer[i]->tag = 'C'; spreads[1]++; } 
			}
			else if (read->gaps > 0 && lq <= 0.5*(float)read->readlength) {	readbuffer[i]->tag = 'G'; spreads[0]++; }
			readbuffer[i]->matequal = b->core.qual; // same as read ??
		}
		else if ((b->core.flag & 4) == 4  && (b->core.flag & 8) == 0 && ( b->core.tid < 0 || b->core.tid == b->core.mtid))  // unmapped read 
		{
			read_stats(read,&lq,&missing); b->core.tid = b->core.mtid; 
			if ((b->core.flag & 64) == 64) // read is first in pair, mate should be next read in list
			{
				if (i+1 < breads && strcmp(bam1_qname(readbuffer[i]->b),bam1_qname(readbuffer[i+1]->b)) == 0)
				{
					readbuffer[i]->tag = 'X'; readbuffer[i]->matequal = readbuffer[i+1]->b->core.qual; 
					if (readbuffer[i]->matequal >= MIN_MQ) spreads[2]++;
				}
			}
			else
			{
				if (i-1 >= 0 && strcmp(bam1_qname(readbuffer[i]->b),bam1_qname(readbuffer[i-1]->b)) == 0)
				{
					readbuffer[i]->tag = 'X'; readbuffer[i]->matequal = readbuffer[i-1]->b->core.qual; 
					if (readbuffer[i]->matequal >= MIN_MQ) spreads[2]++;
					//uint32_t* cigar = bam1_cigar(b); for (int k=0;k<b->core.n_cigar;k++) fprintf(stdout,"%d%c:",cigar[k]>>4,INT_CIGAROP[cigar[k]&0xf]);
					//fprintf(stdout,"%s flag %d %d %d %d %d\n",read->readid,b->core.flag,b->core.pos,b->core.isize,b->core.qual,readbuffer[i]->matequal);
				}
			}
		}
		free_readmemory(read);
	}
	free(read);
	fprintf(stderr,"candidates for split read mapping: gapped %d clipped %d unmapped %d\n",spreads[0],spreads[1],spreads[2]);
}

int splitreadanalysis_direct(struct BAM_record** readbuffer,int breads,REFLIST* reflist,samfile_t *fp)
{
	int i=0,j=0;
	struct alignedread* read = (struct alignedread*)malloc(sizeof(struct alignedread));
        unsigned char* shortseq = (unsigned char*)malloc(kmertable->kmer+1); shortseq[kmertable->kmer] = '\0';
	MATCH* fmlist1 = (MATCH*)malloc(sizeof(MATCH)*MAXMATCHES); for (i=0;i<MAXMATCHES;i++) fmlist1[i].cigarlist = NULL;

	int tablesize = 1<<(2*kmertable->kmer);  int tablesize1 = 1<<(2*kmertable1->kmer);
        uint32_t* fcigarlist = calloc(4096,sizeof(uint32_t));
        int windows_analyzed=0;
        int newmapping =0; // 1 if split read mapping finds a better alignment than original read alignment
	int chromosome = -1;

        int start=0,end=0; int reads_window =0;
	int reads_splitmapped=0;
	
	for (i=0;i<breads;i++) 
	{
		if (chromosome == -1 && (readbuffer[i]->b->core.flag & 4) == 0) 
		{
			chromosome = readbuffer[i]->b->core.tid; 
			fprintf(stderr,"chromosome %d \n",chromosome);
		}

		if (readbuffer[i]->realigned >=0 || readbuffer[i]->tag == '-') continue; 
		// update hashtable to next window if mateposition exceeds 'end'
		// add condition to generate hashtable block if read is unmapped and mate is mapped...
		if (SPLITMAP ==1 && (readbuffer[i]->b->core.flag & 4) == 0 && readbuffer[i]->b->core.pos + MAXDEL + 2*MAX_IS > end && end < reflist->lengths[chromosome])
		{
                        start = readbuffer[i]->b->core.pos -2*MAX_IS -MAXDEL; if (start < 0) start =0;
                        end = start + chunksize; if (end > reflist->lengths[chromosome]) end = reflist->lengths[chromosome];
                        for (j=0;j<tablesize;j++) { kmertable->counts[j] = 0; kmertable->index[j] = -1; }
                        for (j=0;j<chunksize;j++) { kmertable->locarray[j][0] = kmertable->locarray[j][1] = kmertable->locarray[j][2] = -1; }
                        build_kmertable(reflist,kmertable->kmer,kmertable,0,chromosome,start,end);
                        if (HASHTABLES ==2)
                        {
                                for (j=0;j<tablesize1;j++) { kmertable1->counts[j] = 0; kmertable1->index[j] = -1; }
                                for (j=0;j<chunksize;j++) { kmertable1->locarray[j][0] = kmertable1->locarray[j][1] = kmertable1->locarray[j][2] = -1; }
                                build_kmertable(reflist,kmertable1->kmer,kmertable1,kmertable1->gap,chromosome,start,end);
                        }
                        if (windows_analyzed%100 ==0) fprintf(stderr,"reads analyzed %d\n building kmerindex for chrom %s:%d-%d ",reads_window,reflist->names[chromosome],start,end);
                        reads_window =0; windows_analyzed++;
		}
		if (readbuffer[i]->tag == '-' || ((readbuffer[i]->b->core.flag & (BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)) > 0) || (readbuffer[i]->tag == 'X' && readbuffer[i]->matequal <= MIN_MQ)) continue;

		fetch_func(readbuffer[i]->b, fp,read);  // should be coupled to free_readmemory(read) call
		if (readbuffer[i]->tag != 'X') parse_cigar(read,reflist,fcigarlist);
		newmapping =0;
		if (readbuffer[i]->tag != 'G' && SPLITMAP ==1) 
		{
			if (read->mateposition >= reflist->lengths[read->tid]) 
			{
				fprintf(stderr,"error %d %d %d flag %d \n",read->mateposition,read->tid,reflist->lengths[read->tid],read->flag);
			}
			newmapping  = splitreadmapping(read,kmertable,kmertable1,reflist,shortseq,fmlist1);
		}
		if (newmapping ==1) 
		{
			readbuffer[i]->cigs = read->cigs; readbuffer[i]->newcigar = calloc(read->cigs,sizeof(uint32_t));
			for (j=0;j<read->cigs;j++) readbuffer[i]->newcigar[j] = read->cigarlist[j]; 
			// newpos and realigned tag should also be set | different from realignment tag
			reads_splitmapped++;
		}
		//if (newmapping ==0) parse_gapped_read(read,reflist,read->tid,indelreadlist,&indelreads,'G');
		//else parse_gapped_read(read,reflist,read->tid,indelreadlist,&indelreads,'S');
		//create new cigar for mapped read and copy it to data structure 
		free_readmemory(read); 
		reads_window++;
	}
	free(read); free(fmlist1);free(fcigarlist);
	return reads_splitmapped;
}

// some reads that were realigned could be left in buffer after cleanup, avoid processing again using realigned tag
// should PCR duplicate reads be realigned as well -> since if we do not their cigar will differ from the duplicate copy 07/18/13
// chrom is index in hashtable 'ht' used  for VCF index not to be confused with the index for reflist->sequences[] read->tid

int realign_reads_buffer(struct BAM_record** readbuffer,int breads,int chrom,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist,samfile_t *fp)
{
	int i=0,realigned=0,clipped_ins=0,reads_realigned=0,reads_clipped=0;
	bam1_t* b; 
	struct alignedread* read = calloc(1,sizeof(struct alignedread));       uint32_t* fcigarlist = calloc(4096,sizeof(uint32_t));
	struct REALIGNMENT* realignments = calloc(MAX_REALIGNMENTS,sizeof(struct REALIGNMENT)); 
	for (i=0;i<MAX_REALIGNMENTS;i++) realignments[i].cigarlist = calloc(256,sizeof(uint32_t)); 

	for (i=0;i<breads;i++)
	{
		if (readbuffer[i]->realigned >=0) continue; 
		b = readbuffer[i]->b; 	
		if (b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)) continue; // do not realign unmapped/dups
		if (b->core.qual >= MIN_MQ && chrom >= 0 && b->core.tid >=0) 
		{ 
			fetch_func(b, fp,read); 
			clipped_ins =0; if (CLIP_INSERTIONS ==1 || CLIP_INSERTIONS ==0) clipped_ins = modify_insertion_cigar(read);
			//if (read->gaps > 1) combine_indelcigars(read,reflist); // for multiple indels per read
			parse_cigar(read,reflist,fcigarlist); // read->mismatches is calculated in this function
		
			realigned = -1;
			if (read->clipped > read->XC || read->mismatches >= 1 || read->gaps > 1 || read->mismatches + read->gaps >=2) 
			{
				realigned = extract_variants_read(read,chromvars,varlist,chrom,reflist,realignments);
			}
			readbuffer[i]->realigned = 0;
			if (realigned >0)
			{
				//if ((b->core.flag & 8) == 8) fprintf(stderr,"single end read realigned \n");
				// edit cigarlist to add 'H' ops at start and end if needed 
				update_bam_record(b,b->core.pos,realignments[realigned].cigarlist,realignments[realigned].cigs,1); 
				readbuffer[i]->newpos = realignments[realigned].newpos; readbuffer[i]->realigned = 1; 
				reads_realigned +=1;
			}
			else if (clipped_ins > 0) 
			{
				update_bam_record(b,b->core.pos,(uint32_t*)read->cigarlist,read->cigs,0); 
				readbuffer[i]->newpos = read->position; readbuffer[i]->realigned = 1; 
				reads_clipped++;
			}
			free_readmemory(read);
		}
	}
	free(read); free(fcigarlist);
	for (i=0;i<MAX_REALIGNMENTS;i++) free(realignments[i].cigarlist); free(realignments); 
	return reads_realigned;
}

int indel_analysis_buffer(struct BAM_record** readbuffer,int breads,int chrom,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist,samfile_t *fp)
{
	extract_candidatereads_buffer(readbuffer,breads,reflist,fp);
	splitreadanalysis_direct(readbuffer,breads,reflist,fp); 
	//realign_reads_buffer(readbuffer,breads,chrom,chromvars,varlist,reflist,fp);
	return 1;
}
