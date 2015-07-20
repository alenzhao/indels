#include "bamread.h"
char INT_CIGAROP[] = {'M','I','D','N','S','H','P','E','X'};

int STRICT_REFMATCH =0;
int QVoffset = 33;

static unsigned int RCtable[256] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 84, 4, 71,  4, 4, 4, 67,  4, 4, 4, 4,  4, 4, 78, 4, 4, 4, 4, 4,  65, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 84, 4, 71,  4, 4, 4, 67,  4, 4, 4, 4,  4, 4, 78, 4, 4, 4, 4, 4,  65, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

#include "indels/indelcigar.c" // shared indel functions 

// works in-situ on same sequence
void RC(char* sequence,char* revsequence,int l)
{
        int i=0; char ch;  // should work on same sequence as well
	for (i=0;i<l/2;i++) 
	{
		ch = RCtable[(int)sequence[l-i-1]]; sequence[l-i-1] = RCtable[(int)sequence[i]]; sequence[i] = ch; 
		//revsequence[i] = RCtable[(int)sequence[l-i-1]]; 
	}
	revsequence[l] = '\0';
}

// size of sequence and cigarlist are bounded | function is unused
void read_init(struct alignedread* read)
{
	read->sequence = calloc(256,1); read->quality = calloc(256,1); 
	read->cigarlist = calloc(64,sizeof(uint32_t)); read->fcigarlist = calloc(256,sizeof(uint32_t));
}

// if the first or last cigar operation is 'I' (novoalign or stampy) -> change it to softclip  | june 12 2013

// first four bits store cigar operation, remaining 28 bits store the length of operation 
// populate data structure 'read' using bam record 'b', the new cigarlist does not store hard clips 
int fetch_func(const bam1_t *b, void* data, struct alignedread* read)
{
	uint32_t *cigar = bam1_cigar(b);  const bam1_core_t *c = &b->core;
        int i,op,ol,j=0;
        read->cigs =0; read->alignedbases = 0; read->clipped =0; read->span =0; read->gaps=0; read->cflag =0;
	read->readlength= b->core.l_qseq;
	read->sequence = (char*)malloc(b->core.l_qseq+1); read->quality = (char*)malloc(b->core.l_qseq+1);
        uint8_t* sequence = bam1_seq(b); uint8_t* quality = bam1_qual(b);
        for (i=0;i<b->core.l_qseq;i++) read->sequence[i] = bam_nt16_rev_table[bam1_seqi(sequence,i)]; read->sequence[i] = '\0';
        for (i=0;i<b->core.l_qseq;i++) read->quality[i] = (char)(quality[i]+33); read->quality[i] = '\0';
       	
	read->mismatches=0; 
	read->flag = c->flag; read->mquality= c->qual; read->position = c->pos; read->mateposition = c->mpos; read->IS = c->isize;
        read->strand = '+'; if ((read->flag & 16) == 16) read->strand = '-'; // fixed sept 29 2011 
	read->fcigs =0; // important to initialize to 0

	read->cigs =c->n_cigar;  
	read->cigarlist = calloc(read->cigs,sizeof(uint32_t));  // maxsize of cigarlist
	if ((cigar[0]&0xf) == BAM_CHARD_CLIP) read->cigs--; if ((cigar[c->n_cigar-1]&0xf) == BAM_CHARD_CLIP) read->cigs--;
	j=0;
        for (i = 0; i < c->n_cigar; ++i)
        {
                op = cigar[i]&0xf; ol = cigar[i]>>4;
		if (op != BAM_CHARD_CLIP) read->cigarlist[j++] = cigar[i];   // ignore HARD CLIP in list 
		//read->cigarlist[j++] = cigar[i]; // ignore HARD CLIP in list 
                if (op == BAM_CMATCH) 
		{ 
			read->alignedbases += ol;   read->span += ol; 
		}
                else if (op == BAM_CDEL) 
		{ 
			read->gaps +=1; read->span += ol;
		}
                else if (op == BAM_CINS) 
		{ 
			read->alignedbases += ol;  read->gaps += 1; 
		}
                else if (op == BAM_CREF_SKIP) read->span += ol; 
                else if (op == BAM_CSOFT_CLIP) 	read->clipped += ol;
                else if (op == BAM_CHARD_CLIP) {}
		else read->cflag = 1;
        }

	read->readid=(char*)malloc(c->l_qname+1); char* qs = (char*)b->data;
        for (i=0;i<c->l_qname;i++) read->readid[i] = qs[i]; read->readid[i]= '\0';

	// only part that requires fp/data access to be current seek in bamfile
        samfile_t *fp = (samfile_t*)data; 
	if (c->tid >= 0) read->chrom = fp->header->target_name[c->tid]; else read->chrom = NULL;
	if (c->mtid >= 0) read->matechrom = fp->header->target_name[c->mtid]; else read->matechrom = NULL;


        read->tid = c->tid; read->mtid = c->mtid;
	
	read->XC = 0; 
	uint8_t* cmstring = bam_aux_get(b,"XC");
        if (cmstring != NULL) read->XC = read->readlength - *(cmstring+1);
	//if (read->XC < read->readlength) fprintf(stdout,"%s %s %d %d\n",read->chrom,read->readid,read->position,read->XC);
        // for MAQ bam files, mtid is not set resulting in lack of paired-end reads, may 1 2012
        return 0;
}

void free_readmemory(struct alignedread* read)
{
        free(read->readid); free(read->sequence); free(read->quality); free(read->cigarlist);
	if (read->fcigs > 0) free(read->fcigarlist);
}

// need to change the 0-based position, bin (??), n_cigar, change cigar in data, add additional flags: RL:i:1, OC:i:original_cigar, OP:i:originalposition
int update_bam_record(bam1_t* b,int newpos,uint32_t* cigarlist,int32_t cigs,int flag)
{
	// we can update old bam record to make it faster 
	// copy the data portion of record after cigar to a new location and use memcpy...
	uint32_t np = b->core.pos+1; bam_aux_append(b,"OP",'i',4,(uint8_t*)&np); // RL stores the original start position of read
	// create string from original cigar to append to bam record as metadata 
	int slen = 0; int i=0, op=0,ol=0; uint32_t *cigar = bam1_cigar(b); 
	for (i=0;i<b->core.n_cigar;i++) slen += (cigar[i]>>4)/10 +2;
	char* cigarstring = calloc(1,slen+1); slen =0; 
	for (i=0;i<b->core.n_cigar;i++)
	{
                op = cigar[i]&0xf; ol = cigar[i]>>4; 
		slen += sprintf(cigarstring+slen,"%d%c",ol,INT_CIGAROP[op]);
	}
	cigarstring[slen] = '\0';
	//fprintf(stdout,"cigar %s \n",cigarstring);
	bam_aux_append(b,"OC",'Z',slen+1,(uint8_t*)cigarstring); // RL stores the original start position of read
	free(cigarstring);
	int RL =1;
	if (flag ==1) bam_aux_append(b,"RL",'i',4,(uint8_t*)&RL);
	 
	//also add info about indel to which realigned, useful for long insertions.... need to compact it for long deletions...
	// IR:Z:position_ref_alternate from VCF file
	// # of mismatches in new read XM: tag should be updated, MD tag...

	// if original cigar had 'H' ops at first or last position then update cigarlist to reflect this 06/08/13
	if ((cigar[0]&0xf) == BAM_CHARD_CLIP) // move all cigars to right by one position
	{
		for (i=cigs-1;i>=0;i--) cigarlist[i+1] = cigarlist[i]; 
		cigarlist[0] = cigar[0];  
		cigs++;
	}
	if ((cigar[b->core.n_cigar-1]&0xf) == BAM_CHARD_CLIP) cigarlist[cigs++]= cigar[b->core.n_cigar-1];
	/*
	*/
	uint32_t new_data_len = b->data_len - 4*b->core.n_cigar + 4*cigs;
	uint8_t* temp_data = calloc(1,b->data_len-b->core.l_qname-b->core.n_cigar*4); // store copy of previous data 
 	memcpy(temp_data,b->data+b->core.l_qname+b->core.n_cigar*4,b->data_len-b->core.l_qname-b->core.n_cigar*4);
	if (b->m_data < new_data_len)
	{
		b->m_data = new_data_len; kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}

	memcpy(b->data+b->core.l_qname,cigarlist,4*cigs); // cigar copied
	memcpy(b->data+b->core.l_qname + 4*cigs,temp_data,b->data_len-b->core.l_qname-b->core.n_cigar*4);
	b->core.n_cigar = cigs; 
	b->core.pos = newpos; // new values 
	b->data_len = new_data_len; // important to set this variable...	

	free(temp_data);
	return 1;

	// void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data) 
	// function that will append new TAGS to bam record: RL:i:1 original_cigar and position as well...
	// data has format: qname-cigar-seq-qual-aux 
}
//uint32_t* cigar = bold->data+bold->core.l_qname;  int i=0;
//typedef struct { bam1_core_t core; int l_aux, data_len, m_data; uint8_t *data; } bam1_t;
// typedef struct { int32_t tid;  int32_t pos;  uint32_t bin:16, qual:8, l_qname:8;  uint32_t flag:16, n_cigar:16;  int32_t l_qseq;  int32_t mtid;  int32_t mpos; int32_t isize; } bam1_core_t;  
// pointer to start of sequence:  bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)	
// define bam1_aux(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (b)->core.l_qseq + ((b)->core.l_qseq + 1)/2)
//for (i=0;i<c->n_cigar;i++) fprintf(stdout,"%d:%d \t",cigar[i]&0xf,cigar[i]>>4); fprintf(stdout," %d %d %d\n",bold->l_aux,bold->data_len,bold->m_data);

// parse cigar to calculate no of mismatches in alignment and populate fcigarlist
int parse_cigar(struct alignedread* read,REFLIST* reflist,uint32_t* fcigarlist)
{
	int i=0,t=0, l1=0,l2=0; int l=0; 
	int f=0,m=0; int op;
	//int* fcigarlist = (int*)malloc(sizeof(int)*read->cigs+sizeof(int)*(int)(read->readlength)); 
	read->fcigs =0; read->cflag =0; read->mismatches=0;  // l1 and l2 are indexes on read and reference sequence
	for (i=0;i<read->cigs;i++) 
	{
		op = read->cigarlist[i]&0xf; l = read->cigarlist[i]>>4; 
		// ignore hard clips as well (SOLID data can cause BUG) // HARD_CLIP is being ignored
		if (op != BAM_CMATCH && op != BAM_CHARD_CLIP) fcigarlist[f++] = read->cigarlist[i];
		if (op == BAM_CMATCH)
		{
			m=0;
			for (t=0;t<l;t++)
			{
				if (read->sequence[l1+t]  != reflist->sequences[read->tid][read->position+l2+t] && read->sequence[l1+t] != reflist->sequences[read->tid][read->position+l2+t]-32 && read->sequence[l1+t]  !='N')
				{
					read->mismatches++;
					if (m > 0) fcigarlist[f++] = m<<4; // M
					fcigarlist[f++] = 24; m=0; // 1X = 1<<4 + 8 = 24
					//if (m > 0) fcigarlist[f++] = m<<4;fcigarlist[f++] = (read->sequence[l1+t]<<4)+8; m=0;
				}
				else m++; 
			}
			if (m > 0) fcigarlist[f++] = m<<4;
			l1 += l; l2 +=l; 
		}
		else if (op == BAM_CDEL)
		{
			l2 += l;
		}
		else if (op == BAM_CREF_SKIP) l2 += l; 
		else if (op == BAM_CINS)
		{
			l1 += l;
		}
		else if (op == BAM_CSOFT_CLIP) l1 += l;
		else if (op != BAM_CHARD_CLIP) { read->cflag =1; break; }
	}
	read->l1 = l1; read->l2 = l2; 
	if (f > 0)
	{
		read->fcigs = f; read->fcigarlist = calloc(read->fcigs,sizeof(uint32_t));
		for (i=0;i<read->fcigs;i++) read->fcigarlist[i] = fcigarlist[i];
		if (read->clipped >=2 || read->mismatches >=3 || read->gaps >0)
		{
			//fprintf(stdout,"%s %d %s %d %s ",read->sequence,read->position,read->readid,current,reflist->names[current]);
			//for (i=0;i<f;++i) fprintf(stdout,"%d%c ",read->fcigarlist[i]>>4,INT_CIGAROP[read->fcigarlist[i]&0xf]);
			//fprintf(stdout,"\n");
		}
	}
	return 1;
	//free(fcigarlist); return 0;
	//fprintf(stdout,"mismatches %d %d %s %d %d %s\n",read->mismatches,read->NM,read->sequence,current,read->position,read->cigar);
}


int validate_bam_header(char* bamfile,REFLIST* reflist)
{
	bamFile fp1=bam_open(bamfile, "r"); bam_header_t *header;		
	if(fp1==0)
	{
		fprintf(stderr,"Cannot open %s\n",bamfile);
		return -1;
	}
	int i=0; int valid = 1; int mismatches = 0;
	header= bam_header_read(fp1); //fprintf(stderr,"targets %d \n",header->n_targets); 
	if (header->n_targets != reflist->ns) 
	{
		fprintf(stderr,"number of sequences in bam reference header (%d) does not match reference fasta file (%d) \n",header->n_targets,reflist->ns);
		fprintf(stdout,"##########################################\n");
		for (i=0;i<header->n_targets && i < reflist->ns && i < 10;i++) fprintf(stderr,"name %s %d ref %s %d\n",header->target_name[i],header->target_len[i],reflist->names[i],reflist->lengths[i]);
		for (i=0;i<header->n_targets && i < reflist->ns && i < 10;i++) 
		{
			if (strcmp(header->target_name[i],reflist->names[i]) !=0 || header->target_len[i] != reflist->lengths[i]) mismatches++;
		}
		if (mismatches > 0 || STRICT_REFMATCH ==1)
		{
			fprintf(stdout,"##########################################\n");
			fprintf(stderr,"please make sure the reference fasta file is the same as the one used to generate BAM file, mismatches %d\n\n",mismatches);
			//for (i=0;i<reflist->ns;i++) fprintf(stderr,"name %s %d \n",reflist->names[i],reflist->lengths[i]);
			valid = 0;
		}
		else
		{
			fprintf(stdout,"##########################################\n");
		        fprintf(stderr,"partial match of fasta files, continue ?\n\n");
		}
	}
	else 
	{
		for (i=0;i<header->n_targets;i++) 
		{ 
			if (header->target_len[i] != reflist->lengths[i]) valid = 0; 
		} 
	}
	bam_header_destroy(header); bam_close(fp1);
	return valid;
	
	//b1 = bam_init1(); bam_index_t* idx; bam_iter_t iter; 
	//int* ref,int* start,int* end;
	//if (bam_regions != NULL)
	{
		//if ( (idx  = bam_index_load(bamfile)) ==0)
                {
           //                     fprintf(stderr,"unable to load bam index for file %s\n",bamfile); return -1;
                }
                //bam_parse_region(bamfiles_data[i].fp->header,regions,&ref,&beg,&end);
	}

}

