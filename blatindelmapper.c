/* code for indel split-read-mapping similar to BLAT, translated from SAM input to BAM input format on january 5 2013, 
 author Vikas Bansal
*/

#include<math.h>
#include<time.h>
#include<ctype.h>
#include "bamread.h"
#include "kmertable.h"
#include "indelfunctions.h"
#include "../readfasta.h"
#include "print_options.c"

int MAXDEL= 5000;
int HASHTABLES = 2; // 1 table for now
int MAX_IS = 500; // could be different for each read, so allow it to be input from the fastq file
int PRINT_PARTIAL = 0;
int PFLAG =0;
int LOWQV = 2; int MRLENGTH = 50;
int MIN_MQ = 20; float LOWQVfrac = 0.1;
int SPLITMAP = 0; 
int SINGLE_READS = 0;
int TARGETED = 0; 
char* bam_regions;

#include "findhits.c" // functions for finding hits and clustering them
#include "clumphits.c"
#include "splitreadmap.c" // function to do split-read mapping of single read wit mate mapped

void read_stats(struct alignedread* read1,int* lq1,int* missing1) 
{
	*lq1 = 0; *missing1 = 0; int i=0;
	for (i=0;i<read1->readlength;i++) 
	{ 
		if (read1->sequence[i] == 'N') (*missing1)++; if (read1->quality[i] <= LOWQV+QVoffset) (*lq1)++; 
	}
}

// read could have an indel or it could be clipped or both, how to handle both case may 18 2012
//for (i=read->readlength-kmertable->kmer;i<read->readlength;i++) shortseq[i-l+kmertable->kmer] = sequence[i]; 
// standard sorted bam file, mates are not together but unmapped mates should be next to mapped mate 
// the only value we need is the mapping quality of the mates of unmapped reads 
 // CRUCIAL: need to set read->mtid and read->mateposition since these are used to allocate new kmertable
int getnextread(samfile_t *fp,bam1_t* b1,bam1_t* b2,struct alignedread* read1,char* alignment)
{
	int readfound = 0; int i=0,missing1=0,lq1 = 0;
	static int matemq = 0;  // declared static so that value of previous read's mapping quality is retained on next call
	while (readfound ==0 && samread(fp,b1) >=0)
	{
		if ( (b1->core.flag & 1) ==0 && SINGLE_READS ==0) continue; // discard single end reads 
		fetch_func(b1,fp,read1); // this allocates memory that needs to be freed 
	
		alignment[0] = '-'; alignment[1] = '\0'; 	
		if ((read1->flag & 4) ==0 && read1->mquality >= MIN_MQ) // read is mapped and has high mapping quality
		{
			read_stats(read1,&lq1,&missing1);  
			// clipped but treating as unmapped, should enforce constraint on matepair being mapped as well jan 5 2013
			// require the insert size to be within the maxlimit, important
                        if ((read1->clipped-read1->XC >= 5 || read1->gaps >1) && abs(read1->IS) > 0 && abs(read1->IS) < MAX_IS + MAXDEL+read1->readlength)
			{
				//if (read1->tid >= 0) fprintf(stderr,"BUG flag %s %d %s %d\n",read1->readid,read1->position,read1->sequence,read1->IS);
				//read1->mtid = read1->tid; // if the mate is unmapped, don't want to cause error
				if ((float)missing1/read1->readlength <= 0.05)
				{
					strcpy(alignment,"CLIPPED"); readfound = 1; 
				}
			}
			else if (read1->gaps > 0 && lq1 <= 0.5*read1->readlength) 
			{
				//read1->mtid = read1->tid; 
				read1->mateposition = read1->position; // if the mate is unmapped, don't want to cause error
				strcpy(alignment,"GAPPED"); readfound = 1;
				//fprintf(stdout,"%s\t%s\t%c\t%d\t%s\t%d\tGAPPED:%s:%d\n",read1->readid,read1->sequence,read1->strand,read1->chrom,read1->position,read1->mquality,read1->matechrom,read1->mateposition);
			}
		}
		else if ( (read1->flag & 4) ==4 && (read1->flag & 8) ==0 && (read1->tid < 0 || read1->tid == read1->mtid)) // read is unmapped but mate is mapped
		{
			read_stats(read1,&lq1,&missing1);
			// read1 is unmapped and first in pair, the mate should be the next read...read it to get relevant information
			// otherwise matemq should already be set from previous read, should explicity check this
			if ((read1->flag & 64) == 64)  
			{ 
				if (samread(fp,b2) >=0)  matemq = b2->core.qual; 
				else { free_readmemory(read1); return 0; } 
			}
			if ((float)missing1/read1->readlength <= 0.05 && matemq >= MIN_MQ) 
			{
			//	fprintf(stderr,"unmapped read \n");
				read1->strand = '+'; if ((read1->flag & 32) == 32) read1->strand = '-'; // mate strand
				read1->mquality = matemq;  strcpy(alignment,"X"); 	
				readfound =1;
				read1->tid = read1->mtid; // assign its mate 
				//read1->mtid = read1->tid;  // added 06/19/13 due to BUG
				// read has flag 69, -> unmapped but 
				//if (read1->tid >= 0 && read1->tid != read1->mtid) fprintf(stderr,"BUG flag %d %d %d \n",read1->flag,read1->tid,read1->mtid);
			}
		}

		matemq = read1->mquality; //strcpy(matecigar,read1->cigar);
		if (readfound ==1) return 1; else free_readmemory(read1);
	}
	return readfound;
}

// store all gapped aligned reads in array and potential indels found as well.....
// reduce penalty of gap in smith waterman alignment for reads | combine gapped, split-read, realignment, multiple realignment...
// multiple possible candidate cigars for each read 
// use buffer to store bam reads in window -> makes finding pairs easier... | de novo assembly easier... | long insertions 
// do realignment at same time -> do split-read alignment | local de novo assembly at same time.. single pass 

int splitreadanalysis_direct(char* bamfile,REFLIST* reflist,KMERTABLE* kmertable,KMERTABLE* kmertable1,int chunksize,FILE* outfp,FILE* fastafp)
{
	struct alignedread* read = (struct alignedread*)malloc(sizeof(struct alignedread));
	int i=0,j=0,lines=0;
	unsigned char* shortseq = (unsigned char*)malloc(kmertable->kmer+1); shortseq[kmertable->kmer] = '\0';

	int MAX_READS = 10000000; int indelreads =0;
	/* data structures for split-read calling */
	MATCH* fmlist1 = (MATCH*)malloc(sizeof(MATCH)*MAXMATCHES); for (i=0;i<MAXMATCHES;i++) fmlist1[i].cigarlist = NULL; 
	SPINDEL* indelreadlist = (SPINDEL*)malloc(sizeof(SPINDEL)*MAX_READS); 

	int tablesize = 1<<(2*kmertable->kmer);  int tablesize1 = 1<<(2*kmertable1->kmer); 
	int chromosome=-1,start=0,end=0; int reads_window =0;
	char alignment[256];         int* fcigarlist = (int*)malloc(sizeof(int)*4096);
	int windows_analyzed=0;
	int newmapping =0; // 1 if split read mapping finds a better alignment than original read alignment
	

	samfile_t *fp;
        if ((fp = samopen(bamfile, "rb", 0)) == 0)
        {
                fprintf(stderr, "Fail to open BAM file %s\n", bamfile); return -1;
        }
        bam1_t *b1 = bam_init1();  bam1_t *b2 = bam_init1();  
	// need to check that bam header matches that of reference fasta header...

	while (getnextread(fp,b1,b2,read,alignment) == 1)  // get next candidate read for split-read mapping, if read is unmapped, read->tid set using mate
	{
		//if (++lines%1000000 ==0)  fprintf(stderr,"processed %d reads\n",lines);
		if (read->tid != chromosome && indelreads > 0) 
		{
			fprintf(stderr,"new chromosome %s, clustering and printing indels %d\n",reflist->names[chromosome],indelreads);
			cluster_indelreads(indelreadlist,indelreads,reflist,PRINT_PARTIAL);    
			print_clusters(indelreadlist,indelreads,reflist,PRINT_PARTIAL,outfp);
			indelreads = 0;
		}
		if (read->tid != chromosome) 
		{
			// free previous chromosome data 
			if (chromosome >=0) free(reflist->sequences[chromosome]);
			if (chromosome >=0) fprintf(stderr,"free mem from prev. chrom %s new %d:%s\n",reflist->names[chromosome],read->tid,reflist->names[read->tid]);
			read_chromosome(reflist,read->tid,fastafp);
			chromosome =  read->tid; end = -1; // for loop below
		}	
		// if the current read mateposition exceeds range of current tableindex, create new table 
		// this can cause bug if mateposition is too big !! Jan 8 2013
		if (read->mateposition + MAXDEL + MAX_IS > end && end < reflist->lengths[chromosome] && SPLITMAP ==1)
		{
			if (windows_analyzed%50 ==0) fprintf(stderr,"reads analyzed %d ",reads_window);
			start = read->mateposition -MAX_IS -MAXDEL; if (start < 0) start =0;
			end = start + chunksize; if (end > reflist->lengths[chromosome]) end = reflist->lengths[chromosome]; 
			for (i=0;i<tablesize;i++) { kmertable->counts[i] = 0; kmertable->index[i] = -1; }
			for (i=0;i<chunksize;i++) { kmertable->locarray[i][0] = kmertable->locarray[i][1] = kmertable->locarray[i][2] = -1; } 
			build_kmertable(reflist,kmertable->kmer,kmertable,0,chromosome,start,end);
			if (HASHTABLES ==2)
			{
				for (i=0;i<tablesize1;i++) { kmertable1->counts[i] = 0; kmertable1->index[i] = -1; }
				for (i=0;i<chunksize;i++) { kmertable1->locarray[i][0] = kmertable1->locarray[i][1] = kmertable1->locarray[i][2] = -1; } 
				build_kmertable(reflist,kmertable1->kmer,kmertable1,kmertable1->gap,chromosome,start,end);
			}
			if (windows_analyzed%50 ==0) fprintf(stderr,"building kmerindex for chromosome %s length %d window %d-%d \n",reflist->names[chromosome],reflist->lengths[chromosome],start,end);
			reads_window =0; windows_analyzed++;
		}

		//printf("GAP %s:%d:%s %s,%c,%s %s\n",read->matechrom,read->mateposition,read->cigar,read->readid,read->strand,alignment,read->sequence);
		if (alignment[0] == '-' || (read->flag & (BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)) ==1 || read->mquality < MIN_MQ ) { } // filter mark duplicate reads and low-mapping quality 
		else 
		{
			if (alignment[0] != 'X') parse_cigar(read,reflist,fcigarlist); 
			
			newmapping =0;
			if (alignment[0] != 'G' && SPLITMAP ==1) 
			{
				newmapping  = splitreadmapping(read,kmertable,kmertable1,reflist,shortseq,fmlist1);
			}
			// once the read is aligned using split mapping, copy cigar list and parse it 
			if (newmapping ==0) parse_gapped_read(read,reflist,read->tid,indelreadlist,&indelreads,'G'); 
			else parse_gapped_read(read,reflist,read->tid,indelreadlist,&indelreads,'S');
		}
		free_readmemory(read); // free read memory that was allocated in function getnextread
		reads_window++;
	}
	// cluster the list of indels and then print each cluster, last argument is to enable/disable printing of partially aligned reads.
	cluster_indelreads(indelreadlist,indelreads,reflist,PRINT_PARTIAL);
	print_clusters(indelreadlist,indelreads,reflist,PRINT_PARTIAL,outfp);

	if (chromosome >=0) free(reflist->sequences[chromosome]);
        bam_destroy1(b1); bam_destroy1(b2); samclose(fp);
	return 1;
}

int main(int argc, char* argv[])
{
	char outfile[1024]; char bamfile[1024]; char fastafile[1024]; char bedfile[1024]; bam_regions = NULL;
	strcpy(bamfile,"None"); strcpy(fastafile,"None"); strcpy(outfile,"None"); strcpy(bedfile,"None");
	int i=0;
	for (i=1;i<argc;i+=2)
	{
		if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0)        strcpy(bamfile,argv[i+1]);
		if (strcmp(argv[i],"--out") ==0 || strcmp(argv[i],"--VCF") ==0)        strcpy(outfile,argv[i+1]);
		else if (strcmp(argv[i],"--reffile") ==0 || strcmp(argv[i],"--ref") ==0)        strcpy(fastafile,argv[i+1]);
		else if (strcmp(argv[i],"--mmq") ==0)       MIN_MQ = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--maxIS") ==0)       MAX_IS = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--maxdelsize") ==0)       MAXDEL = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--qvoffset") ==0)       QVoffset = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--pflag") ==0)       PFLAG = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--split") ==0)       SPLITMAP = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--singlereads") ==0)       SINGLE_READS = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--minreads") ==0)       MIN_READS_FOR_INDEL = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--bed") ==0)   strcpy(bedfile,argv[i+1]);
		
		else if (strcmp(argv[i],"--regions") ==0)       
		{
			bam_regions = (char*)malloc(strlen(argv[i+1])+1);   strcpy(bam_regions,argv[i+1]);
		}
		else if (strcmp(argv[i],"--complexindels") ==0 || strcmp(argv[i],"--complex") ==0 || strcmp(argv[i],"--complexindel") ==0)       
		{
			COMPLEX_INDELS = atoi(argv[i+1]);
			fprintf(stderr,"complex indels or indels close to each other will be output as single haplotype \n");
		}
	}

	if (strcmp(bamfile,"None") ==0 || strcmp(fastafile,"None") ==0)
	{
		print_options_indelmapper(); return -1;
	}
	REFLIST* reflist = NULL; reflist = init_reflist(fastafile,reflist); 
	int bh = validate_bam_header(bamfile,reflist); if (bh ==0) return 0; 

	TARGETED = 0; if (strcmp(bedfile,"None") != 0 && read_bedfile(bedfile,reflist) != -1) TARGETED= 1;
	

        FILE* fastafp = fopen(fastafile,"r");

	//if (reflist != NULL) read_fasta(fastafile,reflist); 
	fprintf(stderr,"finished reading reference sequence \n");

	FILE* outfp = stdout; if (strcmp(outfile,"None") !=0) outfp = fopen(outfile,"wb"); 

	// maximum of two hashtables. // chunksize can be small since complexity is independent of chunksize !!
        //int kmer = 10; int kmer1 = 8; int chunksize = 200000; 
	int kmer = 9; int kmer1 = 8; int chunksize = 65536; // 4^8
	KMERTABLE* kmertable = (KMERTABLE*)malloc(sizeof(KMERTABLE)); kmertable->gap = 0; kmertable->kmer = kmer;
	KMERTABLE* kmertable1 = (KMERTABLE*)malloc(sizeof(KMERTABLE)); kmertable1->gap = 4; kmertable1->kmer = kmer1;
	//int tablesize = 1<<(2*kmer); double maxhitsforseed = (double)(chunksize+100000)*16/tablesize; MAX_KMER_HITS =(int)maxhitsforseed;
	int padding = 5*(MAX_IS+MAXDEL); // to account for insert size and max deletion size

	if (SPLITMAP ==1)  // perform split read alignment 
	{
		init_kmertable(reflist,kmer,kmertable,chunksize+padding);
		if (HASHTABLES ==2) init_kmertable(reflist,kmer1,kmertable1,chunksize+padding);
	}
	splitreadanalysis_direct(bamfile,reflist,kmertable,kmertable1,chunksize,outfp,fastafp); 
	if (outfp !=0) fclose(outfp);  
	fclose(fastafp);
	return 1;
}

// ./testmapper --bam /projects/stsi3/PROJECTS/T2D-pooledseq/T2D-pooledseq-mergedbams/G1P1.merged.picard.bam --ref /projects/stsi/vbansal/T2D-pooledseq-july2012/target/ncbi37.fa --out G1P1.indels.vcf --complex 1 > G1P1.indels

