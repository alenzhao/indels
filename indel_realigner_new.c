// author: Vikas Bansal, vbansal@scripps, 2011-2013, last modified jan 4 2013, 
// code for indel realignment using bam file and VCF file with list of known indels 
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "../hashtable.h"
#include "../readfasta.h"
#include "bamread.h"
#include "sam.h"
#include "kmertable.h"
#include "indelfunctions.h" 
#include "readvariant.h"
#include "print_options.c"

/*
 ./INDELrealigner1 --bam /oasis/tscc/scratch/vbansal/T2D-pools/phase1/G1pools/G1P1.novoalign.MD.bam --ref /projects/stsi/vbansal/T2D-pooledseq-july2012/target/ncbi37.fa --VCF G1P1.indels.VCF.1 --ob 2 --out G1P1.new.bam > a
 
similar to extracthairs, sort reads by start position of 2nd read in pair | fix pairs | sort again by reverse start position | print from bottom of list until we encounter a read whose mate is not yet observed (unless it is single ended...) | repeat 
much easier for single ended reads... | skip mate pair sorting step 

we want to output bam file that is sorted and has correct mate pair information...

D3NJ6HQ1:258:D12WJACXX:1:2110:12371:78866	163	1	120612030	31	8M9D92M	=	120612156	163	| original cigar was 8S92M and position was 120613047 (17 bp decrease)

do we want to realign reads that do not span entire HP tract: or just clip last base 
D3NJ6HQ1:258:D12WJACXX:3:1109:18993:92327	99	1	66100473	70	89M1I10M	=	66100553	137	TCAGTGCATATATAGTAGAAGCCTTAAGAAA
AAAGAAAATGAGCAAGCAAATATTTGAAGAAATGTATAAAACCATAGATTTCTTTCAGAAAAAAAAAAA	EEFEEEFGFGFGFGFFGEEFEFFGGGGFGGGGGGFGGGGFFGFFGGFFGGGGGDGGFFGDGFFGEGFGGGGHFFHGGEGGGGFFFEFHFGFEFFE
ECC=<	PG:Z:novoalign	AS:i:30	UQ:i:30	NM:i:1	MD:Z:99C0	PQ:i:32	SM:i:70	AM:i:70	OQ:Z:@@CFBDFFGGGHHIJHGGEHIIIJIJGIIJJJIJGIJJJIGGIIIGIIJGIJJAGGGEFAH>DGIIEEIIIGIH
HFEHFEFFCDDAEEEDA>CB@=9@00	OP:i:66100473	OC:Z:100M	RL:i:1

do we want to place I cigar tag at end or beginning of read, not followed by a 'M' | can be useful for long insertions or equivalent to 'S'

// simplest solution: keep all realigned reads in buffer | 2-pass algorithm to fix everything | sorted list of reads (by first in pair, 2nd in pair) 
// invariant, position for unmapped reads in pair should be equal to its mate and they should be together in sorted bam 10/27/13 

// samread = bam_read1 -> only does realloc but no malloc/free | repeated calls do not require memory free
// same approach should be used for fetch_func (own data structure) to save time
	
// we can never leave pointer hanging, only swap, since memory for bam record is not deallocated
// how much buffer size do we need depends on maximum leftshift in alignment of read | alignment can also be rightshifted (possible)
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MINQ = 13; // minimum base quality
int MIN_MQ = 20; // minimum read mapping quality
int MAX_IS =  500; // maximum insert size
int MIN_IS =  0; // maximum insert size
int PEONLY = 0; // if this is set to 1, reads for which only one end is mapped are not considered for hairs 
int BSIZE = 500; int IFLAG = 0; int VARIANTS = 0; int PARSEINDELS =0;
//int QVoffset = 33; declared in samread.h
int PFLAG = 1;
int OUTPUT_BAM =1; 
int MAX_REALIGNMENTS = 100;
int CLIP_INSERTIONS = 0; // parameter for novoalign fullNW alignments, remove excessive 'I' cigars at ends of reads 
int STRINGENCY = 0; // 0/1/2 to control how much realignment is done close to ends of reads...
//int COMPLEX_INDELS = 0;

#include "detect_indels_buffer.c"
#include "clean_sort_buffer.c"
//#include "bam_realign.c" 

int realign_bamfile_sorted(char* bamfile,char* outfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist,FILE* fastafp)
{
	fprintf(stderr,"reading sorted bamfile %s \n",bamfile);
	int i=0; int chrom=0; int reads=0; int prevtid = -1; int max_buffer_size = 5000000; int nc=0;
 
	struct BAM_record* temp;
	bam1_t *b;// = bam_init1();  // bam record 

	struct BAM_record** readbuffer = calloc(max_buffer_size,sizeof(struct BAM_record*)); int breads =0; 
	for (i=0;i<max_buffer_size;i++) 
	{
		readbuffer[i] = calloc(1,sizeof(struct BAM_record)); readbuffer[i]->b = bam_init1();
		readbuffer[i]->realigned = -1; readbuffer[i]->newcigar = NULL; readbuffer[i]->cigs =0;
	}
	
	samfile_t *fp; 
	if ((fp = samopen(bamfile, "rb", 0)) == 0) 
	{ 
		fprintf(stderr, "Fail to open BAM file %s\n", bamfile); return -1; 
	}
	//bam_header_t* bamheader_obj = bam_header_read(fp); // read header of bam file  // char* samstring;
	samfile_t* fout = NULL;
	if (strcmp(outfile,"None") != 0) 
	{
		fout = samopen(outfile,"wb",fp->header); 
		if (fout ==0) fprintf(stderr,"failed to open output bam file %s \n",outfile); 
	}

	// need to handle unmapped reads (pairs) at end of sorted bam file
	while (1)
	{
		if (breads >= max_buffer_size) 
		{
			fprintf(stderr,"buffer size exceeded, need to clean buffer %d\n",breads);
			indel_analysis_buffer(readbuffer, breads,chrom,chromvars,varlist,reflist,fp);
			breads = clean_bam_buffer(readbuffer,breads,fout,0);
			fprintf(stderr,"buffer cleaned, reads left %d\n",breads);
		}
		if (samread(fp,readbuffer[breads]->b) < 0) break;
		b = readbuffer[breads]->b; 
		readbuffer[breads]->realigned = -1; readbuffer[breads]->newpos = b->core.pos;
		readbuffer[breads]->oldpos = (readbuffer[breads]->b)->core.pos; //readbuffer[breads]->mate = -1;
		if (readbuffer[breads]->cigs > 0) free(readbuffer[breads]->newcigar); 
		readbuffer[breads]->newcigar = NULL; readbuffer[breads]->cigs = 0; readbuffer[breads]->tag = '-';
		if (b->core.tid != prevtid && b->core.tid >=0)
		{
			if (prevtid >=0) 
			{
				if (breads >0)
				{ 
					indel_analysis_buffer(readbuffer, breads,chrom,chromvars,varlist,reflist,fp);
					clean_bam_buffer(readbuffer,breads,fout,1);
					temp = readbuffer[0]; readbuffer[0] = readbuffer[breads]; readbuffer[breads] = temp; 
					breads = 0; 
				}
				free(reflist->sequences[prevtid]); reflist->sequences[prevtid] = NULL; 
                        	fprintf(stderr,"free memory from chrom %s.. total reads processed %d ",reflist->names[prevtid],reads);
			} 
                        read_chromosome(reflist,b->core.tid,fastafp); prevtid = b->core.tid;
			chrom = getindex(ht,fp->header->target_name[b->core.tid]);
			//fprintf(stderr,"chrom %d length %d\n",fp->header->target_name[b->core.tid],reflist->lengths[fp->header->target_name[b->core.tid]]);
		}
		breads++;
		//if (b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)) continue;
		if (++reads%1000000 ==0) fprintf(stderr,"processed %d reads\n",reads);
	}
	if (breads > 0)
	{
		fprintf(stderr,"final cleanup of buffer with reads = %d chrom %d\n",breads,prevtid);
		indel_analysis_buffer(readbuffer, breads,chrom,chromvars,varlist,reflist,fp);
		realign_reads_buffer(readbuffer, breads,chrom,chromvars,varlist,reflist,fp);
		breads = clean_bam_buffer(readbuffer,breads,fout,1);
	}

	for (i=0;i<max_buffer_size;i++) { bam_destroy1(readbuffer[i]->b); free(readbuffer[i]); } 
	samclose(fp); if (fout != NULL) samclose(fout); // close the output bam file 
	//if (prevtid >=0 && reflist->sequences[prevtid] != NULL) free(reflist->sequences[prevtid]);
	//bam_destroy1(b);  // deletes memory associated with BAM record 'b'
	return 1;
}

//if (read->clipped == read->XC && ((read->mismatches < 1 && read->gaps < 2) || (read->gaps ==0 && read->mismatches < 2)))
// this function converts bam record to sam string
//samstring = bam_format1(bamheader_obj,b); fprintf(stdout,"sam format read %s \n",samstring);

// main function //
int main (int argc, char** argv)
{
	// use bamfile as integer index into argv and no need to do strcpy operations 
	
	char outfile[1024]; char bamfile[1024]; char variantfile[1024]; char fastafile[1024];
	strcpy(outfile,"None"); strcpy(bamfile,"None"); strcpy(variantfile,"None"); strcpy(fastafile,"None");
	char* sampleid = (char*)malloc(1024); sampleid[0] = '-'; sampleid[1] = '\0';
	int i=0,variants=0;

	for (i=1;i<argc;i+=2)
	{
		if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0)      strcpy(bamfile,argv[i+1]); 
		if (strcmp(argv[i],"--out") ==0 || strcmp(argv[i],"--output") ==0 )      strcpy(outfile,argv[i+1]);  // output bam/sam file
		if (strcmp(argv[i],"--ob") ==0)      OUTPUT_BAM = atoi(argv[i+1]);  // output all reads or only realigned reads 1/2
		else if (strcmp(argv[i],"--reffile") ==0 || strcmp(argv[i],"--ref") ==0)        strcpy(fastafile,argv[i+1]);
		else if (strcmp(argv[i],"--VCF") ==0 || strcmp(argv[i],"--vcf") ==0)    {     strcpy(variantfile,argv[i+1]);  }
		else if (strcmp(argv[i],"--indels") ==0)       PARSEINDELS = atoi(argv[i+1]);  // allow indels in hairs
		else if (strcmp(argv[i],"--pflag") ==0)      IFLAG  = atoi(argv[i+1]);  // allow indels in hairs
		else if (strcmp(argv[i],"--qvoffset") ==0)       QVoffset = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--novoalign")==0 || strcmp(argv[i],"--clip") ==0) 
		{
			CLIP_INSERTIONS = atoi(argv[i+1]);  
			if (CLIP_INSERTIONS ==1) fprintf(stderr,"insertions close to read ends will be soft-clipped before realignment\n");
		}
		else if (strcmp(argv[i],"--mmq") ==0)       MIN_MQ = atoi(argv[i+1]);
                else if (strcmp(argv[i],"--maxIS") ==0)       MAX_IS = atoi(argv[i+1]);
                else if (strcmp(argv[i],"--maxdelsize") ==0)       MAXDEL = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--split") ==0)       SPLITMAP = atoi(argv[i+1]);
                else if (strcmp(argv[i],"--singlereads") ==0)       SINGLE_READS = atoi(argv[i+1]);
                else if (strcmp(argv[i],"--minreads") ==0)       MIN_READS_FOR_INDEL = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--complexindels") ==0 || strcmp(argv[i],"--complex") ==0 || strcmp(argv[i],"--complexindel") ==0)
                {
                        COMPLEX_INDELS = atoi(argv[i+1]);
                        fprintf(stderr,"complex indels or indels close to each other will be output as single haplotype \n");
                }
	}
	if (strcmp(bamfile,"None") ==0 || strcmp(variantfile,"None")==0 || strcmp(fastafile,"None")==0) 
	{
		print_options_realigner(); return 1;
	}


	HASHTABLE ht; ht.htsize = 7919;  init_hashtable(&ht);
	variants = count_variants(variantfile); if (variants < 0) return -1; 
	VARIANT* varlist = (VARIANT*)malloc(sizeof(VARIANT)*variants);
	int chromosomes = read_variantfile(variantfile,varlist,&ht); 
	VARIANTS = variants;  
	CHROMVARS* chromvars  = (CHROMVARS*)malloc(sizeof(CHROMVARS)*chromosomes);
	build_intervalmap(chromvars,chromosomes,varlist,VARIANTS);

	// read reference fasta file for INDELS//
	REFLIST* reflist = NULL; reflist = init_reflist(fastafile,reflist); 
        FILE* fastafp = fopen(fastafile,"r"); //read_fasta(fastafile,reflist);

	int kmer = 9; int kmer1 = 8; chunksize = 65536; // 4^8
	//int tablesize = 1<<(2*kmer); double maxhitsforseed = (double)(chunksize+100000)*16/tablesize; MAX_KMER_HITS =(int)maxhitsforseed;
	int padding = 5*(MAX_IS+MAXDEL); // to account for insert size and max deletion size
	kmertable = (KMERTABLE*)malloc(sizeof(KMERTABLE)); kmertable->gap = 0; kmertable->kmer = kmer;
	kmertable1 = (KMERTABLE*)malloc(sizeof(KMERTABLE)); kmertable1->gap = 4; kmertable1->kmer = kmer1;

        if (SPLITMAP ==1)  // perform split read alignment 
        {
		fprintf(stderr,"initializing kmertable for split read mapping \n");
                init_kmertable(reflist,kmer,kmertable,chunksize+padding);
                if (HASHTABLES ==2) init_kmertable(reflist,kmer1,kmertable1,chunksize+padding);
        }

	realign_bamfile_sorted(bamfile,outfile,&ht,chromvars,varlist,reflist,fastafp);
	fclose(fastafp);
	return 0;
}


