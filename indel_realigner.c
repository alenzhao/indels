
// original code that just does indel realignment using bam file and VCF file 
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
#include "readvariant.h"
#include "print_options.c"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MINQ = 13; // minimum base quality
int MIN_MQ = 20; // minimum read mapping quality
int MAX_IS =  1000; // maximum insert size
int MIN_IS =  0; // maximum insert size
int PEONLY = 0; // if this is set to 1, reads for which only one end is mapped are not considered for hairs 
int BSIZE = 500; int IFLAG = 0; int VARIANTS = 0; int PARSEINDELS =0; int SINGLEREADS =0;
//int QVoffset = 33; declared in samread.h
int PFLAG = 1;
int OUTPUT_BAM =1; 
int MAX_REALIGNMENTS = 100;
int CLIP_INSERTIONS = 0; // parameter for novoalign fullNW alignments, remove excessive 'I' cigars at ends of reads 
int STRINGENCY = 0; // 0/1/2 to control how much realignment is done close to ends of reads...

#include "bam_realign.c"

int realign_bamfile_sorted(char* bamfile,char* outfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist,FILE* fastafp)
{
	fprintf(stderr,"reading sorted bamfile %s \n",bamfile);
	int i=0; int chrom=0; int reads=0;int realigned=0,clipped_ins=0;   int prevchrom=-1; int prevtid = -1; 
	struct alignedread* read = (struct alignedread*)malloc(sizeof(struct alignedread));
        int* fcigarlist = (int*)malloc(sizeof(int)*4096);

	// hold information about the realignment candidates
	struct REALIGNMENT* realignments = calloc(MAX_REALIGNMENTS,sizeof(struct REALIGNMENT)); 
	for (i=0;i<MAX_REALIGNMENTS;i++) realignments[i].cigarlist = calloc(256,sizeof(uint32_t)); 
	
	samfile_t *fp; 
	if ((fp = samopen(bamfile, "rb", 0)) == 0) 
	{ 
		fprintf(stderr, "Fail to open BAM file %s\n", bamfile); return -1; 
	}
	bam1_t *b = bam_init1();  // bam record 
	//bam_header_t* bamheader_obj = bam_header_read(fp); // read header of bam file  // char* samstring;
	samfile_t* fout = NULL;
	if (strcmp(outfile,"None") != 0) 
	{
		fout = samopen(outfile,"wb",fp->header); 
		if (fout ==0) fprintf(stderr,"failed to open output bam file %s \n",outfile); 
	}

	while (samread(fp, b) >= 0)
	{
		// what about unmapped read pairs read->tid = -1 at end of sorted bam file
		// should PCR duplicate reads be realigned as well -> since if we do not their cigar will differ from the duplicate copy 07/18/13
		if (b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP))
		{
			if (OUTPUT_BAM==1 && fout != NULL) samwrite(fout,b); 
			continue;
		}
		fetch_func(b, fp,read);
		if (read->tid != prevtid && read->tid >= 0)
		{
			if (prevtid >=0) { free(reflist->sequences[prevtid]); reflist->sequences[prevtid] = NULL; } 
                        if (prevtid >=0) fprintf(stderr,"freeing memory from previous chrom %s.. reads %d ",reflist->names[prevtid],reads);
                        read_chromosome(reflist,read->tid,fastafp); chrom = getindex(ht,read->chrom); 
		}
		else chrom = prevchrom;
		// chrom is index in hashtable 'ht' used  for VCF index
		// need to discard reads that are marked as duplicates using flag //
		if (read->mquality >= MIN_MQ && chrom >= 0 && read->tid >=0) 
		{
			clipped_ins =0; if (CLIP_INSERTIONS ==1 || CLIP_INSERTIONS ==0) clipped_ins = modify_insertion_cigar(read); 
			//if (read->gaps > 1) combine_indelcigars(read,reflist); // for multiple indels per read
			parse_cigar(read,reflist,fcigarlist);
		
			realigned = -1;
			if (read->clipped > read->XC || read->mismatches >= 1 || read->gaps > 1 || read->mismatches + read->gaps >=2) 
			{
				realigned = extract_variants_read(read,chromvars,varlist,chrom,reflist,realignments);
			}
			//realigned = -1;
			if (realigned >0)
			{
				// edit cigarlist to add 'H' ops at start and end if needed 
				update_bam_record(b,realignments[realigned].newpos,realignments[realigned].cigarlist,realignments[realigned].cigs,1); 
				if (OUTPUT_BAM==2 && fout != NULL) samwrite(fout,b); 
			}
			else if (clipped_ins > 0) update_bam_record(b,read->position,read->cigarlist,read->cigs,0); 

		}
		if (OUTPUT_BAM==1 && fout != NULL) samwrite(fout,b); 
		if (++reads%1000000 ==0) fprintf(stderr,"processed %d reads\n",reads);
		prevchrom = chrom; prevtid = read->tid; free_readmemory(read);
	}
	bam_destroy1(b); samclose(fp); 
	if (fout != NULL) samclose(fout); // close the output bam file 
	for (i=0;i<MAX_REALIGNMENTS;i++) free(realignments[i].cigarlist); free(realignments); 
	free(fcigarlist);
	//if (prevtid >=0 && reflist->sequences[prevtid] != NULL) free(reflist->sequences[prevtid]);
	return 1;
}

//if (read->clipped == read->XC && ((read->mismatches < 1 && read->gaps < 2) || (read->gaps ==0 && read->mismatches < 2)))
// this function converts bam record to sam string
//samstring = bam_format1(bamheader_obj,b); fprintf(stdout,"sam format read %s \n",samstring);

// main function //
int main (int argc, char** argv)
{
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
			if (CLIP_INSERTIONS ==1) fprintf(stderr,"insertions close to ends of reads will be converted to soft-clips before realignment\n");
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
	realign_bamfile_sorted(bamfile,outfile,&ht,chromvars,varlist,reflist,fastafp);
	fclose(fastafp);
	return 0;
}


