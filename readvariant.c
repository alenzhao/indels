
#include "readvariant.h"
#include "../readfasta.h"

// count the # of variants in VCF file to allocate space for VCF variant array
int count_variants(char* vcffile)
{
	FILE* fp = fopen(vcffile,"r"); 
	if (fp == NULL) { fprintf(stderr,"could not open file %s\n\n",vcffile); return -1; }
	int variants =0; char buffer[10000];

	while (fgets(buffer,10000,fp))
	{
		if (buffer[0] != '#') variants++; // this should work for non-VCF files as well.
	}
	fclose(fp);
	fprintf(stderr,"VCF file %s has %d variants \n",vcffile,variants);
	return variants;
}

// in VCF format all data lines are tab-delimited but we still allow spaces (why ??) 
// we do not check the VCF file for consistency with format, assume that it is in correct format 
int parse_variant(VARIANT* variant, char* buffer)
{
	int i=0,j=0,s=0,e=0; int flag =0;
	char* tempstring;
	variant->depth =0; variant->A1 =0; variant->A2 = 0; variant->H1 = 0; variant->H2 =0;

	while (buffer[i] == ' ' || buffer[i] == '\t') i++; s = i; while (buffer[i] != ' ' && buffer[i] != '\t') i++; e= i;
	variant->chrom = (char*)malloc(e-s+1); for (j=s;j<e;j++) variant->chrom[j-s] = buffer[j]; variant->chrom[j-s] = '\0';
	while (buffer[i] == ' ' || buffer[i] == '\t') i++; s = i; while (buffer[i] != ' ' && buffer[i] != '\t') i++; e= i;
	tempstring = (char*)malloc(e-s+1); for (j=s;j<e;j++) tempstring[j-s] = buffer[j]; tempstring[j-s] = '\0';
	variant->position = atoi(tempstring); free(tempstring);

	while (buffer[i] == ' ' || buffer[i] == '\t') i++; s = i; while (buffer[i] != ' ' && buffer[i] != '\t') i++; e= i; // varid
	while (buffer[i] == ' ' || buffer[i] == '\t') i++; s = i; while (buffer[i] != ' ' && buffer[i] != '\t') i++; e = i;
	variant->RA = (char*)malloc(e-s+1); for (j=s;j<e;j++) variant->RA[j-s] = buffer[j]; variant->RA[j-s] = '\0';
	while (buffer[i] == ' ' || buffer[i] == '\t') i++; s = i; while (buffer[i] != ' ' && buffer[i] != '\t') i++; e= i;
	variant->AA = (char*)malloc(e-s+1); for (j=s;j<e;j++) variant->AA[j-s] = buffer[j]; variant->AA[j-s] = '\0';

	while (buffer[i] == ' ' || buffer[i] == '\t') i++; s = i; while (buffer[i] != ' ' && buffer[i] != '\t') i++; e= i;
	while (buffer[i] == ' ' || buffer[i] == '\t') i++; s = i; while (buffer[i] != ' ' && buffer[i] != '\t') i++; e= i;
	while (buffer[i] == ' ' || buffer[i] == '\t') i++; s = i; while (buffer[i] != ' ' && buffer[i] != '\t') i++; e= i;
	while (buffer[i] == ' ' || buffer[i] == '\t') i++; s = i; while (buffer[i] != ' ' && buffer[i] != '\t') i++; e= i;

	variant->type = strlen(variant->AA)-strlen(variant->RA); 
	if (strlen(variant->AA) ==1 || strlen(variant->RA) ==1) variant->simple = '1';
	else variant->simple = '0';
	variant->rightshift =-1; // set to -1
	//if (variant->type != 0) variant->position++; // add one to position for indels 

}

// change this to VCF file now 
int read_variantfile(char* vcffile,VARIANT* varlist,HASHTABLE* ht)  
{
	FILE* fp = fopen(vcffile,"r");
	char buffer[10000];
	int i=0;
	//	char allele1[256]; char allele2[256]; char genotype[256]; int quality; 
	char prevchrom[256];  strcpy(prevchrom,"----"); 
	int chromosomes = 0; //int blocks=0;

	while (fgets(buffer,10000,fp))
	{
		if (buffer[0] == '#') continue;
		else
		{
			//fprintf(stdout,"buffer %s \n",buffer);
			parse_variant(&varlist[i],buffer); 
			//fprintf(stdout,"%s %d %s %s %s %s\n",varlist[i].chrom,varlist[i].position,varlist[i].RA,varlist[i].AA,varlist[i].genotype,prevchrom); 
			if (strcmp(varlist[i].chrom,prevchrom) !=0)
			{
				//	fprintf(stderr,"chromosomes %d %d\n",chromosomes,i);
				// insert chromname into hashtable 
				insert_keyvalue(ht,varlist[i].chrom,strlen(varlist[i].chrom),chromosomes);
				strcpy(prevchrom,varlist[i].chrom); chromosomes++;
			}
			i++;
		}
	}
	fclose(fp); //chromosomes--;
	fprintf(stderr,"vcffile %s chromosomes %d %d\n",vcffile,chromosomes,i);
	return chromosomes;	

}

// build a physical map that maps  intervals on chromosomes to the first variant that precedes the start of that interval
void build_intervalmap(CHROMVARS* chromvars,int chromosomes,VARIANT* varlist,int snps)
{
	int i=0,j=0,k=0,blocks=0;
	chromvars[j].first = 0;  j=0;
	for (i=0;i<snps-1;i++) 
	{
		if (strcmp(varlist[i].chrom,varlist[i+1].chrom) !=0)
		{
			chromvars[j].last = i; chromvars[j].variants = chromvars[j].last-chromvars[j].first+1;
			//fprintf(stderr,"new chrom %d %d %s %s\n",j,chromvars[j].variants,varlist[i].chrom,varlist[i+1].chrom);
			j++; chromvars[j].first = i+1; 		
		}
	}
	chromvars[j].last = i;
	//	int** intervalmap; // map 1000bp of chromosome to first snp in that region indexed by snp_array 
	// first SNP to the right of the given base position including that position
	for (j=0;j<chromosomes;j++)
	{
		blocks = (int)(varlist[chromvars[j].last].position/BSIZE)+2; chromvars[j].blocks = blocks;
		//	fprintf(stderr,"chromosomes %d blocks %d \n",j,blocks);
		chromvars[j].intervalmap = (int*)malloc(sizeof(int)*blocks);
		for (i=0;i<blocks;i++) chromvars[j].intervalmap[i] = -1;
		//fprintf(stderr,"blocks for chrom %d: %d \n",j,blocks);
		k = chromvars[j].first;
		for (i=0;i<blocks;i++)
		{
			while (varlist[k].position <= BSIZE*i && k < chromvars[j].last) k++; 
			if (k == chromvars[j].last) break;
			if (varlist[k].position > BSIZE*i && chromvars[j].intervalmap[i] == -1) chromvars[j].intervalmap[i] = k; 
			//		if (chromvars[j].intervalmap[i] != -1) printf("FSNPtoright chrom %d block %d: %d %d \n",j,BSIZE*i,chromvars[j].intervalmap[i],varlist[chromvars[j].intervalmap[i]].position);
			//			else printf("FSNPtoright chrom %d block %d: %d \n",j,BSIZE*i,intervalmap[j][i]);
		}
	}
}

void print_indelhaplotypes(VARIANT* varlist,int ss,REFLIST* reflist,int tid)
{
	int i=0;
	for (i=0;i<varlist[ss].A1;i++) fprintf(stdout,"%c",varlist[ss].RA[i]); 
	for (i=0;i<varlist[ss].rightshift+1;i++) fprintf(stdout,"%c",reflist->sequences[tid][i+varlist[ss].position+varlist[ss].A1-1]);

	fprintf(stdout,"/");
	for (i=0;i<varlist[ss].A2;i++) fprintf(stdout,"%c",varlist[ss].AA[i]); 
	for (i=0;i<varlist[ss].rightshift+1;i++) fprintf(stdout,"%c",reflist->sequences[tid][i+varlist[ss].position+varlist[ss].A1-1]);
	fprintf(stdout," ");
}

// this will only work for pure insertions and deletions, not for block substitutions, needs to have tid variable set
int calculate_rightshift(VARIANT* varlist,int ss,REFLIST* reflist,int tid)
{
        int i=0,j=0;  int a1=0,a2=0; int shift = 0;
        a1 = strlen(varlist[ss].RA); a2 = strlen(varlist[ss].AA); varlist[ss].A1 = a1; varlist[ss].A2 = a2; 
        if (a1 > a2 && a2 ==1) // deletion
        {
                i= varlist[ss].position; // first base of deletion assuming position is +1 and not previous base 
                j = varlist[ss].position + a1-a2;
                while (i < reflist->lengths[tid] && j < reflist->lengths[tid] && reflist->sequences[tid][i] == reflist->sequences[tid][j])
                { 
                        i++; j++; shift++;
                } 
		varlist[ss].rightshift = shift; 
        }
        else if (a1 ==1 && a2 > a1) // insertion event 
        {
                i = 1; j = varlist[ss].position;
                while (j+1 < reflist->lengths[tid] && varlist[ss].AA[i] == reflist->sequences[tid][j+1] && i < a2) { i++; j++; shift++; } 
                if (i == a2) // covered the full length of the inserted bases
                {
                        i = varlist[ss].position; 
                        while (i < reflist->lengths[tid] && j < reflist->lengths[tid] && reflist->sequences[tid][i] == reflist->sequences[tid][j]) { i++; j++; shift++; } 
                }
		varlist[ss].rightshift = shift; 
        }
	else varlist[ss].rightshift =0; // complex indel
}
