
void print_options_realigner()
{
	fprintf(stderr,"\n PROGRAM TO realign reads in coordinate sorted BAM files \n\n");
	fprintf(stderr,"./INDELrealigner [options] --bam reads.sorted.bam --VCF variants.VCF --out realigned.bam  > output.log \n\n");
	fprintf(stderr,"=============== PROGRAM OPTIONS ======================================== \n\n");
	//fprintf(stderr,"--qvoffset : quality value offset, 33/64 depending on how quality values were encoded, default is 33 \n");
	//fprintf(stderr,"--mbq  : minimum base quality to consider a base for haplotype fragment, default 13\n");
	//fprintf(stderr,"--mmq : minimum read mapping quality to consider a read for phasing, default 20\n");
	fprintf(stderr,"--bam : input sorted bam file\n");
	fprintf(stderr,"--out : realigned bam file (output)\n");
	fprintf(stderr,"--ref : reference sequence file (in fasta format), should be indexed using samtools\n");
	fprintf(stderr,"--VCF : file with candidate indels in VCF format\n");
        fprintf(stderr,"--clip : 0/1, soft-clip insertions close to ends of reads before realignment, default 0\n");

	//fprintf(stderr,"--PEonly 0/1 : do not use single end reads, default is 0 (use all reads)\n");
	//fprintf(stderr,"--indels 0/1 : extract reads spanning INDELS, default is 0, variants need to specified in VCF format to use this option\n\n");
}

void print_options_indelmapper()
{
	fprintf(stderr,"\n###########  PROGRAM TO identify candidate indels from coordinate sorted SAM/BAM files  #############\n");
	fprintf(stderr,"\nIndels are identified using gapped alignments in the original bam file and using split-read alignment of unmapped/poorlymapped mates of paired-end reads (assumed to be Illumina paired-end reads in FR orientation). This program is designed for reads in the 30-150 base pair range \n\n");
	fprintf(stderr,"./indelcaller [options] --bam reads.sorted.bam --ref reference.fasta --out indels.VCF > candidateindels.out \n\n");
	fprintf(stderr,"[options] \n\n");
	//fprintf(stderr,"--qvoffset   : quality value offset for base qualities, 33/64 depending on how quality values were encoded, default is 33 \n");
	fprintf(stderr,"--mmq        : minimum read mapping quality to consider a read or its mate for indel detection, default 20\n");
	fprintf(stderr,"--maxIS      : maximum insert size for a paired-end read (used for searching neighborhood for mates), default 500\n");
	fprintf(stderr,"--maxdelsize : maximum size allowed for a deletion, default 5000\n");
	fprintf(stderr,"--ref  <FILE>      : reference sequence file (in fasta format), should be indexed using samtools\n");
	fprintf(stderr,"--pflag  <0/1>      : verbosity of output, default=0, 1 -> output alignments\n");
	fprintf(stderr,"--complexindels <0/1/2> : output complex indels/block substitutions as single indels, default = 0 \n"); 
	// need to output original cigar as multiple variants in VCF file
	fprintf(stderr,"--minreads <INT> : minimum number of reads supporting short indels (<5 bp) to output, default 3 \n");
	fprintf(stderr,"--singlereads <0/1> : use single-end reads for finding candidate indels, default 0\n");
	fprintf(stderr,"--split  <0/1>   : perform split-read mapping to find additional indels, default is 0\n\n");
	fprintf(stderr,"sample usage: ./testmapper --bam DATA/NA18507.300bp.sorted.bam --ref DATA/chr20.fa --out DATA/NA18507.testmapper.vcf > out \n\n");
}

