#! /usr/bin/env python
import sys, os, glob, string, subprocess,time, math, re, compiler
from array import *

###### PROGRAM TO EXTRACT LIST OF INDELS FROM A SAM FILE (or BAM file using stdin) #########3
###### Author: Vikas Bansal, updated aug 7 2011 for stampy

###### options for this program 1. samfile 2. refseq.fasta 3. minreads_to_report_indel 
###### 4. quality of read to extract indels (BFAST specific) 

######################################################################

### unused function 
def leftshift(bases,lb,type,refsequence,start,l1,l2,read):
	if type == 'D':
		k1 = start+l2-2; k2 = start+l2-1 + lb-1;  shift = 0;
		while refsequence[0][k1] == refsequence[k2] and k1 > 0: k1 -=1; k2 -=1; shift  +=1; 
		return shift; 
	elif type == 'I':
		k1 = start+l2-2; k2 = l1+lb-1; shift = 0;
		while refsequence[0][k1] == read[k2] and k1 > 0 and k2 > 0: k1 -=1; k2 -=1; shift +=1; 
		return shift; 

def read_fasta(seqfile):
	sequences = {}; reflist = []; lines =0; strlist = []; sl =0;
	File = open(seqfile,'r');
	if not File: print >>sys.stderr, 'reference sequence file',seqfile,'not found'; sys.exit();
	for line in File:
		if line[0] == '>': 
			if lines > 0:
				print >>sys.stderr, 'added seq lines',refname,lines;
				sequences[refname][0] += ''.join(strlist); lines =0; 
				while sl >0: strlist.pop(); sl -=1;
				#for i in range(len(strlist)): strlist.pop();
			refname = line.strip().lstrip('>').split()[0];
			sequences[refname] = ['',0];    # refsequence and basecalls 
			reflist.append(refname);
		else: strlist.append(line.strip().upper()); lines +=1; sl +=1;
	sequences[refname][0] += ''.join(strlist); 
	print >>sys.stderr, 'added seq',refname,lines;
	File.close();
	return [sequences,reflist];
	
def samtopileup(readfile,seqfile,minreads,bedfile):
	[sequences,reflist] = read_fasta(seqfile); 
 	overlapindex = {}; indexpresent =0;
	if bedfile != "": 
		File = open(bedfile,'r'); indexpresent = 1;
		for line in File: 
			target = line.strip().split();
			for p in xrange(int(target[1]),int(target[2])):  overlapindex[(target[0],p)] = 1; 
		File.close();
#	for p in overlapindex.keys(): print p;
			
	indeltable = {}; readidtable = {};

	bamfilelist = glob.glob(readfile);
	for bamfile in bamfilelist: print >>sys.stderr, bamfile; 

	for bamfile in bamfilelist: 
		print >>sys.stderr, 'reading',bamfile;
		samtoolsin = subprocess.Popen(["/home/vbansal-scripps/bin/samtools-0.1.18/samtools","view",bamfile],stdout=subprocess.PIPE,bufsize=1);
		File = samtoolsin.stdout;
		
		#if readfile == "sys.stdin" or readfile =='stdin' or readfile == '-': File = sys.stdin;
		#else: File=open(readfile,'r');  
		#if not File: print >>sys.stderr, 'samfile',readfile,'not found'; sys.exit();

		reads =0;
		for s in File:
			if s[0] == '@': continue;
			reads +=1; 
			if reads%1000000 ==0: print >>sys.stderr, 'processed',reads;

			line = s.strip().split();
			ll = len(line); 
			if ll < 11: continue;
			if 'D' in line[5] or 'I' in line[5]: indelread = 1; 
			else: continue; 
			stampy =0;
			if 'STAMPY' in line[12]: stampy = 1;
			flag = int(line[1]); locus = line[2]; start = int(line[3]); mq = int(line[4]); read = line[9]; #quality = line[10]; 

			if flag & 4 == 4: continue;
			if (mq < 30 and stampy ==1) or (mq < 20 and stampy ==0): continue;

			XM = 0;
			for i in xrange(11,ll): 
				if 'XM:i:' in line[i]: XM = int(line[i].split(':')[2]);
			if XM > 4: continue;
			
			if flag & 16 ==16: strand = '-';
			else: strand = '+'; # + strand 


			cigarstring = line[5]; 
			cigarlist = []; prev = 0;  c= 0;
			for i in xrange(len(cigarstring)):
				b = ord(cigarstring[i]);
				if b >= 68 and b < 84: cigarlist.append(cigarstring[i]);  cigarlist.append(int(cigarstring[prev:i])); prev = i+1; c+=2;
				#if b < 48 or b > 57: ml = int(cigarstring[prev:i]); cigarlist.append(cigarstring[i]);  cigarlist.append(ml); prev = i+1; c+=2;
			
			l1 = 0; l2 =0;
			for k in xrange(0,c,2):
				if cigarlist[k] == 'M': l1 += cigarlist[k+1]; l2 += cigarlist[k+1];
				elif cigarlist[k] == 'S':		l1 += cigarlist[k+1];
				elif cigarlist[k] == 'D':
					try: 
						refsequence = sequences[locus];  # index into refnames table
						k1 = start+l2-2; k2 = start+l2-1 + cigarlist[k+1]-1;  shift = 0;
						while refsequence[0][k1] == refsequence[0][k2] and k1 > 0: k1 -=1; k2 -=1; shift  +=1; 
						bases = refsequence[0][start+l2-1-shift:start+l2-1+cigarlist[k+1]-shift];
						start -= shift;
						try: 
							indeltable[(locus,start+l2-1,refsequence[0][start+l2-2]+bases,refsequence[0][start+l2-2])].append(strand + cigarstring);
							readidtable[(locus,start+l2-1,refsequence[0][start+l2-2]+bases,refsequence[0][start+l2-2])].append(line[0].split('.')[0]);
						except KeyError: 
							if strand =='+': s1 = 1; s2 =0;
							else: s1 = 0; s2 =1;
							#indeltable[locus + ','+str(start+l2)+',' + 'D:'+bases+':-'] =[s1,s2];
							indeltable[(locus,start+l2-1,refsequence[0][start+l2-2]+bases,refsequence[0][start+l2-2])] = [strand+cigarstring];
							readidtable[(locus,start+l2-1,refsequence[0][start+l2-2]+bases,refsequence[0][start+l2-2])] = [line[0].split('.')[0]];
					except KeyError: pass;
					l2 += cigarlist[k+1];
				elif cigarlist[k] == 'I' and k != 0 and k!= c-2: 

					### important to push inserted bases left as possible, as stampy doesn't correct job 
					k1 = l1-1; k2 = l1+cigarlist[k+1]-1; shift = 0;
					while read[k1] == read[k2] and k1 > 0 and k2 > 0: k1 -=1; k2 -=1; shift +=1; 
					if l1- shift > 0: bases = read[l1-shift:l1+cigarlist[k+1]-shift]; start -= shift; 
					else: bases = read[l1:l1+cigarlist[k+1]];

					try: 
						refsequence = sequences[locus];  # index into refnames table
						indeltable[(locus,start+l2-1,refsequence[0][start+l2-2],refsequence[0][start+l2-2]+bases)].append(strand+cigarstring);
						readidtable[(locus,start+l2-1,refsequence[0][start+l2-2],refsequence[0][start+l2-2]+bases)].append(line[0].split('.')[0]);
					except KeyError: 
						if strand =='+': s1 = 1; s2 =0;
						else: s1 = 0; s2 =1;
						indeltable[(locus,start+l2-1,refsequence[0][start+l2-2],refsequence[0][start+l2-2]+bases)] = [strand+cigarstring];
						readidtable[(locus,start+l2-1,refsequence[0][start+l2-2],refsequence[0][start+l2-2]+bases)] = [line[0].split('.')[0]];
					l1 += cigarlist[k+1];
				else: break;


		File.close();


	### cluster long insertions of same length that differ by one base  August 26 2011
	print '##fileformat=VCFv4.0'; print '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'; 
	indellist = [];
	for indel in indeltable.iterkeys(): indellist.append([indel[0],indel[1],indel]); 
	indellist.sort();

	### cluster indels at same position, TACAC -> TAC,T 

	for i in xrange(len(indellist)):
		### filter indels that don't overlap target bases from bed file
		indel = indellist[i][2]; 
		if indexpresent ==1: 
			try: 
				#print indel[0],indel[1]+1;
				OV = overlapindex[(indel[0],indel[1]+1)]; 
			except KeyError: continue; 

		cigars = indeltable[indel]; readids = readidtable[indel];
		if len(cigars) >= minreads: 
			cigars.sort(); readids.sort(); 
			unique = 1; samples = 1; readsf =0; readsr=0;
			for i in xrange(len(cigars)):
				if cigars[i][0] == '+': readsf +=1;
				else: readsr +=1; 

			for i in xrange(len(cigars)-1):
				if cigars[i] != cigars[i+1]: unique +=1;
				#if readids[i] != readids[i+1]: samples +=1; 
			if unique >= minreads: 	
		#if count[0] + count[1] >= minreads or (count[0] > 0 and count[1] > 0 and count[0]+ count[1] >= minreads-1): 
				sys.stdout.write('%s\t%d\t.\t%s\t%s\t%d\tPASS\t' %(indel[0],indel[1],indel[2],indel[3],len(cigars)))
				sys.stdout.write('COUNTS=%d,%d;UNIQUE=%d,SAMPLES=%d\n' %(readsf,readsr,unique,samples));
				#for i in xrange(min(len(cigars),5)): print cigars[i],
				#print;


if len(sys.argv) > 4: 
	bedfile = sys.argv[4];
	samtopileup(sys.argv[1],sys.argv[2],int(sys.argv[3]),bedfile);
elif len(sys.argv) > 3: samtopileup(sys.argv[1],sys.argv[2],int(sys.argv[3]),"");
else: 
	print 'python parse_samfile_indels_vcf.py BAMfile refsequence.fa minreads_to_report_indel'; sys.exit();

### also provide option of bedfile to filter out indels outside target region...

### samtools view GATK-realignment-200samples/SRR028812.chr20.forGATK.bam | python parse_samfile_indels_vcf.py stdin /home/vbansal-scripps/GENOMES/1000genomes-huref/single-chromosomes/chr20.fa 3 /home/vbansal-scripps/GENOMES/Nimblegen.exome/Nimblegen.exome.bed.ncbi36.chr20  > 28812.indels.vcf


