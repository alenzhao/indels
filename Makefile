
CC=gcc -Wall -lm
CFLAGS=-c 
SAMTOOLS=../NGS_shared/samtools-0.1.18
NGSshared=../NGS_shared

testmapper: readfasta.o bamread.o kmertable.o indelfunctions.o splitreadmap.c blatindelmapper.c findhits.c
	$(CC) -I$(NGSshared) -I$(SAMTOOLS) -g -O2 bamread.o kmertable.o indelfunctions.o readfasta.o -o testmapper blatindelmapper.c -L$(SAMTOOLS) -lbam -lz -lm

realigner: bamread.o hashtable.o readvariant.o readfasta.o indel_realigner.c bam_realign.c realignment.c
	$(CC) -I$(NGSshared) -I$(SAMTOOLS) -g -O2 bamread.o hashtable.o readfasta.o readvariant.o -o INDELrealigner indel_realigner.c  -L$(SAMTOOLS) -lbam -lz

realigner1: bamread.o hashtable.o readvariant.o readfasta.o kmertable.o indelfunctions.o indel_realigner_new.c bam_realign.c realignment.c 
	$(CC) -I$(NGSshared) -I$(SAMTOOLS) -g -O2 bamread.o hashtable.o readfasta.o readvariant.o kmertable.o indelfunctions.o -o INDELrealigner1 indel_realigner_new.c  -L$(SAMTOOLS) -lbam -lz -lm


indelfunctions.o: indelfunctions.h indelfunctions.c bamread.h bamread.c printindels.c
	$(CC) -I$(NGSshared) -I$(SAMTOOLS) -c indelfunctions.c

readvariant.o: readvariant.c readvariant.h hashtable.h hashtable.c
	$(CC) -I$(NGSshared) -c readvariant.c 

bamread.o:	bamread.h bamread.c ../NGS_shared/readfasta.h ../NGS_shared/readfasta.c
	$(CC) -I$(NGSshared) -I$(SAMTOOLS) -c bamread.c

hashtable.o: ../NGS_shared/hashtable.h ../NGS_shared/hashtable.c
	$(CC) -c ../NGS_shared/hashtable.c

readfasta.o: ../NGS_shared/readfasta.c ../NGS_shared/readfasta.h
	$(CC) -c ../NGS_shared/readfasta.c

kmertable.o:    kmertable.h kmertable.c ../NGS_shared/readfasta.h ../NGS_shared/readfasta.c
	$(CC) -I$(NGSshared) -c kmertable.c

clean:
	rm -f bamread.o readfasta.o readvariant.o hashtable.o indelfunctions.o kmertable.o INDELrealigner  testmapper INDELrealigner1 
