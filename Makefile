
CC=gcc -Wall -lm
CFLAGS=-c 
SAMTOOLS=../samtools

testmapper: readfasta.o bamread.o kmertable.o indelfunctions.o splitreadmap.c blatindelmapper.c findhits.c
	$(CC) -I$(SAMTOOLS) -g -O2 bamread.o kmertable.o indelfunctions.o readfasta.o -o testmapper blatindelmapper.c -L$(SAMTOOLS) -lbam -lz -lm

realigner: bamread.o hashtable.o readvariant.o readfasta.o indel_realigner.c bam_realign.c realignment.c
	$(CC) -I$(SAMTOOLS) -g -O2 bamread.o hashtable.o readfasta.o readvariant.o -o INDELrealigner indel_realigner.c  -L$(SAMTOOLS) -lbam -lz

realigner1: bamread.o hashtable.o readvariant.o readfasta.o kmertable.o indelfunctions.o indel_realigner_new.c bam_realign.c realignment.c 
	$(CC) -I$(SAMTOOLS) -g -O2 bamread.o hashtable.o readfasta.o readvariant.o kmertable.o indelfunctions.o -o INDELrealigner1 indel_realigner_new.c  -L$(SAMTOOLS) -lbam -lz -lm

indelfunctions.o:	indelfunctions.h indelfunctions.c ../readfasta.h ../readfasta.c bamread.h bamread.c printindels.c
	$(CC) -I$(SAMTOOLS) -c indelfunctions.c

readvariant.o: readvariant.c readvariant.h ../hashtable.h ../hashtable.c
	$(CC) -c readvariant.c 
bamread.o:	bamread.h bamread.c ../readfasta.h ../readfasta.c
	$(CC) -I$(SAMTOOLS) -c bamread.c
hashtable.o: ../hashtable.h ../hashtable.c
	$(CC) -c ../hashtable.c
readfasta.o: ../readfasta.c ../readfasta.h
	$(CC) -c ../readfasta.c
kmertable.o:    kmertable.h kmertable.c ../readfasta.h ../readfasta.c
	$(CC) -c kmertable.c

clean:
	rm -f bamread.o readfasta.o readvariant.o hashtable.o indelfunctions.o kmertable.o 
