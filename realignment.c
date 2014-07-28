/* code for generating new cigar using indel: algorithm description 
1. find all potential indels that overlap the aligned/partially aligned read (variant index)
2. for each indel, evaluate if adding the indel into the existing cigar of the read reduces the # of mismatches or extends into the clipped part of the read 
3. the above step is done for both directions where indel is inserted L->R and R->L 
*/

// need to be modified to handle complex indels or block substitutions, done on jan 8 2013, working okay 
// can we use the code to change cigar with multiple indels into single complex indel ?? require handling insertion just before insertion of 'evaluated indel'

// python parse_samfile_indels_vcf.py /projects/stsi3/PROJECTS/T2D-pooledseq/T2D-pooledseq-mergedbams/G1P1.merged.picard.bam /projects/stsi/vbansal/T2D-pooledseq-july2012/target/ncbi37.fa 3 > G1P1.indels.vcf
// ./INDELrealigner --bam /projects/stsi3/PROJECTS/T2D-pooledseq/T2D-pooledseq-mergedbams/G1P1.merged.picard.bam --ref /projects/stsi/vbansal/T2D-pooledseq-july2012/target/ncbi37.fa --VCF G1P1.indels.vcf > G1P1.indels.log
// should count mismatches after event separately... to compare ...

// carryover is the number of inserted bases that are not yet accounted for in an insertion event..
// MM[0] is mismatches before the indel event is added to original cigar: shared between the two cigars 
// (read->mismatches - MM[0])   vs MM[1] is the new comparison (how many differences after the indel event is added) 

//char INT_CIGAROP[] = {'M','I','D','N','S','H','P','E','X'};

int realignread_LR(struct alignedread* read,REFLIST* reflist,VARIANT* varlist,int var,struct REALIGNMENT* rl)
{
	rl->cigs =0; //initilze
	// extend partially mapped reads to calculate total score without clips, to avoid making spurios indel realignments...
	int i=0,j=0,l1=0,l2=0,t=0,matches =0, carryover = 0,dl=0,il=0;
	int l11=0,l21=0; int added = -1; int MM[3] = {0,0,0};  // pair of new mismatches,old mismatches before indel, new mismatches after indel event
	int op=0,ol=0,p = 0; int indeladded =0; int b =0;
	//int delta =0; if (varlist[var].type < 0) delta = -1*varlist[var].type; 
	//printf("position %d start %d end %d \n",read->position+l2,varlist[var].position,read->span+read->position);
	int lastposition = read->readlength; // lastposition in read, if read is clipped due to poor quality bases, lastposition < readlength
	op = read->cigarlist[read->cigs-1]&0xf;  ol = read->cigarlist[read->cigs-1]>>4;
	if ( (read->flag &16) ==0 && op == BAM_CSOFT_CLIP) lastposition = read->readlength-read->XC; 
	
	// check read going from left to right 
	for (i=0;i<read->cigs;i++)
	{
	        op = read->cigarlist[i]&0xf; ol = read->cigarlist[i]>>4; p = ol;
		// only move into soft clip region if previous cigar is not 'M' extra condition added jan 2 2013 
		if (op == BAM_CSOFT_CLIP && i  > 0 && rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0 && (rl->cigarlist[rl->cigs-1]>>4) >=10 && indeladded > 0)  
		{
			//fprintf(stdout,"ignoring soft clip region \t");
			rl->cigarlist[rl->cigs++] = read->cigarlist[i]; l1 += ol; l11 += ol; continue;
		}
		else if (op == BAM_CMATCH || (op == BAM_CSOFT_CLIP && i > 0) || (op == BAM_CINS && indeladded > 0) || (op == BAM_CINS  && read->position + l2 == varlist[var].position && varlist[var].type > 0 && indeladded ==0) )
		{
			//if (op == BAM_CINS && indeladded > 0) MM[1]++;  
			for (t=0;t<ol;t++)
			{
				if (read->position + l2 == varlist[var].position && indeladded ==0)  // indel insertion match in 'M' stretch
				{
					// here 't > 0' is necessary otherwise p=0, indel cannot happen | maybe not check this 06/14/13
					if (varlist[var].type < 0 && varlist[var].simple == '1' && (rl->cigs > 0 || t > 0))
					// position in read matches the known position of deletion event 
					{
						p = ol-t; 
						if (t > 0 && rl->cigs ==0) rl->cigarlist[rl->cigs++] = t<<4; 
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = ((-1*varlist[var].type)<<4) + 2; 
						l2 += -1*varlist[var].type; added = 0; carryover = 0; indeladded = 1;
						printf("LR_del->%d %d ",l1,p);
					}
					else if (varlist[var].type > 0 && varlist[var].simple == '1' && (rl->cigs > 0 || t > 0))
					{ 
						p = ol-t - varlist[var].type; 
						if (p < 0) carryover = -1*p; else carryover = 0;
						if (t > 0 && rl->cigs== 0) rl->cigarlist[rl->cigs++] = t<<4;
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = (varlist[var].type<<4) + 1;
						for (j=0;j<varlist[var].type && j < ol-t;j++)  // count mismatches of inserted bases to known insertion allele
						{
							if (read->sequence[l1+j]  != varlist[var].AA[j+1]) MM[2]++;  
						}
						printf("ins->%d %d %d %c carryover %d ",l1,p,t,INT_CIGAROP[op],carryover);
						l1 += varlist[var].type; t += varlist[var].type;
						// when t is incremented by more than is covered by this segment, then something needs to subtracted from next segment.....
						added = 0; indeladded =1;
					}
					else if (varlist[var].simple == '0' && (rl->cigs > 0 || t > 0)) // complex indel still not fully fixed 
					{ 
						dl = strlen(varlist[var].RA)-1; il = strlen(varlist[var].AA)-1;
						p = ol-t; p -= il; 
						if (p < 0) carryover = -1*p; else carryover = 0;
						if (t > 0 && rl->cigs== 0) rl->cigarlist[rl->cigs++] = t<<4;
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = (dl<<4) + 2; rl->cigarlist[rl->cigs++] = (il<<4) + 1;
						for (j=0;j<il && j < ol-t;j++)  // count mismatches of inserted bases to known insertion allele
						{
							if (read->sequence[l1+j]  != varlist[var].AA[j+1]) MM[2]++;  
						}
						printf("complexindel->%d %d %d %c carryover %d ",l1,p,t,INT_CIGAROP[op],carryover);
						l1 += il; t += il; l2 += dl;
						// when t is incremented by more than is covered by this segment, then something needs to subtracted from next segment.....
						added = 0; indeladded =1;
					}
				}

				if (t >=ol) break;
				if ( (op == BAM_CMATCH || op == BAM_CSOFT_CLIP || op == BAM_CINS) && carryover > 0) { carryover--; p--; continue; } 
				if (l1 >= lastposition) { break;} // important check
				if (added >= 0 && op == BAM_CSOFT_CLIP) added++; // variant has been added and now into clipped region 

				// if t >= ol because of insertion -> these won't be counted jan 1 2013 
				//if (read->sequence[l1]  != reflist->sequences[read->tid][read->position+l2] && indeladded ==0) MM[1]++;
				if (read->sequence[l11+t]  != reflist->sequences[read->tid][read->position+l21+t] && indeladded ==0 && op != BAM_CINS) MM[0]++;

				if (read->sequence[l1]  != reflist->sequences[read->tid][read->position+l2]) 
				{
					if (op == BAM_CSOFT_CLIP && added >= 5 && carryover <=0 && varlist[var].type < 0) 
					{ 
						t = ol; added--;   // we should break here rather than continue since cigar 
					} // get out
					else if (indeladded > 0 || op != BAM_CMATCH) MM[2]++; 
				}
				else matches++; 
				l1++; l2++; 
			}
			if (added < 0)  // no indel was inserted into newcigarlist
			{ 
				if (rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0) rl->cigarlist[rl->cigs-1] += ol<<4;
				else rl->cigarlist[rl->cigs++] = read->cigarlist[i]; 
			}
			else if (p > 0 && (op == BAM_CMATCH || op == BAM_CINS)) // left over 'M' after indel insertion  
			{
				if (rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0) rl->cigarlist[rl->cigs-1] += p<<4;
                                else rl->cigarlist[rl->cigs++] = p<<4; 
			}
			else if (carryover ==0 && op == BAM_CSOFT_CLIP) 
			{ 
				b = rl->cigarlist[rl->cigs-1]&0xf;
				if (added > 0 && b != 0) rl->cigarlist[rl->cigs++] = added<<4; 
				else if (added > 0) rl->cigarlist[rl->cigs-1] += added<<4; 
				if (p-added > 0) rl->cigarlist[rl->cigs++] = ((p-added)<<4) + 4; // remaining portion of soft clip
			}
			else if (carryover > 0 && op == BAM_CSOFT_CLIP) 
			{ 
				rl->cigarlist[rl->cigs++] = (ol-carryover)<<4; carryover =0;
			}
			if (op == BAM_CINS) l11 += ol; else { l11 += ol; l21 += ol;} 
		}
		else if (op == BAM_CDEL) 
		{ 
			if (indeladded > 0) { MM[1]++; } // already added the indel event we are evaluating, so ignore extra deletion event
			else if (read->position + l2 + ol < varlist[var].position ) // added 07/17/13, deletion is ignored if the new variant being evaluated falls within the interval deleted on reference 
			{
				if (carryover > 0) {  rl->cigarlist[rl->cigs-1] -= carryover<<4; carryover = 0; } // probably never invoked
				rl->cigarlist[rl->cigs++] = read->cigarlist[i]; 
				l2 += ol;  
			}
			l21 += ol; 
		}
		else if (op == BAM_CINS) // add filter on insertion being close to variant being evaluated...  07/17/13
		{ 
			rl->cigarlist[rl->cigs++] = read->cigarlist[i]; l1 += ol; l11 += ol; 
		}
		else if (op == BAM_CSOFT_CLIP) 
		{ 
			rl->cigarlist[rl->cigs++] = read->cigarlist[i]; l1 += ol; l11 += ol; 
		}
	}
	// if the previous cigarlength is less than carryover -> negative BUG 
	if (carryover > 0)  
	{
	       	op = rl->cigarlist[rl->cigs-1]&0xf; ol = rl->cigarlist[rl->cigs-1]>>4; 
 		if (op != BAM_CDEL && ol > carryover)  rl->cigarlist[rl->cigs-1] -= carryover<<4; 
		else return 0;
	}

	// MM[0] is number of mismatches prior to indel event, MM[1] -> gaps
	// MM[2] is number of mismatches introduced after indel event is added to new cigar
	//fprintf(stdout,"newciglength %d \n",rl->cigs);
	if ( (added ==0 && MM[2] < read->mismatches-MM[0]+MM[1]-1 && MM[2] <2) || (added == 0 && MM[2] < read->mismatches-MM[0] && MM[2] < 3) || (added >= 2 && MM[2] < 1) || (added >= 5 && MM[2] <2)) 
	{
		rl->newpos = read->position; rl->added = added; rl->mismatches = MM[2]; rl->delta = read->mismatches-MM[0]-MM[2]; 
		if (MM[1] > 0) fprintf(stdout,"newflag ");
		printf("newms:%d:%d:%d MM %d:%d added %d | NEWCIGAR_LR: ",MM[2],read->mismatches-MM[0],MM[1],MM[2],MM[0],added);
		for (i=0;i<rl->cigs;i++) fprintf(stdout,"%d%c ",rl->cigarlist[i]>>4,INT_CIGAROP[rl->cigarlist[i]&0xf]); //fprintf(stdout,"NEW "); 
		return 1;
	}
	return 0;
}

// what about 'D' events that are close to the variant we are going to evaluate, try with ignoring all 'deletions' separaterly..
// once indel event is added to new cigar, insertions are treated as 'M' and 'D' are ignored 
// indel realignment when indel is near beginning of read 
int realignread_RL(struct alignedread* read,REFLIST* reflist,VARIANT* varlist,int var,struct REALIGNMENT* rl)
{
	rl->cigs =0; //initilze
	int i=0,j=0,l1=0,l2=0,t=0,matches =0,carryover = 0,il=0,dl=0;
	int l11=read->l1-1,l21=read->l2-1; int added = -1; int MM[3] = { 0,0,0};
	int op=0,ol=0,p = 0; int indeladded =0; int b =0;  
	l1 = read->l1-1; l2 = read->l2-1; 
	int delta =0; if (varlist[var].type < 0) delta = -1*varlist[var].type; 
	if (varlist[var].simple == '0') delta = strlen(varlist[var].RA)-1;
	
	int firstposition = 0; op = read->cigarlist[0]&0xf;  ol = 0>>4;
	if ( (read->flag &16) ==16 && op == BAM_CSOFT_CLIP) firstposition = read->XC+1;  // if read->XC is non zero then we cannot extend the alignment into the clipped region (XC is # low quality bases clipped by BWA)

	//check from other direction right to left so we move over cigar string in that order
	for (i=read->cigs-1;i>=0;i--)
	{
	       	op = read->cigarlist[i]&0xf; ol = read->cigarlist[i]>>4; p = ol; 
		// only move into soft clip region if previous cigar is not 'M' extra condition added jan 2 2013 
		if (op == BAM_CSOFT_CLIP && i ==0 && rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0 && (rl->cigarlist[rl->cigs-1]>>4) >=10 && indeladded > 0)  
		{
			//fprintf(stdout,"ignoring soft clip region \t");
			rl->cigarlist[rl->cigs] = read->cigarlist[i]; rl->cigs +=1; l1 -= ol; l11 -= ol; continue;
		}
		else if (op == BAM_CMATCH || (op == BAM_CSOFT_CLIP && i ==0) || (op == BAM_CINS && indeladded > 0) || (op == BAM_CINS  && read->position + l2 == varlist[var].position + delta-1 && varlist[var].type > 0 && indeladded ==0) )
		{
			//if (op == BAM_CINS && indeladded > 0) MM[1]++;  
			for (t=0;t<ol;t++) // going from left to right in this case
			{
				if (read->position + l2 == varlist[var].position + delta-1 && indeladded ==0)  // indel insertion match in 'M' stretch
				{
					// changed 't < ol-1' -> 't < ol' 06/14/13
					if (varlist[var].type < 0 && varlist[var].simple == '1' && (rl->cigs > 0 || t < ol )) // do not add prior to start of 'M' 
					{
						fprintf(stdout,"cigs %d %d\n",rl->cigs,op);
						p = ol-t; carryover=0;
						if (t > 0 && rl->cigs ==0) 
						{ 
							rl->cigarlist[rl->cigs] = t<<4; rl->cigs +=1; 
						} 
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4;
							else { rl->cigarlist[rl->cigs] = t<<4; rl->cigs +=1; }
						}
						rl->cigarlist[rl->cigs] = ((-1*varlist[var].type)<<4) + 2; rl->cigs +=1; 
						printf("RL_deletion:%d:%d:%d %d:%d ",l1,t,ol,read->position,l2);
						l2 -= -1*varlist[var].type; added = 0; indeladded =1;
					}
					else if (varlist[var].type > 0 && varlist[var].simple == '1' &&  (rl->cigs > 0 || t < ol))
					{
						p = ol-t - varlist[var].type; carryover =0; if (p < 0) carryover = -1*p; 
						// carryover is used for storing # inserted bases over to next cigar
						if (t > 0 && rl->cigs == 0) rl->cigarlist[rl->cigs++] = t<<4; // 'M' to cigar 
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = (varlist[var].type<<4) + 1;
						for (j=varlist[var].type-1;j>=0 && l1-(varlist[var].type-j)+1 >= firstposition;j--) 
						{
							if (read->sequence[l1-(varlist[var].type-j)+1]  != varlist[var].AA[j+1]) MM[2]++;  
						}
						printf("ins<- %d carry %d ",l1,carryover);
						l1 -= varlist[var].type; t += varlist[var].type;
						added =0; indeladded = 1;
					}
					else if (varlist[var].simple == '0' && (rl->cigs > 0 || t < ol-1)) // complex indel
					{ 
						dl = strlen(varlist[var].RA)-1; il = strlen(varlist[var].AA)-1;
						p = ol-t; p -= il; carryover =0;if (p < 0) carryover = -1*p; 

						if (t > 0 && rl->cigs== 0) rl->cigarlist[rl->cigs++] = t<<4;
						else if (t > 0)
						{
							b = rl->cigarlist[rl->cigs-1]&0xf;
							if (b ==0) rl->cigarlist[rl->cigs-1] += t<<4; else rl->cigarlist[rl->cigs++] = t<<4;
						}
						rl->cigarlist[rl->cigs++] = (il<<4) + 1; rl->cigarlist[rl->cigs++] = (dl<<4) + 2;
						for (j=il-1;j>=0 && l1-(il-j)+1 >= firstposition;j--) 
						{
							if (read->sequence[l1-(il-j)+1]  != varlist[var].AA[j+1]) MM[2]++;  
						}
						printf("RLcomplexindel->%d p %d  %c carryover %d %d %d dl:%d il:%d ",l1,p,INT_CIGAROP[op],carryover,MM[2],t,dl,il);
						l1 -= il; t += il; l2 -= dl;
						added = 0; indeladded =1;
					}
				}
				if (t >=ol) break;
				if ( (op == BAM_CMATCH || op == BAM_CSOFT_CLIP || op == BAM_CINS) && carryover > 0) 
				{ 
					carryover--; p--; continue; 
				}
				if (l1 < firstposition) 
				{
					//fprintf(stdout,"| matches %d l1 %d first %d %d | ",matches,l1,firstposition,t);
					break;  // as soon as the read position l1 is smaller than 0/firstnonclipped base, break
				}
				if (added >= 0 && op == BAM_CSOFT_CLIP) added++;  // BUG FIXED 06/14/13 moved this after previous if condition

				// if t >= ol because of insertion -> these won't be counted jan 1 2013 
				//if (read->sequence[l1]  != reflist->sequences[read->tid][read->position+l2] && indeladded ==0) MM[1]++;
				if (read->sequence[l11-t]  != reflist->sequences[read->tid][read->position+l21-t] && indeladded ==0 && op != BAM_CINS)	MM[0]++;

				if (read->sequence[l1]  != reflist->sequences[read->tid][read->position+ l2]) 
				{
					if (op == BAM_CSOFT_CLIP && added >= 5 && carryover <=0 && varlist[var].type < 0) 
					{
						t = ol; added--; // subtract one since we don't count the last base
					}
					else if (indeladded > 0 || op != BAM_CMATCH) MM[2]++; 
				}
				else matches++; 
				//if (added >= 0) fprintf(stdout,"bbb added %d %c:%c",added,read->sequence[l1], reflist->sequences[read->tid][read->position+ l2]);
				if (t < ol) { l1--; l2--; } // BUG fixed, the start position of new cigar is affected ERR024183.19750915 
			}

			if (added < 0) 
			{ 
				if (rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0) rl->cigarlist[rl->cigs-1] += ol<<4;
				else { 	rl->cigarlist[rl->cigs] = read->cigarlist[i];  rl->cigs +=1; } 
			}
			else if (p > 0 && (op == BAM_CMATCH || op == BAM_CINS)) 
			{ 
				if (rl->cigs > 0 && (rl->cigarlist[rl->cigs-1]&0xf) ==0) rl->cigarlist[rl->cigs-1] += p<<4;
                                else { rl->cigarlist[rl->cigs] = p<<4; rl->cigs +=1; } 
			}
			else if (carryover ==0 && op == BAM_CSOFT_CLIP) 
			{ 
				printf("no carry S....%d %d %d ",p,carryover,added);
				b = rl->cigarlist[rl->cigs-1]&0xf;
                                if (added > 0 && b != 0)
                                {
                                        rl->cigarlist[rl->cigs] = added<<4; rl->cigs +=1;
                                }
                                else if (added > 0) rl->cigarlist[rl->cigs-1] += added<<4;

                                if (p-added > 0)
                                {
                                        rl->cigarlist[rl->cigs] = ((p-added)<<4) + 4; rl->cigs +=1; // remaining portion of soft clip
                                }
			}
			else if (carryover > 0 && op == BAM_CSOFT_CLIP && indeladded ==0) 
			{ 
				rl->cigarlist[rl->cigs] = (ol-carryover)<<4; rl->cigs +=1; carryover =0;
				// need to take care when carrayover is not from last iteration but current iteration...
			}
			if (op == BAM_CINS) l11 -= ol; else { l11 -= ol; l21 -= ol;} 
		}
		else if (op == BAM_CDEL) 
		{ 
			if (indeladded > 0) MM[1]++; // already added the indel event we are evaluating, so ignore extra deletion event
			else if (varlist[var].position < read->position + l2 - ol ) // added 07/17/13, deletion is ignored if the new variant being evaluated falls within the interval deleted on reference 
			{
				//fprintf(stdout,"del var %d %d \t",varlist[var].position,read->position+l2,ol);
				if (carryover > 0) {  rl->cigarlist[rl->cigs-1] -= carryover<<4; carryover = 0; } // probably never invoked
				rl->cigarlist[rl->cigs] = read->cigarlist[i]; rl->cigs +=1;
				l2 -= ol; 
			}
			l21 -= ol;
		}
		else if (op == BAM_CINS) // only invoked when insertion in old cigar is prior to indeladded event 
		{
			rl->cigarlist[rl->cigs] = read->cigarlist[i]; rl->cigs +=1; l1 -= ol; l11 -= ol;
		}
		else if (op == BAM_CSOFT_CLIP && i > 0)  // soft clip at end of read, copy as it is
		{ 
			rl->cigarlist[rl->cigs] = read->cigarlist[i]; rl->cigs +=1; l1 -= ol; l11 -= ol;
		}

	}
	if (carryover > 0)
	{
	       	op = rl->cigarlist[rl->cigs-1]&0xf; ol = rl->cigarlist[rl->cigs-1]>>4; 
 		if (op != BAM_CDEL && ol > carryover)  rl->cigarlist[rl->cigs-1] -= carryover<<4; 
		else return 0;
		//fprintf(stdout,"carry %d %d \t",carryover,rl->cigarlist[rl->cigs-1]);
	}
	
	// condition for passing read as better aligned as original read 
	int flag = 0;
	if (added ==0 && MM[2] < read->mismatches-MM[0]+MM[1]-1 && MM[2] <2) flag =1;
	if (added == 0 && MM[2] < read->mismatches-MM[0] && MM[2] < 3) flag = 2; 
	if ((added >= 2 && MM[2] < 1) || (added >= 5 && MM[2] <2)) flag = 3; 
	//fprintf(stdout,"flag %d MM %d:%d:%d %d %d\n",flag,MM[0],MM[1],MM[2],read->mismatches,added);
	if (flag > 0)
	{
		rl->newpos = read->position+l2+1; rl->added = added; rl->mismatches = MM[2]; rl->delta = read->mismatches-MM[0]-MM[2]; 
		if (MM[1] > 0) fprintf(stdout,"newflag ");
		printf("newms:%d:%d:%d MM %d:%d added %d newpos %d l2 %d NEWCIGAR_RL: ",MM[2],read->mismatches-MM[0],MM[1],MM[2],MM[0],added,read->position+l2+1,l2);
		for (i=0;i<rl->cigs/2;i++)  // reverse the order of cigar list
		{
			j = rl->cigarlist[i]; rl->cigarlist[i] = rl->cigarlist[rl->cigs-1-i]; rl->cigarlist[rl->cigs-1-i] = j; 
		}
		for (i=0;i<rl->cigs;i++) fprintf(stdout,"%d%c ",rl->cigarlist[i]>>4,INT_CIGAROP[rl->cigarlist[i]&0xf]); 
		//for (i=rl->cigs-1;i>=0;i--) fprintf(stdout,"%d%c ",rl->cigarlist[i]>>4,INT_CIGAROP[rl->cigarlist[i]&0xf]); 
		return 1;
	}
	return 0;
}

//fprintf(stdout,"l1 %d l2 %d %c %c \n",l1,read->position+l2,read->sequence[l1],reflist->sequences[read->tid][read->position+l2]);
