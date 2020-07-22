#!/usr/bin/awk -f
##########

#2_collision2pairs_v0.3.awk

#awk -f collision2pairs.awk -v pair_radius=infinity  collisions.txt > merged.txt  #return all  pairs
#awk -f collision2pairs.awk -v pair_radius=1  collisions.txt > merged.txt
#version0.1: input collisions.txt include ref end position as well, will integrate parameter of pair distance and pair radius
#pair radius=3 will count 1<=pairDistance<=3
#pair radius=infinity will count all pairs

#version0.2: will only consider the pairs derived from ligation junction, thus pair_radius=1 and pair_radius=infinity should be discounted
#set pair_radius=1 as default
#for instance, if a read has 3 chimeric fragments, A, B, and C, which maps to the reference start and end, the corresponding reference position would be
#A1,A2; B1,B2; C1,C2; we will only consider A2B1, B2C1

#version0.3: include a 'junctionOnly argument' so that if 'junctionOnly=TRUE' will only consider ligation junction fragments
#otherwise if false,  can consider pair_radius=infinity or pair_radius=1
#awk -f collision2pairs_v0.3.awk -v junctionOnly=1 -v pair_radius=1  collisions.txt  | sed $'s/ /\t/g' > merged.txt
#awk -f collision2pairs_v0.3.awk -v junctionOnly=0 -v pair_radius=1  collisions.txt  | sed $'s/ /\t/g' > merged.txt

function less_than(s1,c1,p1,s2,c2,p2)
{

  if (c1 < c2) return 1;
  if (c1 > c2) return 0;
  # c1 == c2
  if (s1 < s2) return 1;
  if (s1 > s2) return 0;
  # s1 == s2 && c1 == c2
  if (p1 < p2) return 1;
  if (p1 > p2) return 0;
  # all are equal, doesn't matter
  return 1;
}



{
	if (NF<7) {next;}
	#only look at lines with chimeric reads
	else {
		count=NF/12;
		split($0,tmp);
		for (i=0;i<count;i++) {
			name[i]= tmp[i*12+1];
			str[i]= and(tmp[i*12+2],16);
	        chr[i]= tmp[i*12+3];
	        pos[i]= tmp[i*12+4];
	        pos_end[i]= tmp[i*12+10];
	        m[i]= tmp[i*12+5];
	        cig[i]= tmp[i*12+6];	        
		}

		for (x=0;x<count-1;x++) {
			for (y=0;y<count;y++) {
				if (y-x<=0) {pass;}
				else {
					read1=x;
					read2=y;
					if (junctionOnly) {
						#junctionOnly=1 (true), eg: #eg: ABC with A1,A2; B1,B2; C1,C2 -- if pair_radius==1, will consider A2B1,B2C1; if pair_radius=="infinity", will consider will consider A2B1,B2C1, and A2C1 
						if (pair_radius=="infinity") {
							#use short format
							if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
			                	print name[read1],str[read1],chr[read1],pos_end[read1],str[read2],chr[read2],pos[read2], m[read1],m[read2]; #A2B1
			                }
			                else {
			                	print name[read2],str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos_end[read1], m[read1],m[read1]; #B1A2
			                }

						}
						else  {
							if (read2-read1<=pair_radius) {
								#use short format
								if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
				                	print name[read1],str[read1],chr[read1],pos_end[read1],str[read2],chr[read2],pos[read2], m[read1],m[read2]; #A2B1
				                }
				                else {
				                	print name[read2],str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos_end[read1], m[read1],m[read1]; #B1A2
				                }

							}
							else {pass;}

						}
						

					}
					else {
						#junctionOnly=0 (false), eg: #eg: ABC with A1,A2; B1,B2; C1,C2 -- if pair_radius==1, will consider A1B1,A1B2,A2B1,A2B2 PLUS B1C1,B1C2,B2C1,B2C2; if pair_radius=="infinity", will consider will consider A1B1,A1B2,A2B1,A2B2 PLUS B1C1,B1C2,B2C1,B2C2, and A1C1,A1C2,A2C1,A2C2. 

						if (pair_radius=="infinity") {
							#use short format
							if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
			                	print name[read1],str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2], m[read1],m[read2];#A1B1
			                	print name[read1],str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos_end[read2], m[read1],m[read2]; #A1B2
			                	print name[read1],str[read1],chr[read1],pos_end[read1],str[read2],chr[read2],pos[read2], m[read1],m[read2]; #A2B1
			                	print name[read1],str[read1],chr[read1],pos_end[read1],str[read2],chr[read2],pos_end[read2], m[read1],m[read2]; #A2B2
			                }
			                else {
			                	print name[read2],str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos[read1],m[read2],m[read1]; #B1A1
			                	print name[read2],str[read2],chr[read2],pos_end[read2],str[read1],chr[read1],pos[read1], m[read2],m[read1]; #B2A1
			                	print name[read2],str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos_end[read1], m[read1],m[read1]; #B1A2
			                	print name[read2],str[read2],chr[read2],pos_end[read2],str[read1],chr[read1],pos_end[read1], m[read1],m[read1]; #B2A2
			                }

						}

						else {
							if (y-x<=pair_radius) {
								if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
			                		print name[read1],str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2], m[read1],m[read2]; #A1B1
			                		print name[read1],str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos_end[read2], m[read1],m[read2]; #A1B2
			                		print name[read1],str[read1],chr[read1],pos_end[read1],str[read2],chr[read2],pos[read2], m[read1],m[read2]; #A2B1
			                		print name[read1],str[read1],chr[read1],pos_end[read1],str[read2],chr[read2],pos_end[read2], m[read1],m[read2]; #A2B2
			                	}
			                	else {
			                		print name[read2],str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos[read1],m[read2],m[read1]; #B1A1
			                		print name[read2],str[read2],chr[read2],pos_end[read2],str[read1],chr[read1],pos[read1], m[read2],m[read1]; #B2A1
			                		print name[read2],str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos_end[read1], m[read1],m[read1]; #B1A2
			                		print name[read2],str[read2],chr[read2],pos_end[read2],str[read1],chr[read1],pos_end[read1], m[read1],m[read1]; #B2A2
			              		}
							}
							else {pass;}

						}
					}


					

				}



			}
		}
	}
}