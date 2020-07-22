#!/usr/bin/awk -f
##########
#make sure that you have gawk installed so as to use asorti function

#  awk -f collisions_long0.6_v0.1.awk -v resfile="res.txt" c2bb6432-8b47-4279-a3db-8eaed1.sam > collisions_long_plus_end.txt
#version0.2: will generate merged.txt 
#version0.3: use and function for flags 256 and 272 - not primary alignment
#version0.4: remove maxpairdist flag, do not generated merged.txt
#version0.5: sort the alignments by CIGAR values
#version0.6: fix errors pointed out by Neva
#version0.6_v0.1: add the RefEnd position into the collisions.txt file




function get_refEnd(refStart, cigar) {
  refspan=0;
  where=match(cigar,/[0-9]+[M|D|N|=|X]/);
  while (where>0) {
    refspan += substr(cigar,where,RLENGTH-1);
    cigar=substr(cigar,where+RLENGTH);
    where=match(cigar,/[0-9]+[M|D|N|=|X]/);
  }
  refEnd=refStart+refspan-1;
  return refEnd;
}



function get_cigStart(strand,cigar)
{ 
  mystr = and(strand,16)
  if (mystr==0) {
    if (cigar ~/^[0-9]+S/) {
      split(cigar,CIGAR,"S");
      cigstart=CIGAR[1];
      
    }
    else if (cigar ~/^[0-9]+H/) {
      split(cigar,CIGAR,"H");
      cigstart=CIGAR[1];

    }
    else if (cigar ~/^[0-9]+M/) {
      cigstart=0
    }   
  }
  else if (mystr==16) {
    if (cigar ~/[0-9]+S$/) {
      where = match(cigar,/[0-9]+S$/);
      cigstart=substr(cigar,where,RLENGTH-1);
      
    }
    else if (cigar ~/[0-9]+H$/) {
      where = match(cigar,/[0-9]+H$/);
      cigstart=substr(cigar,where,RLENGTH-1);
    }
    else if (cigar ~/[0-9]+M$/) {
      cigstart=0
    }   
  }
  return cigstart
}


BEGIN{
  OFS="\t";
  tottot = 0; # will count first non-group
  prev=""
  count=0;
  count_unmapped=0;
  count_single=0;
  count_2frags=0;
  count_gt2frags=0;
  PROCINFO["sorted_in"]="@val_num_asc";
}


#if count==0, it's the 1st line
{
  #make sure SAM file is sorted by read name

  if (/^@/) {next;}
  else {
    if ($1==prev) {
      #not new record
      flag=and($2,256);
      if (flag==256) {next;}
      else {
        count++;
        c[count] = $0;
      }
      
    }
    else {  
      #dealing with a new record, need to check if count==0
      if (count==0) { #1st line of file
        #initialize variables
        prev=$1;
        count = 1;
        c[count] = $0;
      }
      else {
        tottot++;
        #need to process prev SAM info
        for (i=1; i<=count; i++) {
          split(c[i],tmp,"\t");
          name[i] = tmp[1];
          str[i] = and(tmp[2],16);
          chr[i] = tmp[3];
          pos[i] = tmp[4];
          m[i] = tmp[5];
          cig[i] = tmp[6];
        }
        if (count==1 && tmp[2]==4) {count_unmapped++;}
        else {
          if (count==1) {
            count_single++;
            printf("%s %d %s %d %d %s", name[1], str[1], chr[1], pos[1], m[1], cig[1]);
            printf("\n")
          }
          else if (count==2) {count_2frags++;}
          else if (count>2) {count_gt2frags++;}

          if (count>=2) {
            for (i=1; i<=count; i++) {
              #sort SAM info based on CIGAR values
              start=get_cigStart(str[i],cig[i]);
              #CigStart is a dictionary, 9 zero padding of the indexes
              CigStart[sprintf("%0.9d", start)] = c[i];    
            }

            asorti(CigStart,sorted_cig);
            n = length(sorted_cig);
            for (i=1; i<=n; i++) {
              split(CigStart[sorted_cig[i]],stmp);
              name[i] = stmp[1];
              str[i] = and(stmp[2],16);
              chr[i] = stmp[3];
              pos[i] = stmp[4];
              m[i]   = stmp[5];
              cig[i] = stmp[6];
              refend[i] = get_refEnd(pos[i], cig[i])
              printf("%s %d %s %d %d %s ", name[i], str[i], chr[i], pos[i], m[i], cig[i]);
              printf("%s %d %s %d %d %s ", name[i], str[i], chr[i], refend[i], m[i], cig[i]);
            }   
            printf("\n");
          }
        }
        #reset variables
        delete c;
        delete CigStart;
        delete sorted_cig;
        delete refend;

        prev=$1;
        count = 1;
        c[count]=$0;
      }

    }

  }
  
}

END {
  tottot++;
  #need to process prev SAM info
  for (i=1; i<=count; i++) {
    split(c[i],tmp,"\t");
    name[i] = tmp[1];
    str[i] = and(tmp[2],16);
    chr[i] = tmp[3];
    pos[i] = tmp[4];
    m[i] = tmp[5];
    cig[i] = tmp[6];
  }
  if (count==1 && tmp[2]==4) {count_unmapped++;}
  else {
    if (count==1) {
      count_single++;
      printf("%s %d %s %d %d %s", name[1], str[1], chr[1], pos[1], m[1], cig[1]);
      printf("\n")
    }
    else if (count==2) {count_2frags++;}
    else if (count>2) {count_gt2frags++;}

    if (count>=2) {
      for (i=1; i<=count; i++) {
        #sort SAM info based on CIGAR values
        start=get_cigStart(str[i],cig[i]);
        #CigStart is a dictionary, 9 zero padding of the indexes
        CigStart[sprintf("%0.9d", start)] = c[i];
      }

      asorti(CigStart,sorted_cig);
      n = length(sorted_cig);
      for (i=1; i<=n; i++) {
        split(CigStart[sorted_cig[i]],stmp);
        name[i] = stmp[1];
        str[i] = and(stmp[2],16);
        chr[i] = stmp[3];
        pos[i] = stmp[4];
        m[i]   = stmp[5];
        cig[i] = stmp[6];
        refend[i] = get_refEnd(pos[i], cig[i]);
        printf("%s %d %s %d %d %s ", name[i], str[i], chr[i], pos[i], m[i], cig[i]);
        printf("%s %d %s %d %d %s ", name[i], str[i], chr[i], refend[i], m[i], cig[i]);
      }   
      printf("\n");
    }
  }

  #export the counts to res.txt file
  printf("%s %s %s %s %s\n", "total_reads","unmapped","singleton","two_fragments","gt2fragments") >> resfile;
  printf("%d %d %d %d %d\n", tottot, count_unmapped, count_single, count_2frags, count_gt2frags) >> resfile;

}