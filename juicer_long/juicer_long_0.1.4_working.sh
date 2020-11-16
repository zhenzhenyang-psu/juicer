#!/bin/bash
##########

export LC_ALL=C

topDir=$(pwd)
# genome ID, default to human, can also be set in options
genomeID="hg19"

# do not include fragment delimited maps by default
nofrag=1
#set default minimap2 command
minimap2Path="minimap2"

#set junctionOnly=1 as the default: only look at contacts from ligation junction, eg: read with ABC fragment with reference start and end as A1,A2; B1,B2; C1,C2; will consider A2B1, B2C1 if pair_radius=1

fragDist=100 #use fragDistance of 100bp as default
usegzip=1

#0.1.4 version tries to split alignment and SAM processing to generate pairs into two stages - \
#this is because in many cases, we want to process pairs in different ways and we want to start from alignment (SAM) instead of from scratch

#add two separate stage: early_align and early_fromSam; keep the original early stage (all the way from fastq to merged.txt)
#early_align will do just the alignment process of the early stage - i.e. stop once SAM is produced
#early_fromSam will start form SAM and process SAM to pairwise contacts


## Read arguments 
usageHelp="Usage: ${0##*/} [-d topDir] [-g genomeID] [-m minimap2Path] [-s site]  [-y restriction site file]\n[-z reference genome file] [-D Juicer scripts directory]\n[-i fragDist] [-r pair_radius] [-S stage] [-t threads] [-b fragment_distance] [-h] [-f]"
genomeHelp="* [genomeID] must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \n  \"$genomeID\"); alternatively, it can be defined using the -z command"
excludeHelp="* -f: include fragment-delimited maps in hic file creation"
minimap2Help="* -m: full path of minimap2, default is minimap2" 
refSeqHelp="* -z [reference genome file]: enter path for reference sequence file"
siteHelp="* [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" or \"NlaIII\" \n  (default \"$site\")"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
stageHelp="* [stage]: must be one of \"merge\", \"final\", \"postproc\", or \"early\".\n    -Use \"earlyAln\" for an early exit right after alignment is done \n    -Use \"earlyFromSam\" to pick up from Sam, and process SAM to generate pairs\n    -Use \"early\" for an early exit, including map fastq reads to generate SAM\n     and process into pairs based on junction or distance\n    -Use \"merge\" when each single fastq has generated pair information,\n     but results from these multiple fastq files were not combined.\n    -Use \”final\” to generate .hic files and check .res.txt\n    --Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed."
pairRadiusHelp="-r [pair_radius]: use \"-r 1\" if considering pairDistance of 1 only; use \"-r infinity\" if considering all interactions; \"-r 3\" will consider consider loci with pairDist<=3"
fragDistHelp="-b [fragment_distance]: use \"-b 5000\" if considering fragment distance of 5000bp as the cutoff to count as a valid fragment; \"-b 100\" only considers distance btw mapped positions >100bp as a contact "
topDirHelp="* -d [topDir] is the top level directory (default\n  \"$topDir\")\n  [topDir]/fastq must contain the fastq files\n  make sure that input read has to be in .fastq or .fastq.gz format, .fq format woun't work\n  [topDir]/aligned will be created for the final alignment"
juicerDirHelp="* -D [juicerDir] is the juicer home directory \n  [juicerDir]/references must contain the refseq fasta\n  [juicerDir]/scripts must exist"
threadsHelp="* -t threads "
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
	echo -e "$usageHelp"
	echo -e "$genomeHelp"
	echo -e "$refSeqHelp"
	echo -e "$siteHelp"
	echo -e "$siteFileHelp"
    echo -e "$stageHelp"
    echo -e "$pairRadiusHelp"
    echo -e "$fragDistHelp"
    echo -e "$minimap2Help"
	echo -e "$topDirHelp"
	echo -e "$juicerDirHelp"
	echo -e "$threadsHelp"
	echo "$excludeHelp"
	echo "$helpHelp"
	exit "$1"
}

while getopts "d:g:s:y:b:r:m:z:m:S:D:t:f:h" opt; do
	case $opt in
	d) topDir=$OPTARG ;;
	g) genomeID=$OPTARG ;;
	s) site=$OPTARG ;;
	y) site_file=$OPTARG ;;
	r) pairRadius=$OPTARG ;;
	b) fragDist=$OPTARG ;;
	m) minimap2Path=$OPTARG;;
	z) refSeq=$OPTARG ;;
	S) stage=$OPTARG ;;
	D) juiceDir=$OPTARG ;;
	t) threads=$OPTARG ;;
	f) nofrag=0 ;; #use fragment maps, default nofrag=1 (when not set)
	h) printHelpAndExit 0;;
	[?]) printHelpAndExit 1;;
	esac
done

echo "pair radius is ${pairRadius}"




if [ ! -z "$stage" ]
then
    case $stage in
    	early) earlyexit=1 ;;
		earlyAln) earlyAln=1;;
		earlyFromSam) earlyFromSam=1;;
        merge) merge=1 ;;        
        final) final=1 ;;
	postproc) postproc=1 ;; 
        *)  echo "$usageHelp"
	    echo "$stageHelp"
	    exit 1
    esac
fi

#according to current order of the stage, if ! -z $early (if early==1), earlyAln=0, earlyFromSam=0
if [ ! -z $earlyexit ] 
then
	earlyAln=0
	earlyFromSam=0
fi

#if earlyAln==1, earlyFromSam=0 && earlyexit=0
if [ ! -z $earlyAln ] 
then
	earlyFromSam=0
	earlyexit=0
fi

#if earlyFromSam==1, earlyexit=0 && earlyAln=0
if [ ! -z "$earlyFromSam" ] 
then
	earlyAln=0
	earlyexit=0
fi

## Set reference sequence based on genome ID
if [ -z "$refSeq" ] # if $refseq is empty: if -z not set,  you have to set -g to be hg19, or hg38, or mm9, or mm10
then 
	case $genomeID in
	mm9) refSeq="${juiceDir}/references/Mus_musculus_assembly9_norandom.fasta";;
	mm10) refSeq="${juiceDir}/references/Mus_musculus_assembly10.fasta";;
	hg38) refSeq="${juiceDir}/references/Homo_sapiens_assembly38.fasta";;
	hg19) refSeq="${juiceDir}/references/Homo_sapiens_assembly19.fasta";;
	
	*)  echo "$usageHelp"
			echo "$genomeHelp"
			exit 1
	esac
fi

## Set threads for sending appropriate parameters to cluster and string for BWA call
if [ ! -z "$threads" ]
then
	threadstring="-t $threads"
else
	threads="$(getconf _NPROCESSORS_ONLN)"
	threadstring="-t $threads"
fi


# ## commented this section because long read will always have ligation junctions except non-chimeric alignment, which will be handled by looking at number of fragments for each read
# if [ -z "$ligation" ]
# then
#     case $site in
# 	HindIII) ligation="AAGCTAGCTT";;
# 	DpnII) ligation="GATCGATC";;
# 	MboI) ligation="GATCGATC";;
# 	NcoI) ligation="CCATGCATGG";;
# 	none) ligation="XXXX";;
# 	*)  ligation="XXXX"
# 	    echo "$site not listed as recognized enzyme. Using $site_file as site file"
# 	    echo "Ligation junction is undefined"
#     esac
# fi


## If DNAse-type or MNase-type experiment, no fragment maps
if [ "$site" == "none" ]
then
    nofrag=1;
fi

if [ -z "$site_file" ]
then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [ ! -e "$site_file" ] && [ "$site" != "none" ]
then
    echo "***! $site_file does not exist. It must be created before running this script."
    exit 1
fi


## Directories to be created and regex strings for listing files
fastqdir=${topDir}"/fastq/*.fastq*"
splitdir=${topDir}"/splits"

outputdir=${topDir}"/aligned"
tmpdir=${topDir}"/HIC_tmp"
## Check that fastq directory exists
if [ ! -d "$topDir/fastq" ]; then
	echo "Directory \"$topDir/fastq\" does not exist."
	echo "Create \"$topDir/fastq\" and put fastq files to be aligned there"
	echo "Type \"juicer_long.sh -h\" for help"
	exit 1
else 
	echo "(-: Looking for fastq files...fastq files exist"
fi

## Create output directory, only if not in merge, final, or postproc stages
if [[ -d "$outputdir" && -z "$final" && -z "$merge" && -z "$postproc" ]] 
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    echo "***! Type \"juicer.sh -h \" for help"
    exit 1	
else
    if [[ -z "$final" && -z "$merge" && -z "$postproc" ]]; then
        mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 1; } 
    fi
fi		


## Create split directory
if [ -d "$splitdir" ]; then
    splitdirexists=1
else
    mkdir "$splitdir" || { echo "***! Unable to create ${splitdir}, check permissions." ; exit 1; }
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ] && [ -z "$final" ] && [ -z "$postproc" ]; then
	mkdir "$tmpdir"
	chmod 777 "$tmpdir"
fi




## Arguments have been checked and directories created. Now begins
## the real work of the pipeline
testname=$(ls -l ${fastqdir} | awk 'NR==1{print $9}')
if [ "${testname: -3}" == ".gz" ]
then
	readlong=${splitdir}"/*.fastq.gz"
	usegzip=1
	echo "$readlong is in fastq.gz format, usegzip=1"
else
	readlong=${splitdir}"/*.fastq"
	usegzip=0
	echo "$readlong is in fastq format, usegzip=0"
fi

#if doing postproc but without having output dir and outputdir/header file, exit and report a warning
if [ ! -z $postproc ] && [ ! -d ${outputdir} ]; then
	echo "aligned directory doesn't exist, do merge or final before proceeding with postproc"
	exit 1
else
	if [ ! -z $postproc ] && [ ! -s $outputdir/merged_sort.txt ]; then
		echo "aligned/merged_sort.txt does not exist or is empty, please check merge before postprocessing "
		exit 1
	elif [ ! -z $postproc ] && [ ! -s $outputdir/inter.hic ]; then
		echo "aligned/inter.hic does not exist, please run final to generate .hic file before postprocessing "
		exit 1
	fi
fi

#header records the current time of the start run
headfile=${outputdir}/header
date > $headfile

# Get version numbers of all software
${minimap2Path} 2>&1 | awk '$1=="Version:"{printf(" minimap2 %s; ", $2)}' >> $headfile 
#echo -ne: remove newline character at the end 
echo -ne "$threads threads; " >> $headfile
java -version 2>&1 | awk 'NR==1{printf("%s; ", $0);}' >> $headfile 
## Juicer directory, contains scripts/, references/, and restriction_sites/
${juiceDir}/scripts/juicer_tools -V 2>&1 | awk '$1=="Juicer" && $2=="Tools"{printf("%s; ", $0);}' >> $headfile
echo "$0 $@" >> $headfile
#echo "$0 $@" >> $headfile  ask Neva: what does $@ mean?



## if in earlyexit stage, will align reads and variables for other stages will be all 0
if [ -z $merge ] && [ -z $final ] && [ -z $postproc ]
then
	#if in earlyAln or early stage--earlyFromSam will =0, do following minimap2 alignment
	if [ -z $earlyFromSam ]
	then

		if [ ! $splitdirexists ] #if split directory not exists
		then
			echo -e "(-: Aligning files matching $fastqdir\n to genome $genomeID"
			echo "(-: Created $splitdir and $outputdir."
			ln -s $fastqdir $splitdir/.
		else
			echo -e "--- Using already created files in $splitdir\n"
		fi

		for i in ${readlong}
		do
			#i=/public/home/yangzhzh/projects/0_aiden_lab/6_brian/2_run_juicerLong/COLA12236/splits/COLA12236_ont.fastq
			fqName=$(basename $i) #fqName: single.fastq
			name=${i} #/public/home/yangzhzh/projects/0_aiden_lab/6_brian/2_run_juicerLong/COLA12236/splits/COLA12236_ont.fastq
			jname=$(basename $i) #basename: single.fastq
			#Align readlong
			#fqPrefix=${fqName%.fastq*}
			#ext: refers to the extension "fastq" or "gz"
			ext=${fqName##*.} 
			
			echo -e "ext (extension) is: $ext; basename is: $jname; site_file is: $site_file" >> $headfile
			if [ $ext == "gz" ]
			then
				usegzip=1
			fi

			if [ $ext == "fastq" ]
			then
				usegzip=0
			fi
			#source ${juiceDir}/scripts/common/countligations_long.sh
			echo "readname is $fqName"
			samName=$i".sam"
			if [ -s ${samName} ]; then
				echo "sam alignment file already generated, move onto the next step of creating collision and pairs."
			else
				echo "minimap2 path is ${minimap2Path}"
				echo "Running command ${minimap2Path} -ax map-ont -L $threadstring $refSeq $i > $samName"
				echo "Running command ${minimap2Path} -ax map-ont -L $threadstring $refSeq $i > $samName" >> $headfile 
				${minimap2Path} -ax map-ont -L $threadstring $refSeq $i > $samName
			fi
		done
	
	else 
		#if in earlyexit or earlyFromSam stage -- earlyAln will =0, process SAM alignment to pairs
		if [ -z $earlyAln ]
		then

			if [ -s $splitdir/$jname".collisions.txt" ] && [ -s $splitdir/$jname".res.txt"  ] ; then
				echo "collilsions and res.txt for $jname, move onto the next step of resolving pairs."
			else
				awk -f ${juiceDir}/juicer_long/1_collisions_long0.6_v0.2.awk \
					-v resfile=$i".res.txt" $samName > $i".collisions.txt" 
				echo "The command to generate collisions.txt from SAM file is:" >> $headfile
				echo "awk -f ${juiceDir}/juicer_long/1_collisions_long0.6_v0.2.awk -v resfile=$i".res.txt" $samName > $i".collisions.txt" " >> $headfile
			fi
			
			echo "pair radius is ${pairRadius}, fragment distance value is ${fragDist}" 
			echo "pair radius is ${pairRadius}, fragment distance value is ${fragDist}" >> $headfile
			echo "The command to generate merged pairs from collision file is :" >> $headfile
			echo "awk -f ${juiceDir}/juicer_long/2_collision2pairs_v0.4.awk -v pair_radius=${pairRadius}  -v fragmentDistance=${fragDist} $i".collisions.txt"  >  $i".merged.txt" " >> $headfile
			awk -f ${juiceDir}/juicer_long/2_collision2pairs_v0.4.awk \
			-v pair_radius=${pairRadius} -v fragmentDistance=${fragDist} $i".collisions.txt"  >  $i".pairs.txt" 

			if [ $? -ne 0 ]; then
				echo "***! Failure during chimera handling of $i"
				exit 1
			else
				echo "(-: succesfully converted SAM to pairs"
			fi

			if [ "$site" == "none" ]; then
				awk '{printf("%s %d %s %d %d %d %s %d %d %d %d\n",$1, $2, $3, $4, 0, $5, $6, $7, "1", $8, $9)}' $i".pairs.txt" > $i".pair.frag.txt"
			fi

			if [ $? -ne 0 ]; then
				echo "***! Failure during fragment assignment of ${i}.pairs.txt "
				exit 1
			else
				echo "(-: succesfully assigned fragment"
			fi
		fi
	
	fi
fi

#if in the merge stage, final and postprocess will be 0
if [ -z $earlyexit ] || [ -z $earlyAln ] || [ -z $earlyFromSam ]
then # if in merge stage, final and postproc will be 0
	if [ -z $final ] &&  [ -z $postproc ]
	then
		if [ -s "${outputdir}/merged_sort.txt" ]; then
			echo "merged_sort.txt file already created, move onto the next step of creating stats."
		else
			# sort by chromosome,  strand, and position   if no restriction site;
			#otherwise, sort by chromosome, frag, strand, and position                 
			cat ${splitdir}/*.pair.frag.txt > ${outputdir}/merged_frag.txt
	    	sort -T $tmpdir -k3,3d -k7,7d -k2,2n -k6,6n -k4,4n -k8,8n  $splitdir/*.pair.frag.txt  > $outputdir/merged_sort.txt
	    	echo "sort -T $tmpdir -k3,3d -k7,7d -k2,2n -k6,6n -k4,4n -k8,8n  $splitdir/*.pair.frag.txt  > $outputdir/merged_sort.txt" >> $headfile
		    if [ $? -ne 0 ]; then
				echo "***! Some problems occurred somewhere in creating sorted align files."
				exit 1
		    else
		        echo "(-: Finished sorting all sorted files into a single merge."
		        rm -rf ${tmpdir}
		    fi
		fi
	fi
fi

#STATISTICS and create .hic

#if in early stage, we don't do the following; otherwise earlyexit=0, will do .hic file creation and postprocessing 
if [ -z "$earlyexit" ] || [ -z "$earlyAln" ] || [ -z "$earlyFromSam" ] #if in earlyexit stage, $earlyexit=1, this will be false, condition isn't meet and exit
then 
	#if in other stages, like merge, final, postproc=0, check statistics and create .hic file 
	if [ -z $postproc ]  #not early exit, not postprocess, but in merge (final)
	then 
		export _JAVA_OPTIONS=-Xmx16384m
	    export LC_ALL=en_US.UTF-8 
	    head -n2 $headfile | tail -n1| awk '{printf"%-1000s\n", $0}' > $outputdir/inter.txt;
	    cat $splitdir/*.res.txt | awk -f ${juiceDir}/juicer_long/stats_sub_long_v2.awk >> $outputdir/inter.txt
	    cp $outputdir/inter.txt $outputdir/inter_30.txt
		${juiceDir}/scripts/common/juicer_tools pre -s $outputdir/inter.txt -q 1 $outputdir/merged_sort.txt $outputdir/inter.hic ${genomeID}
		${juiceDir}/scripts/common/juicer_tools pre -s $outputdir/inter_30.txt  -q 30 $outputdir/merged_sort.txt $outputdir/inter_30.hic ${genomeID}
	fi
	#postprocessing will be done in all stages, except in early exit stage
	${juiceDir}/scripts/common/juicer_postprocessing.sh -j ${juiceDir}/scripts/common/juicer_tools -i ${outputdir}/inter_30.hic  -g ${genomeID}
fi



#CHECK THAT PIPELINE WAS SUCCESSFUL
export early=$earlyexit
export splitdir=$splitdir
source ${juiceDir}/juicer_long/check_long_v3.sh


