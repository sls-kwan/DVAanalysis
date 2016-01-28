#!/bin/bash

string1="_1.fastq"
string2="_2.fastq"
string3=".fastq"
uscore="_"
countfile=".count"
starOut="_Aligned.sortedByCoord.out.bam"
currDir=`pwd`/
outpdf="output.pdf"
outtxt="output.txt"
while IFS='' read -r line || [[ -n "$line" ]]; do
	if [ ! -f $line$string3 -o ! -f $line$string1 ]
	   then
	      echo "Files not download"
	fi
	if [ -f $line$string2 ]
	   then
	   echo "Doing paired read alignment"
	   /tmp/notbackedup/software/STAR-STAR_2.4.0g1/bin/Linux_x86_64/STAR  --runMode alignReads   --readFilesIn $line$string1 $line$string2  --genomeDir ~/mm10ReferenceDir --outReadsUnmapped fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattributes All --runThreadN 8 --outFileNamePrefix $line$uscore  --chimSegmentMin 15 --chimJunctionOverhangMin 15
	else
	    /tmp/notbackedup/software/STAR-STAR_2.4.0g1/bin/Linux_x86_64/STAR  --runMode alignReads   --readFilesIn $line$string3 --genomeDir ~/mm10ReferenceDir --outReadsUnmapped fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattributes All --runThreadN 8 --outFileNamePrefix $line$uscore  --chimSegmentMin 15 --chimJunctionOverhangMin 15
	fi
	htseq-count -q -f bam $line$starOut ~/Documents/Mus_musculus.GRCm38.83.gtf > $line$countfile
	#java -jar /tmp/notbackedup/software/picard-tools-1.96/CollectRnaSeqMetrics.jar I=$line$starOut REFERENCE_SEQUENCE="/mnt/smb/ngs/References/mm10/mm10_ERCC.fa" STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND CHART_OUTPUT="$currDir$line$uscore$outpdf" METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT="$currDir$line$uscore$outtxt" STOP_AFTER=100000000 REF_FLAT="/mnt/smb/ngs/References/mm10/mm10_refseq.refFlat"

done < "$1"
