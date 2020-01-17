# THROUGHOUT THIS SCRIPT change all paths to your input and reference database files
# Current script works on both microbe and ixodes mothur v. 1.43.0
# to run in interactive mode in unix/linux in a shell: mothur; windows double click the mothur.exe file
# batch mode: mothur mothur_processing.sh
# in case make.file() or make.contigs() dosen't work run bash: gunzip *.gz  
# make tab deliminated 16S.files (https://www.mothur.org/wiki/File_File)

##############
#OUTPUT FILES
##############
# If completed successfully, the following files will be generated
# 16S.final.fasta :: Chimera filtered representative sequences
# 16S.final.count_table :: Frequency OTU table for each sample
# 16S.final.taxonomy :: Representative OTU taxonomy assignments
# raw.count :: Number of raw reads per sample
# merged.count :: Number of reads that were successfully merged (into contigs)
# filter.count :: Number of merged reads that pass quality filtering thresholds

################################
#MERGE, QUALITY FILTER RAW DATA
################################

#merge, quality filter, and classify your reads
#set log file name for current session
set.logfile(name=test.log, append=T) 
#set raw file path, this is where your files will be written to
system(ls *fastq.gz | parallel 'gzip -d {}') # version 1.43.0 cannot handle gzipped files, unzip before continuing
make.file(inputdir=/home/lymelab/Desktop/mann/test, prefix=16S) # CHANGE ME to your raw fastq path
#get number of raw reads in each sample
system(ls *R1* | sed 's/_.*//' > sample.ids)
system(ls *R1* | while read line; do cat $line | grep "^@" | wc -l; done > raw.count)
system(paste raw.count sample.ids > temp)
system(mv temp raw.count)
system(rm sample.ids)
#merge paired end reads
make.contigs(file=16S.files, processors=12, oligos=primers.txt) # merge reads and trim off primers
system(ls *fastq | parallel 'gzip {}') # rezip files
#count reads that were sucessfully merged
system(awk '{print $2}' 16S.contigs.groups | sort | uniq -c | sed 's/^[ \t]*//' > merged.count)
#summarize quality of sequences in merged reads
summary.seqs(fasta=16S.trim.contigs.fasta) 
#filter out any reads after quality filtering that are not at least 275 bp long
screen.seqs(fasta=16S.trim.contigs.fasta, group=16S.contigs.groups, summary=16S.trim.contigs.summary, maxambig=0, maxlength=275) 
summary.seqs(fasta=16S.trim.contigs.good.fasta)
#number of reads post filter
system(awk '{print $2}' 16S.contigs.good.groups | sort | uniq -c | sed 's/^[ \t]*//' > filter.count)
#dereplicate sequences so they are easier to classify
unique.seqs(fasta=16S.trim.contigs.good.fasta)
count.seqs(name=16S.trim.contigs.good.names, group=16S.contigs.good.groups)
summary.seqs(count=16S.trim.contigs.good.count_table)

#######################################
#OPTIONAL: GENERATE REFERENCE DATABASE 
#######################################
#CHANGE paths below to where your reference database is
#generate reference database for the V4 region, only needs to be run once
summary.seqs(fasta=/home/lymelab/Desktop/referenceDB/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_full_align.fasta)
#extract V4 region, allow for 3 mismatches in either primer
#this first round will drop many sequences so use the summary to get an estimate of the actual base position of the V4 region and use this in the next round
pcr.seqs(fasta=/home/lymelab/Desktop/referenceDB/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_full_align.fasta, oligos=primers.txt, processors=12, pdiffs=3, rdiffs=3)
summary.seqs(fasta=/home/lymelab/Desktop/referenceDB/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_full_align.pcr.fasta)
#now get actual position, use position info from above command output (ezbiocloud_full_align.pcr.summary)
pcr.seqs(fasta=/home/lymelab/Desktop/referenceDB/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_full_align.fasta, start=11969, end=14973, keepdots=F, processors=12)
system(mv /home/lymelab/Desktop/referenceDB/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_full_align.pcr.fasta /home/lymelab/Desktop/referenceDB/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_full_align.v4.fasta)
summary.seqs(fasta=/home/lymelab/Desktop/referenceDB/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_full_align.v4.fasta)

#########################################
#SAMPLE CLUSTERING, TAXONOMIC ASSIGNMENT
#########################################

#align reads to reference database, clean up those that do not align
align.seqs(fasta=16S.trim.contigs.good.unique.fasta, reference=/home/lymelab/Desktop/referenceDB/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_full_align.v4.fasta, flip=T, processors=12)
summary.seqs(fasta=16S.trim.contigs.good.unique.align, count=16S.trim.contigs.good.count_table)
#filter out alignments with more than 8 homopolymers
screen.seqs(fasta=16S.trim.contigs.good.unique.align, count=16S.trim.contigs.good.count_table, summary=16S.trim.contigs.good.unique.summary, end= 3004, maxhomop=8)
summary.seqs(fasta=16S.trim.contigs.good.unique.good.align, count=16S.trim.contigs.good.good.count_table)
#filter alignment
filter.seqs(fasta=16S.trim.contigs.good.unique.good.align, vertical=F)
summary.seqs(fasta=16S.trim.contigs.good.unique.good.filter.fasta, count=16S.trim.contigs.good.good.count_table)
unique.seqs(fasta=16S.trim.contigs.good.unique.good.filter.fasta, count=16S.trim.contigs.good.good.count_table)
summary.seqs(fasta=16S.trim.contigs.good.unique.good.filter.unique.fasta, count=16S.trim.contigs.good.unique.good.filter.count_table)
#cluster aligned reads (2bp error threshold)
pre.cluster(fasta=16S.trim.contigs.good.unique.good.filter.unique.fasta, count=16S.trim.contigs.good.unique.good.filter.count_table, diffs=2, processors=2)
summary.seqs(fasta=current, count=current)
#flag chimeras for removal, this step will take a while
chimera.uchime(fasta=16S.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=16S.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t, processors=12)
#filter out chimeras
remove.seqs(fasta=16S.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=16S.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos, count=16S.trim.contigs.good.unique.good.filter.unique.precluster.count_table)
#assign otu taxonomy
classify.seqs(fasta=16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table, template=/home/lymelab/reference_databases/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_full_align.fasta, taxonomy=/home/lymelab/reference_databases/EzBioCloud_16S_database_for_MOTHUR/ezbiocloud_id_taxonomy.tax, cutoff=80, processors=2)
#remove unwanted taxonomic groups 
remove.lineage(fasta=16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table, taxonomy=16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.ezbiocloud_id_taxonomy.wang.taxonomy, taxon=Chloroplast-f__mitochondria-unknown-Archaea-Eukaryota)
summary.seqs(fasta=current, count=current)
#rename representative sequence file
system(mv 16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta 16S.final.fasta)
#rename frequency table file
system(mv 16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.count_table 16S.final.count_table)
#rename taxonomy file
system(mv 16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.ezbiocloud_id_taxonomy.wang.taxonomy 16S.final.taxonomy)
#remove intermediate files
system(ls | grep "16S" | grep -v "final" | while read line; do rm $line; done)
#build reference tree
system(fasttree -nt 16S.final.fasta > 16S.final.tre)

