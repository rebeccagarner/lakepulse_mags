# LakePulse MAGs workflow

_Notes on file names_

Metagenomes were labelled with six-digit LakePulse lake identification codes in which the first two digits (**[ECOZONE_NUMBER]**) typically denote the ecozone and the last three digits after the hyphen were randomly generated (e.g., 06-095 is the LakePulse lake ID for Lac Simoncouche in the Boreal Shield ecozone).

Metagenome co-assemblies were labelled with ecozone abbreviations (**[ECOZONE_ALPHA]**; AH: Atlantic Highlands, AM: Atlantic Maritime, BCTC: Boreal/Taiga Cordilleras, BP: Boreal Plains, BS: Boreal Shield, MC: Montane Cordillera, MP: Mixedwood Plains, P: Prairies, PM: Pacific Maritime, SAP: Semi-Arid Plateaux, TP: Taiga Plains). Bins were also labelled with co-assembly ecozone abbreviations.

## Create MAG set

### Trim reads

Trim raw metagenome reads using Trimmomatic v. 0.36.

```sh
for i in *_R1.fastq.gz
do
	base="${i%*_*.fastq.gz}"
	echo "${base}"
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE "${base}"_R1.fastq.gz "${base}"_R2.fastq.gz "${base}"_p_R1.fq.gz "${base}"_u_R1.fq.gz "${base}"_p_R2.fq.gz "${base}"_u_R2.fq.gz ILLUMINACLIP:NovaSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```
Count number of reads in trimmed fastq files.

```sh
for metagenome in PATH/TO/trimmed/[[:digit:]][[:digit:]]-[[:digit:]][[:digit:]][[:digit:]]*_p_R1.fq.gz
do
	sampleid=`basename $metagenome _p_R1.fq.gz`
	zcat PATH/TO/trimmed/${sampleid}_p_R1.fq.gz | echo $((`wc -l`/4)) > PATH/TO/lpmetagenomes_nreads/nreads_${sampleid}_p_R1.txt
done
```

### Co-assemble metagenomes

Generate metagenome co-assemblies by ecozone from paired trimmed reads using MEGAHIT v. 1.2.7.

```sh
megahit -t 40 -1 [METAGENOME1]_p_R1.fq.gz,[METAGENOME2]_p_R1.fq.gz,[METAGENOME3]_p_R1.fq.gz,[...],[METAGENOMEn]_p_R1.fq.gz -2 [METAGENOME1]_p_R2.fq.gz,[METAGENOME2]_p_R2.fq.gz,[METAGENOME3]_p_R2.fq.gz,[...],[METAGENOMEn]_p_R2.fq.gz -o ~/scratch/[ECOZONE_ALPHA] --k-list 27,37,47,57,67,77,87 --min-count 2
```

If MEGAHIT run is interrupted, continue with:

```sh
megahit --continue -o ~/scratch/[ECOZONE_ALPHA] 
```

### Bin contigs

#### Calculate coverage

Index co-assemblies using BWA 0.7.17.

```sh
for i in [ECOZONE_ALPHA].final.contigs.fa; do
	variable_name=$(echo "$i" | sed -e 's/.final.contigs.fa//g')
	bwa index -p "$variable_name"-index -a bwtsw $i; wait
	pigz $i
done
```

Calculate depth of coverage using BWA 0.7.17.

```sh
for i in *_p_R1.fq.gz
do
	variable_name=$(echo $i | cut -f 1 -d '_')
	echo $variable_name;
	bwa mem -t 40 ~/scratch/[ECOZONE_ALPHA]/[ECOZONE_ALPHA]-index ~/scratch/LP2019/trimmed/"$variable_name"_p_R1.fq.gz ~/scratch/LP2019/trimmed/"$variable_name"_p_R2.fq.gz | samtools view -h -o ~/scratch/coassemblies/mapped/"$variable_name".sam; wait
	samtools flagstat ~/scratch/coassemblies/mapped/"$variable_name".sam; wait
	samtools view -F 4 -bS ~/scratch/coassemblies/mapped/"$variable_name".sam | samtools sort -o ~/scratch/coassemblies/mapped/"$variable_name".bam; wait
	rm ~/scratch/coassemblies/mapped/"$variable_name".sam
done
```

Summarize depth of coverage using MetaBAT v. 2.12.1.

```sh
jgi_summarize_bam_contig_depths --outputDepth [ECOZONE_ALPHA]_depth.txt [ECOZONE_NUMBER]-*.bam
```

#### Perform binning

Bin ecozone co-assembly contigs based on depths of coverage using MetaBAT v. 2.12.1.

```sh
metabat -i [ECOZONE_ALPHA].final.contigs.fa -a [ECOZONE_ALPHA]_depth.txt -o [ECOZONE_ALPHA]_bin
```

### Filter bins

Filter bins based on quality assessments generated using CheckM v. 1.0.7.

```sh
# Information on genome completeness, contamination, strain heterogeneity
checkm lineage_wf PATH/TO/lpmags/ PATH/TO/lpmags_checkm/ -x fa -t 8

# Phylogenetic markers found in each bin
checkm tree_qa -o 2 -f lpmags_checkm_tree_qa.txt --tab_table PATH/TO/lpmags_checkm/
```

### Dereplicate MAGs

Dereplicate MAGs at 95% ANI using dRep v. 3.2.0.

```sh
dRep dereplicate -g lpmags_paths.txt -comp 50 --S_algorithm gANI -sa 0.95 -p 8 PATH/TO/lpmags_drep/
```

## Annotate MAGs

### Annotate genome features

Annotate MAG genome features with Prokka v. 1.12.

```sh
for filename in *.fa
	do
	mag=`basename $filename .fa`
	prokka --outdir ${mag}_prokka --prefix ${mag} $filename
done
```

### Annotate gene functions

Annotate MAG gene functions with KofamScan.

```sh
for file in *.faa
do
	exec_annotation --tmp-dir $file.tmp -f detail-tsv --report-unannotated -o $file.kofam.tsv $file
done
```

Filter KO table.

```python
#!/usr/bin/env python2.7

import sys
import os
import argparse
import pandas as pd
pd.options.mode.chained_assignment = None 

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',action="store",help='KOFAMSCAN_output after removal of first two symbols, include path',dest="input")
parser.add_argument('-t', '--threshold',action="store",help='Threshold: between 0 and 1',dest="threshold")
parser.add_argument('-d', '--dash',default=False,help='Include gene families without threshold: True or False',dest='dash')
args = parser.parse_args()

#print(args)

path=args.input
T_set=args.threshold
data = pd.read_csv(path, sep="\t")
data1=data.iloc[1:]
data2=data1.drop(data1.columns[0],axis=1)
dataframe = data2.drop_duplicates(subset=['gene name'])
dataframe_dash = dataframe.loc[dataframe['thrshld']=='-']
dataframe_dash['T-value'] = '0'
dataframe_clean = dataframe.loc[dataframe['thrshld']!='-']
dataframe_clean.loc[:,'score']=dataframe_clean.loc[:,'score'].apply(lambda x: float(x))
dataframe_clean.loc[:,'thrshld']=dataframe_clean.loc[:,'thrshld'].apply(lambda x: float(x))
dataframe_clean.loc[:,'T-value'] = dataframe_clean.loc[:,'score']/dataframe_clean.loc[:,'thrshld']
filt=dataframe_clean.loc[dataframe_clean['T-value']>=float(T_set)]
length=filt.shape[0]
if args.dash == 'True':
	#print("dash")
	frames=[filt,dataframe_dash]
	df_comp=pd.concat(frames)
	df_comp.to_csv(str(path)+"_filter_"+str(T_set)+".txt",header=True,index=False,sep='\t',mode='a')
else:
	filt.to_csv(str(path)+"_filter_"+str(T_set)+".txt",header=True,index=False,sep='\t',mode='a')
print(str(length)+" genes pass the filter criteria at threshold "+str(T_set))
```

### Annotate RNA genes

Annotate MAG RNA genes using Infernal v. 1.1.2.

First, compress the Rfam database.

```sh
cmpress Rfam.cm
```

Next, run Infernal cmscan.

```sh
for filename in *.fa
do
	mag=`basename $filename .fa`
	cmscan --rfam --cut_ga --nohmmonly --tblout ${mag}_cmscan_rfam.tblout --fmt 2 --clanin Rfam.clanin --cpu 12 Rfam.cm $filename > ${mag}_cmscan_rfam.out
done
```

### Annotate transporters

Make BLAST database for the Transporter Classification Database (TCDB; release 2021/01/27).

```sh
makeblastdb -in tcdb_20210127.faa -dbtype prot -out tcdb_20210127_blastdb
```

Perform protein BLASTs against the TCDB.

```sh
for folder in PATH/TO/*prokka/
do
	cd $folder
	mag=*.faa
	prefix=`basename $mag _prokka.faa`
	blastp -query $mag -db tcdb_20210127_blastdb -out PATH/TO/lpmags_tcdb/${prefix}_tcdb.txt -evalue 0.00000000000000000001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"
done
```

### Annotate CAZymes

Perform carbohydrate-active enzyme (CAZyme) HMM searches with HMMER 3.1b2.

```sh
for folder in PATH/TO/*_prokka/
do
	cd $folder
	mag=*.faa
	prefix=`basename $mag _prokka.faa`
	hmmsearch -o PATH/TO/lpmags_cazymes/${prefix}_cazymehmm.out --tblout PATH/TO/lpmags_cazymes/${prefix}_cazymehmm_tblout.txt --domtblout PATH/TO/lpmags_cazymes/${prefix}_cazymehmm_domtblout.txt -E 0.000000000000001 --cpu 12 dbCAN-HMMdb-V9.txt $mag
done
```

## Classify MAGs

Run GTDB-Tk v. 1.3.0.

```sh
# Step 1. Identify (marker genes)
gtdbtk identify --batchfile lpmags_batchfile_allbins.txt --prefix lpmags --cpus 12 --out_dir PATH/TO/lpmags_gtdbtk_identify

# Step 2. Align
gtdbtk align --identify_dir PATH/TO/lpmags_gtdbtk_identify/ --prefix lpmags --cpus 12 --out_dir PATH/TO/lpmags_gtdbtk_align/

# Step 3. Infer (phylogenetic tree)
gtdbtk infer --msa_file PATH/TO/lpmags_gtdbtk_align/lpmags.bac120.user_msa.fasta --prot_model JTT --prefix lpmags --cpus 12 --out_dir PATH/TO/lpmags_gtdbtk_infer/

# Step 4. Root (tree)
# Unnecessary to root tree

# Step 5. Classify
gtdbtk classify --batchfile lpmags_batchfile_allbins.txt --align_dir PATH/TO/lpmags_gtdbtk_align/ --prefix lpmags --cpus 8 --out_dir PATH/TO/lpmags_gtdbtk_classify/
```

## Perform fragment recruitment

### Mask RNA genes

#### Write rRNA and tRNA gene coordinates

Write ribosomal and transfer RNA gene coordinates in BED format (refer to R script/s).

#### Mask rRNA and tRNA genes

Mask ribosomal and transfer RNA genes (with *N* nucleotides) using Bedtools v. 2.26.0 maskfasta.

```sh
for filename in *_rna.bed
do
	prefix=`basename $filename _rna.bed`
	bedtools maskfasta -fi ${prefix}.fa -bed $filename -fo ${prefix}_maskrna.fasta
done
```

### Map unassembled reads to MAGs

#### Index MAGs

Index MAGs (with masked ribosomal and transfer RNA genes) using BBMap v. 37.36

```sh
for mag in *_maskrna.fasta
do
	prefix=`basename $mag _maskrna.fasta`
	bbmap.sh t=40 ref=$mag path=PATH/TO/${prefix}_ref
done
```

#### Map reads

Map unassembled metagenome reads to MAGs (with masked ribosomal and transfer RNA genes) at 95% sequence identity threshold using BBMap v. 37.36 and format mapping files with Samtools v. 1.9.

```sh
cd PATH/TO/lpmags_ref/
for folder in [ECOZONE_ALPHA]*.[[:digit:]][[:digit:]][[:digit:]]_ref
do
	for metagenome in PATH/TO/lpmags_trimmed/[ECOZONE_NUMBER]-[[:digit:]][[:digit:]][[:digit:]]_p_R1.fq.gz
	do
		sampleid=`basename $metagenome _p_R1.fq.gz`
		
		# cd into _ref directory
		cd PATH/TO/lpmags_ref/${folder}
		
		# Define prefix by removing _ref from directory names
		prefix=`basename $folder _ref`
		
		# Recruit reads with BBMap
		bbmap.sh t=40 maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t usejni=t minid=0.95 covstats=PATH/TO/lpmags_bbmap/${prefix}_recruit_${sampleid}_reads_95id_covstats.txt idhist=PATH/TO/lpmags_bbmap/${prefix}_recruit_${sampleid}_reads_95id_idhist.txt in=PATH/TO/lpmags_trimmed/${sampleid}_p_R1.fq.gz in2=PATH/TO/lpmags_trimmed/${sampleid}_p_R2.fq.gz outm=PATH/TO/lpmags_bbmap/${prefix}_recruit_${sampleid}_reads_95id.sam
		cd PATH/TO/lpmags_bbmap/
		
		# Convert SAM to BAM
		samtools view -S -b ${prefix}_recruit_${sampleid}_reads_95id.sam > ${prefix}_recruit_${sampleid}_reads_95id.bam
		
		# Sort BAM
		samtools sort ${prefix}_recruit_${sampleid}_reads_95id.bam -o ${prefix}_recruit_${sampleid}_reads_95id_sorted.bam
		
		# Delete unsorted BAM
		rm ${prefix}_recruit_${sampleid}_reads_95id.bam
		
		# Zip files
		pigz ${prefix}_recruit_${sampleid}_reads_95id*
	done
done
```

### Format and count mapped reads

#### Retain mapped reads

Retain only mapped reads in sorted BAM files.

```sh
cd PATH/TO/bam/
while read p
do
	filename="$p"
	prefix=`basename $filename _sorted.bam`
	samtools view -b -F 4 $filename > PATH/TO/bam_mappedonly/${prefix}_mappedonly.bam
done <PATH/TO/bbmap_server_unzipped.txt
```

#### Filter mapping files

Filter mapping files at a strict 96% sequence identity threshold.

```python
#!/usr/local/bin/python3.6

import pysam
import sys

threshold = 0.96

samfile = pysam.AlignmentFile(sys.argv[1], "rb")
with pysam.AlignmentFile(sys.argv[1] + ".filtered", "wb", header=samfile.header) as outf:
	for aligned_segment in samfile.fetch():
		#print(aligned_segment)
		cgs = aligned_segment.get_cigar_stats()
		# first array contains the number of operations for each CIGAR op
		# see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats\
		match_len= cgs[0][7]
		mismatch_len=cgs[0][8]
		tot_len=match_len+mismatch_len
		if tot_len>0:
			identity = match_len/tot_len 
			if identity >= threshold:
				outf.write(aligned_segment)
```

#### Count mapped reads

Count number of filtered mapped reads.

```sh
# Create a separate file for each BAM containing the number of reads
cd PATH/TO/bam_mappedonly/
while read p
do
	filename="$p"
	prefix=`basename $filename _mappedonly.bam.filtered`
	samtools view -c $filename > PATH/TO/lpmags_nreads/${prefix}_mappedonly_filtered96_nreads.txt
done <PATH/TO/lpmags_mappedonly_filtered_bam_ls.txt

# Combine nreads files into a single file where column 1 is the file name and column 2 is the nreads
while read p
do
	awk '{nreads+=$1} END {print FILENAME"\t"nreads}' "$p" >> lpmags_bam_filtered96_nreads_all.txt
done <lpmags_nreads_ls.txt
```

## Create TAD80 matrix

### Count genome equivalents

Count genome equivalents in metagenomes using MicrobeCensus v. 1.1.0.

### Create TAD80 table

Produce BedGraph file with bedtools genomecov (considering zero-coverage positions with -bga flag).

```sh
for i in *.filtered
do
	if test -f ../bed_out/$i.txt
	then
		echo $i
	else
		bedtools genomecov -bga -ibam $i >> ../bed_out/$i.txt
	fi
done
```

Estimate the central 80% truncated average sequencing depth (TAD80) from BedGraph file.

```sh
for i in *.txt
do
	base="${i%.*}"
	if test -f ../tad_out/$base.txt
	then
		echo $base
	else
		BedGraph.tad.rb -i $i -r 0.8 >> ../tad_out/$base.txt
	fi
done
```

### Calculate MAG size

Count number of ACTG nucleotides (excludes Ns) in RNA gene-masked MAG fasta files.

```sh
for mag in *_maskrna.fasta
do
	prefix=`basename $mag _maskrna.fasta`
	grep -o '[ACTG]' $mag | wc -l | tr -d '[:space:]' > PATH/TO/lpmags_maskrna_size_bp/${prefix}_size_bp.txt
done
```

Count number of ACTG nucleotides (excludes Ns) in MAG fasta files

```sh
for mag in *.fa
do
	prefix=`basename $mag .fa`
	grep -o '[ACTG]' $mag | wc -l | tr -d '[:space:]' > PATH/TO/lpmags_size_bp/${prefix}_size_bp.txt
done
```
