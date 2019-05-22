#!/bin/bash

dataDir="/home/clovis/Dropbox/Chumbo/Data/"
sampleDir=${dataDir}"samples/"
sampleList=${dataDir}"sampleList.txt"
genRefDir=${dataDir}"genRef/"
genIndexDir=${genRefDir}"index/"
starResultDir=${dataDir}"alings/star/"
curDir=$(pwd)

if [[ ! -d $sampleDir ]]; then
    mkdir -p $sampleDir
fi

let samples=$(wc -l $sampleList | cut -f 1 -d " ")-1

echo -e "Processando $samples arquivos de amostras...\n"

cd ${sampleDir}
for i in $(tail -$samples $sampleList|cut -f 1 -d ",");do 
 echo $i; 
 #scp -P 4422 cfreis@10.7.40.56:/storages/10TB-79/iarasouza/lead/${i}.* ${sampleDir}
 fastq-dump --gzip --split-files -A ${i}
done;
cd ${curDir}
 #SRR721944${i};done

 
#download da referencia
#wget -c ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
#wget -c ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#gunzip *

#criando os indices
if [[ ! -d $sampleDir ]]; then
    mkdir -p ${genIndexDir}
fi

cd ${genIndexDir}

#hisat buil index
hisat2-build 
#extract splicing sites 4 hisat
extract_splice_sites.py ../Homo_sapiens.GRCh38.93.gtf > ../Homo_sapiens.GRCh38.93.txt
#alinhamento
hisat2 -x ${genIndexDir} --known-splicesite-infile ${genIndexDir}"/Homo_sapiens.GRCh38.93.txt"
-p 2
-U SRR453566_yeast_rnaseq_trimmed.fq.gz
| samtools view -bS - > SRR453566_yeast_rnaseq.hisat.bam

#STAR --runThreadN {number of cores} --runMode genomeGenerate --genomeDir /path/to/resulting/STAR/genome/ --genomeFastaFiles /path/to/genome/fasta/file --sjdbGTFfile /path/to/GTF/or/GFF --sjdbOverhang {read length - 1}
#
    STAR --runThreadN 8 \
        --runMode genomeGenerate \
        --genomeDir /home/clovis/Dropbox/Chumbo/Data/genRef/indexStar \
        --genomeFastaFiles ${genRefDir}"Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
        --sjdbGTFfile ${genRefDir}"Homo_sapiens.GRCh38.93.gtf" \
        --sjdbOverhang 100

for i in $(tail -$samples $sampleList|cut -f 1 -d ",");do 
 echo $i; 
 i="SRR3944315"
    STAR \
    --genomeDir /home/clovis/Dropbox/Chumbo/Data/genRef/indexStar \
    --runThreadN 8 \
    --readFilesIn ${sampleDir}${i}".fastq" \
    --outFileNamePrefix ${starResultDir}${i}".fastq" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes Standard \
    --outSAMunmapped Within \
    --runMode alignReads  genome \
    --quantMode GeneCounts --twopassMode Basic 
    #--readFilesCommand zcat #para arquivos fasq.gz
    
    ;done


STAR --genomeDir /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/ensembl38_STAR_index/ \
--runThreadN 6 \
--readFilesIn Mov10_oe_1.subset.fq \
--outFileNamePrefix ../results/STAR/Mov10_oe_1_ \
--outSAMtype BAM SortedByCoordinate \
 \

