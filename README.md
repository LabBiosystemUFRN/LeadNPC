#LeadNPC
## Systems Biology-Based Analysis Indicates Global Transcriptional Impairment in Lead-Treated Human Neural Progenitor Cells

In order, to run this pipeline you will need to download the  raw data sequence reads from Sequence Read Archive, accession SRP079342, Gene Expression Omnibus, accession GSE84712, get the Ensembl GRCh38 Human genome reference and annotation (release 91), and have the software Hisat2 and FeatureCounts installed.
##### Create hisat2 index
> hisat2-build -p 8 ${genRefDir}"/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ${genIndexDir}

##### Extract the splice sites
> extract_splice_sites.py ${genRefDir}"/Homo_sapiens.GRCh38.94.gtf" > ${genRefDir}"/Homo_sapiens.GRCh38.94.txt"

##### Generate the aligns
For each sample file do:
> hisat2 -x ${genIndexDir} --known-splicesite-infile ${genRefDir}"/Homo_sapiens.GRCh38.94.txt" -p 8 -1 ${sampleDir}/${file}_1".fastq.gz" -2 ${sampleDir}/${file}_2".fastq.gz"| samtools view -bS - > ${resultDir}/${file}".hisat.bam"; 

##### Count aligned reads 
First, generate the bam files list in a single line
> ls *.hisat.bam | tr '\n' ' '> bamList.txt

Then generate the count's file
> featureCounts -T 32  -t gene -g gene_id -a ${genRefDir}"/Homo_sapiens.GRCh38.94.gtf" -o ./allCountsHisat.txt $(cat bamList.txt)

Now you can use the R scripts as follow:
##### Create logCPM file
> 00ProcessCounts.R

##### PCA analysis and create transcriptogramer objects
> 1analisePCA.R

##### Plot transcriptogramer graphics
> 2plotTrancript.R

##### Cluster superposition analysis
> 4intersecClusters.R

##### Create clusters graphos
This is not a completely automatic process. You will need to use Cytoscape manually...
> 5graphTemposManual.R

##### Generate Figure 2 components
> 6BarraTempos.R

##### Generate Figures 3 and 4
> 8graphoClustersRCy3.R

##### Generate Markers figures
> 14Marquers.R

