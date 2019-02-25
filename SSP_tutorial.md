---


---

<h2 id="protocol-1-small-secret-peptide-gene-discovery-from-genomic-sequences">Protocol #1: Small Secret Peptide Gene discovery from genomic sequences</h2>
<p>Traditional genome annotation policy is biased to discover long gene; leading to missing of some small secret peptide (SSP) genes. The following workflow was optimized to identify SSP gene from assembled genomic sequences utilizing specific RNA-seq data as expression evidence and conserved SSP motif.</p>
<h3 id="prerequisite">1.1 Prerequisite</h3>
<p><strong>1.1.1 Suggestion:</strong></p>
<ol>
<li>We recommend an X-windows desktop (such as gnome/XFCE/MATE) stead of SSH terminal because it is more convenient to edit file.</li>
<li>All commands below are typed  under Linux terminal.</li>
<li>A line start with <code>#</code> in Linux command line indicates that this is explanatory information only.</li>
<li>You need a no-root user with sudo privilege in host system to install docker packages and enable docker service. Asking your system administrator to install docker service and add your user as a member of docker group if you can’t have <code>sudo</code> privileges.</li>
<li>You need to be a sudo user or a member of <code>docker</code> group in host system to start Docker container and attach to container terminal.</li>
<li>Default user in Docker container is <code>test</code>.</li>
</ol>
<p><strong>1.1.2 Computer:</strong><br>
A  high-performance computer (I7/Xeron processor and &gt;16GB RAM) with CentOS 7, Ubuntu 16.04 or higher as operation system.</p>
<p><strong>1.1.3 Work folder</strong><br>
Work folder is the place for raw input data (genomic sequences, gff, RNA-seq and protein, etc ) and analysis result. It usually can be <code>work</code> under current directory (your home directory).</p>
<p><strong>1.1.4 Input data</strong></p>
<ul>
<li>Genomics sequences in FASTA format</li>
<li>Reference annotation in gff format if available</li>
<li>SSP gene expression specific RNA-seq data in compress fastq format</li>
<li>Protein sequence of known SSP genes other related protein sequences in FASTA format</li>
<li>Other EST/transcript sequences from the same species.</li>
</ul>
<p><strong>1.1.5 Demo data</strong><br>
The demo data is available for <a href="">download</a>. Copy it to <code>work</code> folder and type the following command to unzip it.</p>
<pre><code>cd work
tar -xzvf ssp-demo.tar.gz
</code></pre>
<p>The above command will generate <code>ssp</code> folder under <code>work</code></p>
<p><strong>1.1.6 Software installation:</strong><br>
All software have been configured and packed as a docker image hosted in <a href="https://hub.docker.com/">Docker Hub</a>. Firstly, install docker packages and enable/start docker service:</p>
<pre><code>#under CentOS 7, install docker packages:
sudo yum install docker
#or if you are using Ubuntu, install docker packages as below:
sudo apt install docker.io

#enable and start docker service for CentOS/ubuntu
sudo systemctl enable docker
sudo systemctl start docker    
</code></pre>
<p>Then, start a container of SSP-mining image to input Linux command line:</p>
<pre><code>#start a Docker container with name `bioinfo` using xxx as template image
#this step will take a while depend on your network speed
sudo docker run -d -it -e "uid=$(id -u)" -e "gid=$(id -g)" --name bioinfo -v $(pwd)/work:/work centos7bioperl:test6 shell
sudo docker attach bioinfo
</code></pre>
<p><code>$(pwd)/work</code> indicates that the <code>work</code> folder under current directory will be mounted on <code>/work</code> in Docker container. This is the folder you exchange data between Host computer (<code>$(pwd)/work</code>) and container virtual machine (<code>/work</code>).</p>
<p>The <code>attach</code> subcommand will link your current Linux terminal to the running docker container (<code>bioinfo</code> in this case).  Tips:  Hold Ctrl key and press P,Q to detach the container terminal and get back to host OS.</p>
<p>Type the following command to enter demo data folder work folder in attached Docker container terminal</p>
<pre><code>cd /work/ssp
</code></pre>
<p>All Linux commands below should be typed in this terminal.</p>
<h3 id="optimized-genome-annotation-procedure-for-mining-ssp-gene-using-maker-pipeline">1.2 Optimized genome annotation procedure for mining SSP gene using Maker pipeline</h3>
<p>General genome annotation procedure can be optimized to identify more SSP gene through including SSP-specific expression evidence and conserved known SSP domain.</p>
<h4 id="preparing-rna-sed-based-gene-expression-evidence-for-maker-pipeline">1.2.1 Preparing RNA-sed based gene expression evidence for Maker pipeline</h4>
<p>Some plant SSP genes may only express under limited condition or tissue, such as nutrient deficient root. Related RNA-seq data will help to improve the performance of SSP gene mining. The following sample code will perform  reference-based transcriptome assembly and generate a gff file for downstream Maker analysis.</p>
<pre><code>cd /work/ssp
mkdir transcriptome
cd transcriptome/
gffread data/ref.gff3 -T -o ref.gtf
python /opt/hisat2/hisat2_extract_splice_sites.py ref.gtf &gt; splicesites.txt
time /opt/hisat2/hisat2 -p 20 -x genome_hisat2 --known-splicesite-infile splicesites.txt --dta --dta-cufflinks -1 data/RNA-seq/root_R1.fq.gz,data/RNA-seq/nod_R1.fq.gz,data/RNA-seq/bud_R1.fq.gz -2 data/RNA-seq/root_R2.fq.gz,data/RNA-seq/nod_R2.fq.gz,data/RNA-seq/bud_R2.fq.gz -U data/RNA-seq/SRR1377073.fastq.gz,data/RNA-seq/SRR1377076.fastq.gz | samtools view -bS - &gt; all_runs.bam
sambamba sort -m 40G --tmpdir tmp/ -o all_runs.sorted.bam -p -t 20 all_runs.bam
stringtie all_runs.sorted.bam -o transcriptome_models.gtf -p 20
cufflinks2gff3  transcriptome/transcriptome_models.gtf &gt; transcriptome/transcriptome_models.gff3
</code></pre>
<p><code>-p 20</code> or <code>-t 20</code> is the number of CPU cores assigned to the program. Type <code>nproc</code> to check the maximum number in your computer.<br>
<code>-m 40G</code> is max RAM size assigned to your computer. Type <code>free</code> to check your computer RAM size.<br>
<code>transcriptome/transcriptome_models.gff3</code> contains transcriptome data which is expression evidence in Maker genome annotation (next step)</p>
<h4 id="identification-of-ssp-gene-using-maker-pipeline">1.2.2  Identification of SSP gene using Maker pipeline</h4>
<p>The protocol for genome annotation using Maker has been well documented [ref]. We installed and tested Maker pipeline in the Docker image xxx. Users need to generate three maker_opts files: <code>maker_opts_1.ctl</code>, <code>maker_opts_2.ctl</code> and <code>maker_opts_3.ctl</code>. In addition, Maker also need maker_bopts.ctl and maker_exe.ctl files. These files include input data files path and other settings for genome annotation. Maker will take these files  as input to generate final gff file with genome annotation information. The annotation procedure will be typically invoked for three rounds to generate optimized result.</p>
<p>The gff file for transcriptome generated by previous step and known SSP protein sequences (as of 01/2019, under /work/ssp/data) will be put in above three maker_opts files (line 18 and 23). The extra information will help maker to identify novel SSP genes.</p>
<p>We included all Maker configuration files in demo data. Users should be able to test Maker for the first round using the following command:</p>
<pre><code>cd /work/ssp
time /usr/lib64/mpich/bin/mpiexec -n 20 maker -fix_nucleotides maker_opts_1.ctl maker_bopts.ctl maker_exe.ctl  1&gt;&amp;2 2&gt;log
</code></pre>
<p><code>-n 20</code> is the number of CPU cores. Please check <code>log</code> file if pipeline failed.</p>
<p>To generate gff file, run this command:</p>
<pre><code>gff3_merge -d genome.maker.output/genome_master_datastore_index.log -s &gt; maker_all.gff
(head -1 maker_all.gff;cat maker_heavy.gff|grep -P $'\t(CDS|contig|exon|five_prime_UTR|gene|mRNA|three_prime_UTR)\t') &gt; maker.gff
</code></pre>
<h3 id="ssp-genome-annotation-using-spada-pipeline">1.3 SSP genome annotation using SPADA pipeline</h3>
<p>SPADA pipeline typically utilizes conserved domain (in HMM format) information of known SSP families to identify SSP gene from genomic sequences[ref]. We included a copy of SPADA and a comprehensive HMM dataset from <a href="">PlantSSPDB</a> database and our curation in Docker image.</p>
<p>Here is an example of SPADA analysis:</p>
<pre><code>perl /opt/spada_soft/spada/spada.pl --cfg /opt/spada_soft/spada/cfg.txt -d sspanno -p /opt/spada_soft/spada/CRP_PlantSSPv1_Noble -f data/genome.fa -t 20 -o arabidopsis 1&gt;&amp;2 2&gt;log
</code></pre>
<p>In this example,<br>
<code>/opt/spada_soft/spada/</code>: SPADA pipeline is installed here.<br>
<code>/opt/spada_soft/spada/CRP_PlantSSPv1_Noble</code>: location of above-mentioned comprehensive HMM dataset. You are welcome to change it to any your favorite HMM dataset.<br>
<code>data/genome.fa</code>: is the genomic sequences to be analyzed.<br>
<code>arabidopsis</code>: one of folder name under <code>/opt/spada_soft/augustus/config/species/</code>. Change it to the closest species in your case.</p>
<p>The gff  file is available at <code>sspanno/31_model_evaluation/61_final.gff</code></p>
<h3 id="merge-the-annotation-result-from-maker-and-spada">1.4 Merge the annotation result from Maker and SPADA</h3>
<p>There are (partially) duplicate genes between Maker and SPADA output. Additionally, users also need to consider the duplicate gene model in reference annotation if you are working on a model plant. Here we describe the policy of  identifying duplicated gene list:</p>
<ul>
<li>
<p>generate cds sequence from different annotation using gff and genomic sequences as described in xxx</p>
<pre><code> #generate CDS sequences for SPADA gene models
 gffread sspanno/31_model_evaluation/61_final.gff -g data/genome.fa -x spada_cds.fa
 #replace transcript id with gene id
 sed -ri 's/^&gt;\S+\s+gene=/&gt;/' spada_cds.fa     
</code></pre>
</li>
<li>
<p>Run NCBI BLASTN between two generated cds files and only keep the query-hit gene pairs with &gt; 50% overlapped region.</p>
</li>
<li>
<p>Check the coordinates of the query-hit pairs in gff file. If both genes share overlapped region and the same direction on chromosome, choose one of them as redundant gene.</p>
</li>
<li>
<p>Remove duplicate genes from corresponding gff file.  <code>nr_gene_in_spada.txt</code> includes gen list to be removed. <code>new_spada.gff</code> is new gff file.</p>
<pre><code>  gffremove.py --infile sspanno/31_model_evaluation/61_final.gff --outfile new_spada.gff --genefile nr_gene_in_spada.txt
</code></pre>
</li>
<li>
<p>Merge two gff files (taking <code>maker.gff</code> and <code>new_spada.gff</code> as example):</p>
<pre><code>  head -1 maker.gff  &gt; all.gff
  grep -v '^#' maker.gff new_spada.gff  &gt;&gt; maker_spada.gff
</code></pre>
</li>
<li>
<p>Generate protein and transcript files:</p>
<pre><code> gffread all.gff -g data/genome.fa -y all_protein.fa
 sed -ri 's/^&gt;\S+\s+gene=/&gt;/' all_protein.fa
 gffread all.gff -g data/genome.fa -w all_transcript.fa
</code></pre>
<p>The protein and transcript files will be used to further annotate gene function in protocol #2.</p>
</li>
</ul>
<h2 id="protocol-2-functional-annotation-and-family-classification-of-ssp-genes">Protocol #2: Functional annotation and family classification of SSP genes</h2>
<p>Due to the short conserved region and less homologous among the member of the same gene family, it is less efficient to search and identify SSP protein using standard NCBI BLASTP.  Here we introduce a comprehensive  annotation procedure to identify SSP gene from candidate genes.</p>
<h3 id="prerequisite-1">2.1 Prerequisite</h3>
<p>The Prerequisite is almost the same as protocol #1 (see 1.1) with only exception on input data.</p>
<p><strong>2.1.1 Input data</strong></p>
<ul>
<li>SSP gene expression specific RNA-seq data in compress fastq format</li>
<li>Protein sequence (with gene ID as protein ID) of SSP gene candidates in FASTA format</li>
<li>Transcript sequences of SSP gene  candidates and a two-columns mapping file between gene and transcript id.</li>
</ul>
<h3 id="smith-waterman-search-against-known-ssp-protein">2.2 Smith-waterman search against known SSP protein</h3>
<p>Smith-Waterman alignment is more accurate way for sequence homolog search compared to BLAST search. The wrapped shell script <code>swsearch</code> use <a href="https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml">FASTA</a> software to perform Smith-waterman search against known SSP proteins (<code>target.fa</code> in demo data) and take top 2 hits (Expection&lt;0.01) as output. The script will remove hit sequence ID and only output family name.</p>
<pre><code>swsearch  all_protein.fa /work/ssp/data/ssp_family.fa &gt; sw.txt
</code></pre>
<h3 id="hmm-search-against-hmms-of-known-ssp-families">2.3 HMM search against HMMs of known SSP families</h3>
<p>Generate HMM library for all known SSP families from SPADA pipeline</p>
<pre><code>cat /opt/spada_soft/spada/CRP_PlantSSPv1_Noble/15_hmm/*.hmm &gt; all.hmm
</code></pre>
<p>compile HMM library</p>
<pre><code>hmmpress  all.hmm 
</code></pre>
<p>search your protein sequence ‘all_protein.fa’ against HMM library ‘all.hmm’</p>
<pre><code>hmmscan --cpu 20 -E 0.01 --tblout hmm_output.txt all.hmm all_protein.fa &gt; /dev/null
</code></pre>
<p>The expectation cutoff for <code>hmmscan</code> is <code>0.01</code>. <code>all_protein.fa</code> is input protein file. The above commands will generate a tab-delimited table file <code>hmm_output.txt</code>.  In the result table, column 1 is HMM family name and column 2 is gene ID in input protein sequence file.</p>
<h3 id="signal-peptide-detection">2.4 Signal peptide detection</h3>
<p>We choose SignalP to predict signal peptide from SSP candidate peptides. We recommend to use “No TM” and long output format. D-score thresholds are usually 0.45 or 0.5, depending on the type of network chosen (with or without transmembrane segments), but we recommend to use a D-score of ≥ 0.45 for known SSPs or ≥ 0.25 for putative SSPs.</p>
<pre><code>/opt/spada_soft/signalp-4.1/signalp -t euk -f long -s notm all_protein.fa &gt; signalp_long.txt
cat signalp_long.txt | singalP_parser  &gt; sp.txt
</code></pre>
<p>In this case, <code>/opt/spada_soft/signalp-4.1/signalp</code> is Signalp program which output prediction result file <code>signalp_long.txt</code>.  <code>singalP_parser</code> is a script to parse long format output of SignalP. The final signal peptide prediction result <code>sp.txt</code> includes four columns:  gene ID,  start coordinates,  end coordinates,  D-value, cut-off and conclusion(YES/NO)</p>
<h3 id="identification-of-novel-ssp-gene-families-using-mcl-analysis">2.5 Identification of novel SSP gene families using MCL analysis</h3>
<p>SSP candidate (signalp D-value&gt;0.25)   can be clustered into candidate SSP families using Markov Chain Cluster (MCL). The procedure should be performed on last 50 a.a. and candidate  pipetide should be less than 200 a.a…</p>
<pre><code>#create index to retrieve sequence by sequnece ID
cdbfasta all_protein.fa
#only will take the protein with D-value&gt;0.25
cat sp.txt | awk '{if($4&gt;0.25) print $1}' | cdbyank all_protein.fa.cidx  &gt; all_putative_ssp.fa
# the final candidates should be shorter than 200 a.a. and only take the last 50 a.a..
shortseqtail all_putative_ssp.fa 200 50 &gt; peptide-tail.fa
</code></pre>
<p>Generate protein vs protein relationship pairs:</p>
<pre><code>search35_t -T 20 -Q -H -m 9 -b 100 -d 100 peptide-tail.fa peptide-tail.fa &gt; sw_for_peptidetail
bioparser -t ssearch -m sw_for_peptidetail |awk '{print $3,$6,$14}' FS="\t" |sort | uniq &gt; protein-protein-rel.txt
</code></pre>
<p><code>protein-protein-rel.txt</code> is a three column file with two protein/gene IDs and e-value.</p>
<h3 id="gene-expression-analysis-of-rna-seq-data">2.6 Gene expression analysis of RNA-seq data</h3>
<h3 id="comprehensive-table-of-gene-annotation-evidence">2.7 Comprehensive table of gene annotation evidence</h3>
<p>The eviden</p>

