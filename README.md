---


---

<h2 id="protocol-1-small-secret-peptide-gene-discovery-from-genomic-sequences">Protocol #1: Small Secret Peptide Gene discovery from genomic sequences</h2>
<p>Traditional genome annotation policy is biased to discover long gene; leading to missing of some small secret peptide (SSP) genes. The following workflow was optimized to identify SSP gene from assembled genomic sequences utilizing specific RNA-seq data as expression evidence and conserved SSP motif.</p>
<h3 id="prerequisite">1.1 Prerequisite</h3>
<p><strong>1.1.1 Suggestion:</strong></p>
<ol>
<li>We recommend an X-windows desktop (such as gnome/XFCE/MATE) instead of SSH terminal because it is more convenient to edit file.</li>
<li>All commands below are typed  under Linux terminal.</li>
<li>A line start with <code>#</code> in Linux command line indicates that this is explanatory information only.</li>
<li>You need a no-root user with sudo privilege in host system to install docker packages and enable docker service. Asking your system administrator to install docker service and add your user as a member of docker group if you can’t have <code>sudo</code> privileges.</li>
<li>You need to be a sudo user or a member of <code>docker</code> group in host system to start Docker container and attach to container terminal.</li>
<li>Default user in Docker container is <code>test</code>.</li>
</ol>
<p><strong>1.1.2 Computer:</strong><br>
A  high-performance computer (I7/Xeon processor and &gt;16GB RAM) with CentOS 7, Ubuntu 16.04 or higher as your host operation system(OS).</p>
<p><strong>1.1.3 Work folder</strong><br>
Work folder is the place for raw input data (genomic sequences, gff, RNA-seq and protein, etc ) and analysis result in your host OS. It is recommend to be <code>work</code> under your home directory.  For example, if your username is <code>test</code> in host OS, the recommended work folder will be <code>/home/test/work</code> in your host OS. To make work in your home directory of host OS:</p>
<pre><code>cd ~   # ~ means your home directory, e.g. /home/test
mkdir work
</code></pre>
<p><strong>1.1.4 Input data</strong></p>
<ul>
<li>Genomics sequences in FASTA format</li>
<li>Reference annotation in gff format if available</li>
<li>SSP gene expression specific RNA-seq data in compress fastq format</li>
<li>Protein sequence of known SSP genes other related protein sequences in FASTA format</li>
<li>Other EST/transcript sequences from the same species.</li>
</ul>
<p><strong>1.1.5 Demo data</strong><br>
The demo data is available for <a href="">download</a>. <strong>In host OS</strong>, copy it to your work folder and type the following command to unzip it.</p>
<pre><code>cd ~/work
tar -xzvf ssp-demo.tar.gz
</code></pre>
<p>The above command will generate <code>ssp</code> folder under <code>work</code>.</p>
<p>In the <code>~/work/ssp/data</code> folder, <code>ssp_family.fa</code> is a protein sequences of known SSP genes. The known SSP file is used in Maker genome annotation (Protocol #1) and SSP gene annotation (Protocol #2).</p>
<p><strong>1.1.6 Software installation:</strong><br>
All software have been configured and packed as a docker image hosted in <a href="https://hub.docker.com/">Docker Hub</a>. Firstly, install docker packages and enable/start docker service in your host OS:</p>
<p>Under CentOS 7, install docker packages:</p>
<pre><code>sudo yum install docker
</code></pre>
<p>If you are using Ubuntu, install docker packages as below:</p>
<pre><code>sudo apt install docker.io
</code></pre>
<p>Enable and start docker service for CentOS/ubuntu:</p>
<pre><code>sudo systemctl enable docker
sudo systemctl start docker    
</code></pre>
<p>Then, start a container of SSP-mining image to input Linux command line:</p>
<pre><code>sudo docker run -d -it -e "uid=$(id -u)" -e "gid=$(id -g)" --name sspvm -v $(pwd)/work:/work docker.io/noblebioinfo/sspgene
sudo docker attach sspvm
</code></pre>
<p>The above commands will start a Docker container named <code>sspvm</code> using  <code>docker.io/noblebioinfo/sspgene</code> as template image. The step will take a while depend on your network download speed.</p>
<p>In  <code>-v $(pwd)/work:/work</code>:  <code>$(pwd)/work</code>, the path of work folder in your host OS, is <code>work</code> under your current directory. Here,<code>$(pwd)</code> will be converted to your current folder, e.g. home folder by Linux Bash interpreter. The work folder in host OS will be mounted on <code>/work</code> in Docker container. Thus, the folder makes it possible to exchange data between Host computer (<code>$(pwd)/work</code>) and Docker container  (<code>/work</code>). You can copy your demo data or other research data to the work folder in hosts OS (<code>$(pwd)/work</code>) and access them in <code>/work</code> in Docker container.</p>
<p>The <code>attach</code> subcommand will link your current Linux terminal to the running docker container (<code>bioinfo</code> in this case).  Tips:  Hold Ctrl key and press P,Q to detach the container terminal and get back to host OS.</p>
<p>Type the following command to enter demo data folder work folder in attached Docker container terminal</p>
<pre><code>cd /work/ssp
</code></pre>
<p>All Linux commands below should be typed in this container terminal.</p>
<h3 id="preparing-rna-sed-based-gene-expression-evidence-for-maker-pipeline">1.2 Preparing RNA-sed based gene expression evidence for Maker pipeline</h3>
<p>Some plant SSP genes may only express under limited condition or tissue, such as nutrient deficient root. Related RNA-seq data will help to improve the performance of SSP gene mining. The following sample code will perform  reference-based transcriptome assembly and generate a GFF file for Maker genome annotation.</p>
<p><strong>1.2.1 Prepare work folder</strong></p>
<pre><code>cd /work/ssp
mkdir transcriptome
cd transcriptome/
</code></pre>
<p><strong>1.2.2 Compile the genomics sequences using HISAT2</strong></p>
<pre><code>hisat2-build data/genome.fa genome_hisat2
</code></pre>
<p><strong>1.2.3 Extract splicing sites (if reference annotation is available) using HISAT2</strong></p>
<pre><code>gffread data/maker/ref.gff3 -T -o ref.gtf
hisat2_extract_splice_sites.py ref.gtf &gt; splicesites.txt
</code></pre>
<p><strong>1.2.4 Map RNA-seq read on genomic sequences</strong></p>
<pre><code>time hisat2 -p 20 -x genome_hisat2 --known-splicesite-infile splicesites.txt --dta --dta-cufflinks -1 data/RNA-seq/root_R1.fq.gz,data/RNA-seq/bud_R1.fq.gz -2 data/RNA-seq/root_R2.fq.gz,data/RNA-seq/bud_R2.fq.gz | samtools view -bS - &gt; all_runs.bam
</code></pre>
<p><code>-1</code> and <code>-2</code> are input parameters for pair-end libraries<br>
<code>-U</code> is for single end libraries</p>
<p>Mapping result is <code>all_runs.bam</code> in this example.</p>
<p><strong>1.2.5 Sort BAM file using sambamba</strong></p>
<pre><code>sambamba sort -m 40G --tmpdir tmp/ -o all_runs.sorted.bam -p -t 20 all_runs.bam
</code></pre>
<p>Sorted BAM result is <code>all_runs.sorted.bam</code></p>
<p><strong>1.2.6 Reference-based transcriptome assemble</strong></p>
<pre><code>stringtie all_runs.sorted.bam -o transcriptome_models.gtf -p 20
cufflinks2gff3  transcriptome/transcriptome_models.gtf &gt; transcriptome/transcriptome_models.gff3
</code></pre>
<p>In above commands, <code>-p 20</code> or <code>-t 20</code> is the number of CPU cores assigned to the program. Type <code>nproc</code> to check the maximum number in your computer. <code>-m 40G</code> is max RAM size assigned to your computer. Type <code>free</code> to check your computer RAM size.<br>
<code>transcriptome/transcriptome_models.gff3</code> contains transcriptome data which is expression evidence in Maker genome annotation (next step)</p>
<h3 id="optimized-genome-annotation-procedure-for-mining-ssp-gene-using-maker-pipeline">1.3 Optimized genome annotation procedure for mining SSP gene using Maker pipeline</h3>
<p>General genome annotation procedure can be optimized to identify more SSP gene through including SSP-specific expression evidence and conserved known SSP domain.</p>
<p><strong>1.3.1. prepare Maker configuration file</strong><br>
The protocol for genome annotation using Maker has been well documented [ref]. We installed and tested Maker pipeline in the Docker image. Users need to generate three maker_opts files: <code>maker_opts_1.ctl</code>, <code>maker_opts_2.ctl</code> and <code>maker_opts_3.ctl</code>. In addition, Maker also need maker_bopts.ctl and maker_exe.ctl files. These files include input data files path and other settings for genome annotation. Maker will take these files  as input to generate final gff file with genome annotation information. The annotation procedure will be typically invoked for three rounds to generate optimized result.</p>
<p>The GFF file for transcriptome generated by previous step and known SSP protein sequences (as of 01/2019, under /work/ssp/data) will be put in above three maker_opts files (line 18 and 23). The extra information will help maker to identify novel SSP genes.</p>
<p><strong>1.3.2 Run MAKER pipeline</strong></p>
<p>We included all Maker configuration files in demo data. Users should be able to test Maker for the first round using the following command:</p>
<pre><code>cd /work/ssp
time /usr/lib64/mpich/bin/mpiexec -n 20 maker -fix_nucleotides maker_opts_1.ctl maker_bopts.ctl maker_exe.ctl  1&gt;&amp;2 2&gt;log
</code></pre>
<p><code>-n 20</code> is the number of CPU cores. Please check <code>log</code> file if pipeline failed.</p>
<p>To generate GFF file, run this command:</p>
<pre><code>gff3_merge -d genome.maker.output/genome_master_datastore_index.log -s &gt; maker_all.gff
(head -1 maker_all.gff;cat maker_heavy.gff|grep -P $'\t(CDS|contig|exon|five_prime_UTR|gene|mRNA|three_prime_UTR)\t') &gt; maker.gff
</code></pre>
<h3 id="ssp-genome-annotation-using-spada-pipeline">1.4 SSP genome annotation using SPADA pipeline</h3>
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
<h3 id="merge-the-annotation-result-from-maker-and-spada">1.5 Merge the annotation result from Maker and SPADA</h3>
<p>There are (partially) duplicate genes between Maker and SPADA output. Additionally, users also need to consider the duplicate gene model in reference annotation if you are working on a model plant. Here we describe the policy of  identifying duplicated gene list:</p>
<p><strong>1.5.1 Generate cds sequence from different annotation using gff and genomic sequences</strong></p>
<p>Generate CDS sequences for SPADA gene models:</p>
<pre><code>gffread sspanno/31_model_evaluation/61_final.gff -g data/genome.fa -x spada_cds.fa
</code></pre>
<p>Replace transcript id with gene id:</p>
<pre><code>  sed -ri 's/^&gt;\S+\s+gene=/&gt;/' spada_cds.fa     
</code></pre>
<p><strong>1.5.2 Run NCBI BLASTN between two generated cds files and only keep the query-hit gene pairs with &gt; 50% overlapped region.</strong></p>
<p><strong>1.5.3 Check the coordinates of the query-hit pairs in gff file</strong></p>
<p>If both genes share overlapped region and the same direction on chromosome, choose one of them as redundant gene**</p>
<p><strong>1.5.4 Remove duplicate genes from corresponding GFF file</strong></p>
<pre><code>gffremove.py --infile sspanno/31_model_evaluation/61_final.gff --outfile new_spada.gff --genefile nr_gene_in_spada.txt
</code></pre>
<p><code>nr_gene_in_spada.txt</code> includes gen list to be removed. <code>new_spada.gff</code> is new gff file.</p>
<p><strong>1.5.5 Merge two GFF files</strong></p>
<p>Take <code>maker.gff</code> and <code>new_spada.gff</code> as example</p>
<pre><code>    head -1 maker.gff  &gt; all.gff
    grep -v '^#' maker.gff new_spada.gff  &gt;&gt; maker_spada.gff
</code></pre>
<p><strong>1.5.6 Generate protein and transcript files</strong></p>
<p>Protein sequences:</p>
<pre><code>gffread all.gff -g data/genome.fa -y all_protein.fa
sed -ri 's/^&gt;\S+\s+gene=/&gt;/' all_protein.fa
</code></pre>
<p>Transcript sequences:</p>
<pre><code>gffread all.gff -g data/genome.fa -w all_transcript.fa
</code></pre>
<p>The protein and transcript files will be used to further annotate gene function in protocol #2.</p>
<h2 id="protocol-2-functional-annotation-and-family-classification-of-ssp-genes">Protocol #2: Functional annotation and family classification of SSP genes</h2>
<p>Due to the short conserved region and less homologous among the member of the same gene family, it is less efficient to search and identify SSP protein using standard NCBI BLASTP.  Here we introduce a comprehensive  annotation procedure to identify SSP gene from candidate genes.</p>
<h3 id="prerequisite-1">2.1 Prerequisite</h3>
<p>The Prerequisite is almost the same as protocol #1 (see 1.1) with only exception on input data.</p>
<p><strong>2.1.1 Input data</strong></p>
<ul>
<li>SSP gene expression specific RNA-seq data in compress fastq format</li>
<li>Protein sequence (with gene ID as protein ID) of SSP gene candidates in FASTA format</li>
<li>Transcript sequences of SSP gene  candidates and a two-columns mapping file between gene and transcript id.</li>
<li>Known SSP protein sequences and HMM library; Both data are available in Docker image and demo data.</li>
</ul>
<h3 id="only-keep-short-sequences-250-a.a.">2.2 Only keep short sequences (&lt;250 a.a.)</h3>
<pre><code>keepshortseq all_protein.fa 250 &gt; short-seq.fa
</code></pre>
<h3 id="smith-waterman-search-against-known-ssp-protein">2.3 Smith-waterman search against known SSP protein</h3>
<p>Smith-Waterman alignment is more accurate way for sequence homolog search compared to BLAST search. The wrapped shell script <code>swsearch</code> use <a href="https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml">FASTA</a> software to perform Smith-waterman search against known SSP proteins (<code>target.fa</code> in demo data) and take top 2 hits (Expection&lt;0.01) as output. The script will remove hit sequence ID and only output family name.</p>
<pre><code>swsearch short-seq.fa /work/ssp/data/ssp_family.fa 0.01 &gt; sw.txt
</code></pre>
<p><strong>E-value option to be added</strong><br>
The final result file is <code>sw.txt</code> in this case.</p>
<h3 id="hmm-search-against-hmms-of-known-ssp-families">2.4 HMM search against HMMs of known SSP families</h3>
<p><strong>2.4.1  Generate HMM library for all known SSP families from SPADA installation</strong></p>
<pre><code>cat /opt/spada_soft/spada/CRP_PlantSSPv1_Noble/15_hmm/*.hmm &gt; all.hmm
</code></pre>
<p><strong>2.4.2  compile HMM library</strong></p>
<pre><code>hmmpress  all.hmm 
</code></pre>
<p><strong>2.4.3  search your protein sequence against HMM library</strong></p>
<pre><code>hmmscan --cpu 20 -E 0.01 --tblout hmm_output.txt all.hmm short-seq.fa &gt; /dev/null
</code></pre>
<p>The expectation cutoff <code>-E</code> for <code>hmmscan</code> is <code>0.01</code>. <code>short-seq.fa</code> is input protein file and <code>all.hmm</code> is HMM library files.</p>
<p>The above commands will generate a tab-delimited table file <code>hmm_output.txt</code>.  In the result table, column #1 is HMM family name and column #2 is gene ID in input protein sequence file.</p>
<h3 id="signal-peptide-detection">2.5 Signal peptide detection</h3>
<p>We choose SignalP to predict signal peptide from SSP candidate peptides. We recommend to use “No TM” and long output format. D-score thresholds are usually 0.45 or 0.5, depending on the type of network chosen (with or without transmembrane segments), but we recommend to use a D-score of ≥ 0.45 for known SSPs or ≥ 0.25 for putative SSPs.</p>
<pre><code>/opt/spada_soft/signalp-4.1/signalp -t euk -f long -s notm short-seq.fa &gt; signalp_long.txt
cat signalp_long.txt | singalP_parser  &gt; sp.txt
</code></pre>
<p>In this case, <code>/opt/spada_soft/signalp-4.1/signalp</code> is Signalp program which output prediction result file <code>signalp_long.txt</code>.  <code>singalP_parser</code> is a script to parse long format output of SignalP. The final signal peptide prediction result <code>sp.txt</code> includes four columns:  gene ID,  start coordinates,  end coordinates,  D-value, cut-off and conclusion(YES/NO).</p>
<h3 id="identification-of-novel-ssp-gene-families-using-mcl-analysis">2.6 Identification of novel SSP gene families using MCL analysis</h3>
<p>SSP candidate (signalp D-value&gt;0.25)   can be clustered into candidate SSP families using Markov Chain Cluster (MCL)[PMID:11917018]. The procedure should be performed on last 50 a.a. and candidate  pipetide should be less than 230 a.a…</p>
<p><strong>2.6.1. Create index to retrieve protein sequence by sequence ID.</strong></p>
<pre><code>cdbfasta short-seq.fa
</code></pre>
<p><strong>2.6.2 Only will take the protein with D-value&gt;0.45.</strong></p>
<pre><code>cat sp.txt | awk '{if($4&gt;0.45) print $1}' | cdbyank short-seq.fa.cidx  &gt; all_putative_ssp.fa
</code></pre>
<p><strong>2.6.3 The candidate proteins should be shorter than 230 a.a. and only take the last 50 a.a…</strong></p>
<pre><code>shortseqtail all_putative_ssp.fa 230 50 &gt; peptide-tail.fa
</code></pre>
<p><strong>2.6.4 Generate protein vs protein relationship file.</strong></p>
<pre><code>search35_t -T 20 -Q -H -m 9 -b 100 -d 100 peptide-tail.fa peptide-tail.fa &gt; sw_for_peptidetail

bioparser -t ssearch -m sw_for_peptidetail | awk '{print $3,$6,$14}' FS="\t" | sort | uniq | awk '{ if($1!=$2&amp;&amp;$3 &lt; 0.01)) print $0; }' FS=" " &gt; protein-protein-rel.txt
</code></pre>
<p><code>protein-protein-rel.txt</code> is a three columns file with two protein/gene IDs and their relation in e-value. All protein-protein pair with e-vakue&gt;0.01 have been removed.</p>
<p><strong>2.6.5 Cluster  protein into clusters</strong></p>
<pre><code>mcxload -abc protein-protein-rel.txt --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o protein-protein.mci -write-tab protein-protein.tab
mcl last50seq.mci -I 1.4 -use-tab protein-protein.tab
</code></pre>
<p>In this case, <code>mcl</code> command will generate <code>out.last50seq.mci.I14</code>, in which all gene/protein belonging to the same cluster will be put in one line.</p>
<p>Referring to <a href="https://micans.org/mcl/">MCL manual</a> to adjust <code>-I</code>. The current value for <code>-I</code> generate larger clusters. You may consider to increase the value to get smaller clusters.</p>
<p>Next, type following command to generate a cluster name vs gene/protein ID table file.</p>
<pre><code>cat out.last50seq.mci.I14 | awk '{print "Cluster_" NR "\t" $0}' |awk '{ OFS="\n" $1 "\t";$1="";print $0;}'|grep -Pv '^\s*$' &gt; mclcluster_protein.txt
</code></pre>
<p><code>mclcluster_protein.txt</code> is the result file containing cluster-protein mapping.</p>
<h3 id="gene-expression-analysis-of-rna-seq-data">2.7 Gene expression analysis of RNA-seq data</h3>
<p>Use the published software RSEM [PMID:21816040]  to generate gene expression table using RNA-seq  data. In generated table, each gene will occupy one row and each RNA-seq sample will take one column. We recommend transcript per million transcripts (TPM) value in this table.</p>
<h3 id="perform-transmembrane-helix-tmh-prediction">2.8 Perform transmembrane helix (TMH) prediction</h3>
<p>TMH prediction is a criterion used to classify the putative SSPs (de Bang et al., 2017) because a gene harboring TMHs cannot be considered as an SSP.</p>
<p><strong>2.8.1  Remove the N-terminal signal peptide from the input protein sequence(s)</strong></p>
<p>Signal peptides tend to display high hydrophobicity and are often reported as false positive prediction. [to-be-complete]</p>
<p><strong>2.8.2  Perform TMH prediction</strong><br>
We recommend the TMHMM server because it is easy-to-use and accurate due to its based on hidden Markov model search (Krogh et al., 2001). TMHMM is user-friendly web-based analysis tool. Just submit your protein sequence(s) input file, and select the output format with graphics. The output result provides a list of the predicted TMHs and their locations along with statistical information such as protein length. In our plant SSP pipeline , TMH prediction is a criterion used to classify the putative SSPs.</p>
<h3 id="comprehensive-table-of-gene-annotation-evidence">2.9 Comprehensive table of gene annotation evidence</h3>
<p>The result file generated from Section 2.3 to 2.8 should be merged  into a comprehensive data table in which each gene will use one row and each annotation information, such as RNA-seq sample, family name from Smith-Waterman search and HMM search, D-value from signalP and cluster ID from MCL analysis will be take one column. Microsoft Excel is a nice tool to merge and generate such comprehensive table.</p>
<p>Further curation based on the table will help to screen known and putative SSP genes. The genes with Smith-Waterman or HMM hits will be considered as known SSP genes. Other genes with D-value(&gt;0.25) and shorter than 230 a.a. will be  considered as putative SSP genes. The gene expression value is also helpful to identify gene with high confidence.</p>
<h2 id="trouble-shooting">Trouble shooting</h2>
<p><strong>You cannot attached to a stopped container, start it first</strong></p>

