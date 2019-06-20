---


---

<h2 id="protocol-1-small-secret-peptide-gene-discovery-from-genomic-sequences">Protocol #1: Small Secret Peptide Gene discovery from genomic sequences</h2>
<p>Traditional genome annotation policy is biased to discover long genes; leading to missing of some small secret peptide (SSP) genes. The following workflow was optimized to identify SSP genes from assembled genomic sequences utilizing specific RNA-seq data as expression evidence and conserved SSP motifs.</p>
<h3 id="prerequisites">1.1. Prerequisites</h3>
<p><strong>1.1.1. Suggestions:</strong></p>
<ol>
<li>We recommend an X-windows desktop (such as gnome/XFCE/MATE) instead of SSH terminal because it is more convenient to edit files.</li>
<li>All commands below are typed under Linux terminal.</li>
<li>A line start with <code>#</code> in Linux command line indicates that this is explanatory information only.</li>
<li>You need a no-root user with sudo privilege in host system to install docker packages and enable docker service. Asking your system administrator to install docker service and add your user as a member of docker group if you can’t have <code>sudo</code> privileges.</li>
<li>You need to be a sudo user or a member of <code>docker</code> group in host system to start Docker container and attach to container terminal.</li>
<li>Default user in the Docker container is <code>test</code>.</li>
</ol>
<p><strong>1.1.2. Computer:</strong><br>
A  high-performance computer (I7/Xeon processor and &gt;16GB RAM) with CentOS 7, Ubuntu 16.04 or higher as your host operation system(OS).</p>
<p><strong>1.1.3. Work folder</strong><br>
Work folder is the place for all raw input data (genomic sequences, gff, RNA-seq and protein, etc ), and analysis results in your host OS. It is recommend to be <code>work</code> under your home directory.  For example, if your username is <code>test</code> in host OS, the recommended work folder will be <code>/home/test/work</code> in your host OS. To create the work folder in your home directory of host OS:</p>
<pre><code>cd ~   # ~ means your home directory, e.g. /home/test
mkdir work
</code></pre>
<p><strong>1.1.4. Input data</strong></p>
<ul>
<li>Genomics sequences in FASTA format</li>
<li>Reference annotation in GFF format if available</li>
<li>SSP gene expression specific RNA-seq data in compress FASTQ format</li>
<li>Protein sequence of known SSP genes other related protein sequences in FASTA format</li>
<li>Other EST/transcript sequences from the same species.</li>
</ul>
<p><strong>1.1.5. Demo data</strong><br>
The demo data is available for <a href="http://bioinfo.noble.org/manuscript-support/ssp-protocol/ssp-demo.tar.gz">download</a>. <strong>In host OS</strong>, copy it to your work folder and type the following command to unzip it:</p>
<pre><code>cd ~/work
wget http://bioinfo.noble.org/manuscript-support/ssp-protocol/ssp-demo.tar.gz
tar -xzvf ssp-demo.tar.gz
</code></pre>
<p>The above command will generate <code>ssp</code> folder under <code>work</code>, download the demo file <code>ssp-demo.tar.gz</code>, and uncompress it.</p>
<p>In the <code>~/work/ssp/data</code> folder, <code>ssp_family.fa</code> is a protein sequences of known SSP genes. The known SSP file is used in Maker genome annotation (Protocol #1) and SSP gene annotation (Protocol #2).</p>
<p><strong>1.1.6. Software installation</strong><br>
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
<p>The above commands will start a Docker container named <code>sspvm</code> using  <code>docker.io/noblebioinfo/sspgene</code> as template image. This step will take a while depend on your network download speed.</p>
<p>In  <code>-v $(pwd)/work:/work</code>:  <code>$(pwd)/work</code>, the path of work folder in your host OS, is <code>work</code> under your current directory. Here,<code>$(pwd)</code> will be converted to your current folder, e.g. home folder by Linux Bash interpreter. The work folder in host OS will be mounted on <code>/work</code> in Docker container. Thus, the folder makes it possible to exchange data between Host computer (<code>$(pwd)/work</code>) and Docker container  (<code>/work</code>). You can copy your demo data or other research data to the work folder in hosts OS (<code>$(pwd)/work</code>) and access them in <code>/work</code> in Docker container.</p>
<p>The <code>attach</code> subcommand will link your current Linux terminal to the running docker container (<code>bioinfo</code> in this case).<br>
Tip:  to detach the container terminal and get back to host OS <code>hold Ctrl key and press P,Q</code>.</p>
<p>Type the following command to enter demo data folder work folder in attached Docker container terminal:</p>
<pre><code>cd /work/ssp
</code></pre>
<p>All Linux commands below should be typed in this container terminal.</p>
<h3 id="prepare-rna-sed-based-gene-expression-evidence-for-maker-pipeline">1.2. Prepare RNA-sed based gene expression evidence for MAKER pipeline</h3>
<p>Some plant SSP genes may only express under a specific condition or tissue, such as nutrient deficiency or root tissue. Related RNA-seq data will help to improve the performance of SSP gene mining. The following sample code will perform  reference-based transcriptome assembly and generate a GFF file for MAKER genome annotation.</p>
<p><strong>1.2.1. Prepare work folder</strong></p>
<pre><code>cd /work/ssp
mkdir transcriptome
cd transcriptome/
</code></pre>
<p><strong>1.2.2. Compile the genomics sequences using HISAT2</strong></p>
<pre><code>hisat2-build /work/ssp/data/genome.fa genome_hisat2
</code></pre>
<p><strong>1.2.3. Extract splicing sites (if reference annotation is available) using HISAT2</strong></p>
<pre><code>gffread /work/ssp/data/maker/ref.gff3 -T -o ref.gtf
hisat2_extract_splice_sites.py ref.gtf &gt; splicesites.txt
</code></pre>
<p><strong>1.2.4. Map RNA-seq read on genomic sequences</strong></p>
<pre><code>time hisat2 -p 20 -x genome_hisat2 --known-splicesite-infile splicesites.txt --dta --dta-cufflinks -1 /work/ssp/data/RNA-seq/root_R1.fq.gz,/work/ssp/data/RNA-seq/bud_R1.fq.gz -2 /work/ssp/data/RNA-seq/root_R2.fq.gz,/work/ssp/data/RNA-seq/bud_R2.fq.gz | samtools view -bS - &gt; all_runs.bam
</code></pre>
<p><code>-1</code> and <code>-2</code> are input parameters for paired-end libraries, and <code>-U</code> is the input parameter for single-end libraries.<br>
<code>all_runs.bam</code> file is the mapping result file.</p>
<p><strong>1.2.5. Sort BAM file using sambamba</strong></p>
<pre><code>sambamba sort -m 40G --tmpdir tmp/ -o all_runs.sorted.bam -p -t 20 all_runs.bam
</code></pre>
<p><code>all_runs.sorted.bam</code> is the sorted BAM file.</p>
<p><strong>1.2.6. Generate reference-based transcriptome file</strong></p>
<pre><code>stringtie all_runs.sorted.bam -o transcriptome_models.gtf -p 20
cufflinks2gff3  /work/transcriptome/transcriptome_models.gtf &gt; /work/transcriptome/transcriptome_models.gff3
</code></pre>
<p>In the above commands, <code>-p 20</code> or <code>-t 20</code> is the number of CPU cores assigned to the program. Type <code>nproc</code> to check the maximum number in your computer. <code>-m 40G</code> is max RAM size assigned to your computer. Type <code>free</code> to check your computer RAM size.</p>
<p><code>transcriptome_models.gff3</code> is the output file and contains transcriptome data. This file will be used as expression evidence in MAKER genome annotation (step 1.3.2.).</p>
<h3 id="genome-annotation-procedure-for-mining-ssp-genes-using-maker-pipeline">1.3. Genome annotation procedure for mining SSP genes using MAKER pipeline</h3>
<p>General genome annotation procedure can be optimized to identify more SSP genes through including SSP-specific expression evidence and conserved known SSP domains.</p>
<p><strong>1.3.1. Prepare MAKER configuration file</strong></p>
<p>The protocol for genome annotation using MAKER has been well documented (<a href="http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Main_Page">http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Main_Page</a>). We installed and tested MAKER pipeline in the Docker image.</p>
<p>Users need to generate three MAKER configuration files called <code>maker_opts</code>: <code>maker_opts_1.ctl</code>, <code>maker_opts_2.ctl</code> and <code>maker_opts_3.ctl</code>. In addition, MAKER also needs <code>maker_bopts.ctl</code> and <code>maker_exe.ctl</code> configuration files. These files include paths for the input data files and other settings for the genome annotation. MAKER will take these files as inputs to generate the final GFF file with genome annotation information. The annotation procedure will be done for three rounds to generate optimized results.</p>
<p>The GFF file for transcriptome generated in the previous step (1.2.6.) and known SSP protein sequences (as of 01/2019, under /work/ssp/data) will be included in above three <code>maker_opts</code> files. This additional information will help MAKER to identify novel SSP genes.</p>
<p><strong>1.3.2. Run MAKER pipeline</strong></p>
<p>We generated the optimized gene models using SNAP gene predictior for <em>Medicago truncatula</em>. If you want to use these optimized gene models, skip to Round 3. But if you are going to generate these files (for another species, for example), run rounds 1, 2 and 3.</p>
<p><strong>1.3.2.1. Round 1</strong></p>
<p>We included all MAKER configuration files in the demo data. Users should be able to test MAKER for the first round using the following command:</p>
<pre><code>cd /work/ssp
time /usr/lib64/mpich/bin/mpiexec -n 20 maker -fix_nucleotides /work/ssp/data/maker/maker_opts_1.ctl /work/ssp/data/maker/maker_bopts.ctl /work/ssp/data/maker/maker_exe.ctl 1&gt;&amp;2 2&gt;log_round1
</code></pre>
<p><code>-n 20</code> is the number of CPU cores. Please check the <code>log</code> file if pipeline failed.</p>
<p>Then, generate the GFF file using the following command:</p>
<pre><code>gff3_merge -d genome.maker.output/genome_master_datastore_index.log -s &gt; maker_all_round1.gff
(head -1 maker_all_round1.gff;cat maker_all_round1.gff|grep -P $'\t(CDS|contig|exon|five_prime_UTR|gene|mRNA|three_prime_UTR)\t') &gt; maker_round1.gff
</code></pre>
<p>Then, generate optimized gene models with SNAP gene predictior using the following commands:</p>
<pre><code>mkdir snap
cd snap
maker2zff  -d /work/ssp/genome.maker.output/genome_master_datastore_index.log
/opt/maker/exe/snap/fathom -categorize 1000 genome.ann genome.dna
/opt/maker/exe/snap/forge export.ann export.dna
/opt/maker/exe/snap/hmm-assembler.pl mt . &gt; mt.hmm
cd ..
</code></pre>
<p><strong>1.3.2.2. Round 2</strong></p>
<p>In the second round, users need to include an HMM file generated by the SNAP gene predictor on Round 1. This HMM file (<code>mt.hmm</code>) is included in the demo data, and also the <code>maker_opts_2.ctl</code> file.</p>
<p>Run MAKER again:</p>
<pre><code>time /usr/lib64/mpich/bin/mpiexec -n 20 maker -fix_nucleotides /work/ssp/data/maker/maker_opts_2.ctl /work/ssp/data/maker/maker_bopts.ctl /work/ssp/data/maker/maker_exe.ctl  1&gt;&amp;2 2&gt;log_round2
</code></pre>
<p>Then, generate the GFF file again using the following command:</p>
<pre><code>gff3_merge -d genome.maker.output/genome_master_datastore_index.log -s &gt; maker_all_round2.gff
(head -1 maker_all_round2.gff;cat maker_all_round2.gff|grep -P $'\t(CDS|contig|exon|five_prime_UTR|gene|mRNA|three_prime_UTR)\t') &gt; maker_round2.gff
</code></pre>
<p>Then, generate optimized gene models again with SNAP gene predictior using the following commands:</p>
<pre><code>mkdir snap2
cd snap2
maker2zff  -d /work/ssp/genome.maker.output/genome_master_datastore_index.log
/opt/maker/exe/snap/fathom -categorize 1000 genome.ann genome.dna
/opt/maker/exe/snap/forge export.ann export.dna
/opt/maker/exe/snap/hmm-assembler.pl mt2 . &gt; mt2.hmm
cd ..
</code></pre>
<p><strong>1.3.2.3. Round 3</strong></p>
<p>In the third round, users need to include an HMM file generated by the SNAP gene predictor in the Round 2. This HMM file (<code>mt2.hmm</code>) is included in the demo data, and also the <code>maker_opts_3.ctl</code> file.</p>
<p>Run MAKER again:</p>
<pre><code>time /usr/lib64/mpich/bin/mpiexec -n 20 maker -fix_nucleotides /work/ssp/data/maker/maker_opts_3.ctl /work/ssp/data/maker/maker_bopts.ctl /work/ssp/data/maker/maker_exe.ctl  1&gt;&amp;2 2&gt;log_round3
</code></pre>
<p>Then, generate the GFF file again using the following command:</p>
<pre><code>gff3_merge -d genome.maker.output/genome_master_datastore_index.log -s &gt; maker_all_round3.gff
(head -1 maker_all_round3.gff;cat maker_all_round3.gff|grep -P $'\t(CDS|contig|exon|five_prime_UTR|gene|mRNA|three_prime_UTR)\t') &gt; maker_round3.gff
</code></pre>
<h3 id="genome-annotation-procedure-for-mining-ssp-genes-using-spada-pipeline">1.4. Genome annotation procedure for mining SSP genes using SPADA pipeline</h3>
<p>SPADA pipeline typically utilizes conserved domains (in HMM format) information of known SSP families to identify SSP genes from genomic sequences (Zhou et al., 2013). We included a copy of SPADA and a comprehensive HMM dataset from PlantSSP database (Ghorbani et al., 2015), and our curated known SSPs in Docker image.</p>
<p>Here is an example of SPADA analysis:</p>
<pre><code>perl /opt/spada_soft/spada/spada.pl --cfg /opt/spada_soft/spada/cfg.txt -d sspanno -p /opt/spada_soft/spada/CRP_PlantSSPv1_Noble -f data/genome.fa -t 20 -o arabidopsis 1&gt;&amp;2 2&gt;spada_log
</code></pre>
<p>In this example:<br>
<code>/opt/spada_soft/spada/</code>: SPADA pipeline is installed here.<br>
<code>/opt/spada_soft/spada/CRP_PlantSSPv1_Noble</code>: location of above-mentioned comprehensive HMM dataset. You can change it to any of your favorite HMM dataset.<br>
<code>data/genome.fa</code>: is the genomic sequences to be analyzed.<br>
<code>/opt/spada_soft/augustus/config/species/</code>: change it to the closest species in your case.<br>
<code>arabidopsis</code>: folder name under <code>/opt/spada_soft/augustus/config/species/</code><br>
The GFF  file is available at <code>sspanno/31_model_evaluation/61_final.gff</code></p>
<h3 id="merge-the-annotation-result-from-maker-and-spada">1.5. Merge the annotation result from MAKER and SPADA</h3>
<p>There are (partially) duplicate genes between MAKER and SPADA outputs. Additionally, users also need to consider which gene model to be removed. If you are working with a model plant, it is better to retain the official annotation.<br>
Here we describe the procedure to identify duplicated genes:</p>
<p><strong>1.5.1. Generate cds sequences from different annotations using GFF and genomic sequences</strong></p>
<p>Generate CDS sequences for SPADA gene models:</p>
<pre><code>gffread sspanno/31_model_evaluation/61_final.gff -g data/genome.fa -x spada_cds.fa
</code></pre>
<p>Replace transcript id with gene id:</p>
<pre><code>sed -ri 's/^&gt;\S+\s+gene=/&gt;/' spada_cds.fa     
</code></pre>
<p><strong>1.5.2. Run NCBI BLASTN between two generated CDS files and only keep the query-hit gene pairs with &gt; 50% overlapped region.</strong></p>
<p><strong>1.5.3. Check the coordinates of the query-hit pairs in the MAKER and SPADA GFF files</strong></p>
<p>If both genes have an overlapped region and same coordinates, choose one of them as redundant gene**</p>
<p>First create bed files from the GFF files for MAKER and SPADA:</p>
<pre><code>cat sspanno/31_model_evaluation/61_final.gff | awk 'BEGIN{FS="\t"}; BEGIN{OFS="\t"}; $3=="gene"' | tr '; ' '\t' | awk 'BEGIN{FS="\t"}; BEGIN{OFS="\t"};{print $1,$4,$5,$9}' | sed 's/ID=//' &gt; spada.bed
cat maker_round3.gff | awk 'BEGIN{FS="\t"}; BEGIN{OFS="\t"}; $3=="gene"' | tr '; ' '\t' | awk 'BEGIN{FS="\t"}; BEGIN{OFS="\t"};{print $1,$4,$5,$9}' | sed 's/ID=//' &gt; maker.bed
</code></pre>
<p>Then, check overlaps between MAKER and SPADA, and create duplicate regions from SPADA:</p>
<pre><code>bedtools intersect -wa -wb -a maker.bed -b spada.bed -filenames | cut -f8 &gt; dup_gene_in_spada.txt
</code></pre>
<p><strong>1.5.4. Remove duplicate genes from the corresponding GFF file</strong></p>
<pre><code>gffremove.py --infile sspanno/31_model_evaluation/61_final.gff --outfile new_spada.gff --genefile dup_gene_in_spada.txt
</code></pre>
<p><code>dup_gene_in_spada.txt</code> includes the gene list to be removed. <code>new_spada.gff</code> is the new GFF file.</p>
<p><strong>1.5.5. Merge MAKER and SPADA GFF files without redundant genes</strong></p>
<pre><code>head -1 maker_round3.gff  &gt; header
cat maker_round3.gff new_spada.gff | grep -v "#" &gt; maker_spada.gff
cat header maker_spada.gff &gt; all.gff
</code></pre>
<p><strong>1.5.6. Generate protein and transcript files</strong></p>
<p>The protein and transcript files will be used to further annotate gene functions in the Protocol #2.</p>
<p>Protein sequences:</p>
<pre><code>gffread all.gff -g data/genome.fa -y all_protein.fa
sed -ri 's/^&gt;\S+\s+gene=/&gt;/' all_protein.fa
</code></pre>
<p>Transcript sequences:</p>
<pre><code>gffread all.gff -g data/genome.fa -w all_transcript.fa
</code></pre>
<h2 id="protocol-2-functional-annotation-and-family-classification-of-ssp-genes">Protocol #2: Functional annotation and family classification of SSP genes</h2>
<p>Due to the short conserved regions and less homologous among the members of the same gene family, it is less efficient to search and identify SSP proteins using standard NCBI BLASTP.  Here we introduce a comprehensive annotation procedure to identify SSPs from candidate genes.</p>
<h3 id="prerequisite">2.1. Prerequisite</h3>
<p>The Prerequisite is almost the same as protocol #1 (see 1.1), except for the input data.</p>
<p><strong>2.1.1. Input data</strong></p>
<ul>
<li>SSP gene expression evidence (RNA-seq data) in compressed FASTQ format files</li>
<li>Protein sequence (with gene ID as protein ID) of SSP gene candidates in FASTA format</li>
<li>Transcript sequences of SSP gene candidates and a two-columns mapping file between gene and transcript id.</li>
<li>Known SSP protein sequences and HMM library; both data are available in the Docker image and demo data.</li>
</ul>
<h3 id="only-keep-short-sequences--250-a.a.">2.2. Only keep short sequences (&lt; 250 a.a.)</h3>
<pre><code>keepshortseq all_protein.fa 250 &gt; short-seq.fa
</code></pre>
<h3 id="smith-waterman-search-against-known-ssp-protein">2.3. Smith-waterman search against known SSP protein</h3>
<p>Smith-Waterman alignment is more accurate way for sequence homolog search compared to BLAST search. The wrapped shell script <code>swsearch</code> use <a href="https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml">FASTA</a> software to perform Smith-waterman search against known SSP proteins (<code>target.fa</code> in demo data) and take two top hits (e-value &lt; 0.01) as output. The script will remove hit sequence IDs and will only output the family names.</p>
<pre><code>swsearch short-seq.fa /work/ssp/data/ssp_family.fa 0.01 &gt; sw.txt
</code></pre>
<h3 id="hmm-search-against-hmms-of-known-ssp-families">2.4. HMM search against HMMs of known SSP families</h3>
<p><strong>2.4.1.  Generate HMM library for all known SSP families from SPADA installation</strong></p>
<pre><code>cat /opt/spada_soft/spada/CRP_PlantSSPv1_Noble/15_hmm/*.hmm &gt; all.hmm
</code></pre>
<p><strong>2.4.2 . Compile HMM library</strong></p>
<pre><code>/opt/spada_soft/hmmer/bin/hmmpress all.hmm
</code></pre>
<p><strong>2.4.3  Search your protein sequence against HMM library</strong></p>
<pre><code>/opt/spada_soft/hmmer/bin/hmmscan --cpu 20 -E 0.01 --tblout hmm_output.txt all.hmm short-seq.fa &gt; /dev/null
</code></pre>
<p>The expectation cutoff <code>-E</code> for <code>hmmscan</code> is <code>0.01</code>. <code>short-seq.fa</code> is the input protein file and <code>all.hmm</code> is the HMM library file.</p>
<p>The above command will generate a tab-delimited table file <code>hmm_output.txt</code>.  In the result table, column #1 is HMM family name and column #2 is gene ID in the input protein sequence file.</p>
<h3 id="signal-peptide-detection">2.5. Signal peptide detection</h3>
<p>We use <a href="http://www.cbs.dtu.dk/services/SignalP-4.0/">SignalP</a> to predict the signal peptides from SSP candidates. We recommend to use “No TM” and long output format. <em>D-score</em> thresholds are usually 0.45 or 0.5, depending on the type of network chosen (with or without transmembrane segments), but we recommend a <em>D-score</em> of ≥ 0.45 for known SSPs or ≥ 0.25 for putative SSPs.</p>
<pre><code>/opt/spada_soft/signalp-4.1/signalp -t euk -f long -s notm short-seq.fa &gt; signalp_long.txt
cat signalp_long.txt | singalP_parser  &gt; sp.txt
</code></pre>
<p><code>signalp_long.txt</code>: is the output prediction result file<br>
<code>singalP_parser</code>: is a script to parse long format outputs from SignalP.<br>
The final signal peptide prediction result <code>sp.txt</code> includes four columns:  gene ID,  start coordinate,  end coordinate,  <em>D-score</em>, cut-off, and SSP prediction (YES/NO).</p>
<h3 id="identification-of-novel-ssp-gene-families-using-mcl-analysis">2.6. Identification of novel SSP gene families using MCL analysis</h3>
<p>SSP candidates (SignalP <em>D-score</em> &gt; 0.25)  can be clustered into candidate SSP families using Markov Chain Cluster (MCL). The procedure should be performed on the last 50 a.a. and candidate peptides should be less than 230 a.a…</p>
<p><strong>2.6.1. Create index to retrieve protein sequences by sequence ID.</strong></p>
<pre><code>cdbfasta short-seq.fa
</code></pre>
<p><strong>2.6.2 Select proteins with <em>D-score</em> &gt; 0.45.</strong></p>
<pre><code>cat sp.txt | awk '{if($4&gt;0.45) print $1}' | cdbyank short-seq.fa.cidx  &gt; all_putative_ssp.fa
</code></pre>
<p><strong>2.6.3. Select proteins shorter than 230 a.a. and only take the last 50 a.a…</strong></p>
<pre><code>shortseqtail all_putative_ssp.fa 230 50 &gt; peptide-tail.fa
</code></pre>
<p><strong>2.6.4. Generate protein vs protein relationship file</strong></p>
<pre><code>/opt/bin/ssearch35_t -T 20 -Q -H -m 9 -b 100 -d 100 peptide-tail.fa peptide-tail.fa &gt; sw_for_peptidetail

bioparser -t ssearch -m sw_for_peptidetail | awk '{print $3,$6,$14}' FS="\t" | sort | uniq | awk '{ if($1!=$2&amp;&amp;$3 &lt; 0.01) print $0; }' FS=" " &gt; protein-protein-rel.txt
</code></pre>
<p><code>protein-protein-rel.txt</code> file has three columns with two protein/gene IDs and their relation in e-values. All protein-protein pairs with e-value &gt; 0.01 have been removed.</p>
<p><strong>2.6.5. Cluster protein into clusters</strong></p>
<pre><code>mcxload -abc protein-protein-rel.txt --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o protein-protein.mci -write-tab protein-protein.tab
mcl protein-protein.mci -I 1.4 -use-tab protein-protein.tab
</code></pre>
<p>In this case, <code>mcl</code> command will generate the output file <code>out.protein-protein.mci.I14</code>, in which all gene/protein belonging to the same cluster will be put in one line.</p>
<p>Referring to <a href="https://micans.org/mcl/">MCL manual</a> to adjust <code>-I</code>. The current value for <code>-I</code> generate larger clusters. You may consider to increase the value to get smaller clusters.</p>
<p>Next, type the following command to generate a cluster name vs gene/protein ID table file:</p>
<pre><code>cat out.protein-protein.mci.I14 | awk '{print "Cluster_" NR "\t" $0}' |awk '{ OFS="\n" $1 "\t";$1="";print $0;}'|grep -Pv '^\s*$' &gt; mclcluster_protein.txt
</code></pre>
<p><code>mclcluster_protein.txt</code> is the result file containing cluster-protein mapping.</p>
<h3 id="gene-expression-analysis-of-rna-seq-data">2.7. Gene expression analysis of RNA-seq data</h3>
<p>Use RSEM software (Li &amp; Dewey, 2011) to generate gene expression table using RNA-seq  data. In the generated table, each gene will occupy one row, and each RNA-seq sample will take one column. We recommend to use transcript per million transcripts (TPM) value in this table.</p>
<h3 id="perform-transmembrane-helix-tmh-prediction">2.8. Perform transmembrane helix (TMH) prediction</h3>
<p>TMH prediction is a criterion used to classify the putative SSPs (de Bang et al., 2017) because a gene harboring TMHs cannot be considered as an SSP.</p>
<p><strong>2.8.1.  Generate putative SSP sequences without the signal peptide regions</strong></p>
<p>Run the following command to generate the putative sequences without the signal peptide regions from the putative SSP input protein sequences:</p>
<pre><code>/opt/spada_soft/signalp-4.1/signalp -t euk -f short -m processed_putative.fa -u 0.449 -s notm all_putative_ssp.fa &gt; signalp_putative.txt
</code></pre>
<p><code>processed_putative.fa</code> is the output file that will be used in the next step.</p>
<p><strong>2.8.2.  Perform TMH prediction</strong></p>
<p>We recommend the TMHMM server because it is easy-to-use and accurate due to its based on hidden Markov model search (Krogh et al., 2001). Use the processed FASTA sequences generated in the previous step (2.8.1.).<br>
The output result file (<code>tm_putative.txt</code>) provides a list of the genes and their number of predicted transmembrane helices.</p>
<p>Run the following command:</p>
<pre><code>/opt/tmhmm-2.0c/bin/tmhmm --short processed_putative.fa &gt; tmhmm.out
cat tmhmm.out | cut -f1,5 | sed 's/PredHel=//' &gt; tm_putative.txt
</code></pre>
<h3 id="perform-er-retention-signal-search-kdel-hdel-or-kqel-in-the-c-terminal-region">2.9. Perform ER-retention signal search (KDEL, HDEL or KQEL) in the C-terminal region</h3>
<p>Only take the last 50 amino acids for this step. Run the following command to obtain the peptides with ER-retention signal peptides:</p>
<pre><code>cat peptide-tail.fa | sed -r 's/^\s+|\s+$//g'|awk '{if($0~/^&gt;/) print  "\n" $0 "\t"; else print $0;}' ORS=""|grep -Pv '^\s*$' | grep -i  'KDEL\|HDEL\|KQEL'  
</code></pre>
<p>Then, generate a FASTA file with the ER retention peptides to be excluded:</p>
<pre><code>cat peptide-tail.fa | sed -r 's/^\s+|\s+$//g'|awk '{if($0~/^&gt;/) print  "\n" $0 "\t"; else print $0;}' ORS=""|grep -Pv '^\s*$' | grep -i  'KDEL\|HDEL\|KQEL' | awk 'BEGIN{FS="\t"}; BEGIN{OFS="\t"}; {print $1,"\n"$2}' &gt; peptide-ERretention.fa
</code></pre>
<h3 id="comprehensive-table-of-gene-annotation-evidence">3.0. Comprehensive table of gene annotation evidence</h3>
<p>The result file generated from Section 2.3 to 2.8 should be merged  into a comprehensive data table in which each gene will use one row and each annotation information, such as - RNA-seq sample, family name from Smith-Waterman search and HMM search, D-score from SignalP, cluster ID from MCL analysis, and predicted transmembrane helices - will use one column. Microsoft excel is a good option to merge and generate such comprehensive table.</p>
<p>Further curation based on the table will help to screen known and putative SSP genes. The genes with Smith-Waterman or HMM hits will be considered as known SSP genes. Other genes with D-score (&gt; 0.25) and shorter than 230 a.a. will be  considered as putative SSP genes. The gene expression values is also helpful to identify gene with high confidence.</p>
<h3 id="trouble-shooting">Trouble shooting</h3>
<p>The Docker image was only tested on Linux CentOS 7 and Ubuntu 16.04 LTS. We don’t recommend to use it on Windows and Mac, although it can be installed on Windows and Mac systems. As a Linux user, installing Docker software and starting backend service require root or sudo privileges. Downloading Docker image and starting a container for the image only require to be a member of Docker user group or have root or sudo privileges. Contact your Linux administrator if you are using a virtual Linux machine without root or sudo privileges in Data center and having permission or privileges issues for running Docker container.</p>
<p>MAKER or SPADA usage errors can be found at:<br>
<a href="https://groups.google.com/forum/#!forum/maker-devel">https://groups.google.com/forum/#!forum/maker-devel</a>  <a href="https://groups.google.com/forum/#!forum/SPADA">https://groups.google.com/forum/#!forum/SPADA</a></p>
<p>HMMER help page can be found at: <a href="https://www.ebi.ac.uk/Tools/hmmer/help">https://www.ebi.ac.uk/Tools/hmmer/help</a><br>
Be advised that the HMM libraries compiled by different versions of HMMER are incompatible each other. The installed HMMER in Docker image is version 3.0. If user plans to compile HMM library using their SSP family data, please use the same version.</p>
<p>SignalP frequently asked questions can be found at:<br>
<a href="http://www.cbs.dtu.dk/services/SignalP/faq.php">http://www.cbs.dtu.dk/services/SignalP/faq.php</a></p>
<p>TMHMM instruction information can be found at <a href="http://www.cbs.dtu.dk/services/TMHMM/TMHMM2.0b.guide.php">http://www.cbs.dtu.dk/services/TMHMM/TMHMM2.0b.guide.php</a></p>
<p>MCL frequently asked questions can be found at: <a href="https://micans.org/mcl/">https://micans.org/mcl/</a></p>
<h3 id="references">References</h3>
<ul>
<li>
<p>de Bang, T. C., Lundquist, P. K., Dai, X., Boschiero, C., Zhuang,<br>
Z. , Pant, P., … Scheible, W.-R. (2017). Plant Physiology, 175<br>
(4), 1669-1689.</p>
</li>
<li>
<p>Ghorbani, S., et al. (2015). Expanding the repertoire of secretory<br>
peptides controlling root development with comparative genome<br>
analysis and functional assays. J Exp Bot, 66(17): p. 5257-69.</p>
</li>
<li>
<p>Krogh, A., et al. (2001). Predicting transmembrane protein<br>
topology with a hidden Markov model: application to complete<br>
genomes. J Mol Biol, 305(3): p. 567-80.</p>
</li>
<li>
<p>Li B., &amp; Dewey CN. (2011). RSEM: accurate transcript quantification<br>
from RNA-Seq data with or without a reference genome. BMC<br>
Bioinformatics, 12, 323.</p>
</li>
<li>
<p>Zhou P, Silverstein KA, Gao L, Walton JD, Nallu S, Guhlin J and<br>
Young ND (2013). Detecting small plant peptides using SPADA<br>
(Small Peptide Alignment Discovery Application). BMC<br>
Bioinformatics 14: 335.</p>
</li>
</ul>

