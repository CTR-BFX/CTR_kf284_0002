# Damiano G. Barone




# Experiments

## CTR_kf284_0002

Implants vs Nerve Injury Mouse RNA-seq

## CTR_kf284_0003

RNA-seq comparing expression profile of peripheral nerve damage and microchannel device implantation repair model: control vs local delivery of MCC950 (NLRP3inh) vs  standard FBR treatment dexamethasone (Dex).


# Methods

## Bioinformatics Methods: Processing of RNA-Seq

Data were aligned to mm10 mouse genome (Ensembl Release GRCm38.p5) with STAR (v020201)(Dobin et al, 2012). Alignments and QC were processed using ClusterFlow (v0.5dev)(Ewels et al, 2017) pipelines (FASTQC, Trim_galore) and summarised using MultiQC (0.9.dev0)(Ewels et al, 2016). Gene quantification was determined with HTSeq-Counts (v0.6.1p1)(Anders et al, 2014). Additional quality control was performed with featureCounts (v 1.5.0-p2)(Liao et al, 2013) and qualimap (v2.2)(García-Alcalde et al, 2012). Differential gene expression was performed with DESeq2 package (v1.18.1, R v3.4.0)(Love et al, 2014) and with the same package read counts were normalised on the estimated size factors. Technical replicates run on separate lanes were collapsed using DESeq2.

UpSetR is an alternative for plotting sets of data to visualise overlaps as a more intuitive replacement for Euler/Venn Diagrams (Conway et al, 2017).

## Bioinformatics Methods: Deconvolving RNA-Seq

Samples were deconvoluted for the fractions of immune cell types present, using DeconRNASeq and a signature matrix that distinguish 25 mouse hematopoietic cell types, including six major cell types, granulocytes, B cells, T cells, natural killer cells, dendritic cells and mono/macrophages (see: Chen, Z. et al. (2017) Inference of immune cell composition on the expression profiles of mouse tissue. Scientific Reports 7, Article number: 40508).


## Software

Resource              | URL
--------------------- | --------------
GRCm38                | [Link](http://mar2016.archive.ensembl.org/index.html)
FastQC                | [Link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
Trim_galore           | [Link](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
STAR                  | [DOI](http://dx.doi.org/10.1093/bioinformatics/bts635)
HTSeq-counts          | [DOI](http://dx.doi.org/10.1093/bioinformatics/btu638)
Feature_counts        | [DOI](http://dx.doi.org/10.1093/bioinformatics/btt656)
Qualimap              | [DOI](http://dx.doi.org/10.1093/bioinformatics/bts503)
ClusterFlow           | [DOI](http://dx.doi.org/10.12688/f1000research.10335.2)
MultiQC               | [DOI](http://dx.doi.org/10.1093/bioinformatics/btw354)
DESeq2                | [DOI](http://dx.doi.org/10.1186/s13059-014-0550-8)
UpSetR                | [DOI](http://dx.doi.org/10.1093/bioinformatics/btx364)
DeconRNASeq           | [DOI](http://dx.doi.org/10.1093/bioinformatics/btt090)
Immune cell Reference | [DOI](http://dx.doi.org/10.1038/srep40508)


### CTR_kf284_0002 :: Experimental design: 

Samples are divided into groups: 
-control samples, 
-nerve injury samples at day 1, 4, 7, 14 and 28 
-implants samples at day 1, 4, 7, 14 and 28.

The interest of the project was focused on involvement of immune response genes, fibrosis and inflammatory response to nerve injury and implants. The hypothesis is that nerve injury will have greatest response up to day 7 and than come back to control state, while implants will trigger continues response, with greatest difference between implants and nerve injury at later days (14, 28).


### CTR_kf284_0003 :: Experimental design: 

Samples are divided into groups: 
-nerve implanted with untreated conduit device – NoDrug, 
-nerve implanted with conduit device impregnated with dexamethasone – Dex
-nerve implanted with conduit device impregnated with MCC950 – NLRP3inh)
Samples were explanted at day 28 post-surgery.









## Scripts

Script                      |   Description
----                        | ----
CTR_kf284_0002.R  [CTR_kf284_0002.R](R_scripts/CTR_kf284_0002.R) |  Analysis and QC of Implants vs Nerve Injury Mouse RNA-seq.
CTR_kf284_0003.R  [CTR_kf284_0003.R](R_scripts/CTR_kf284_0003.R) |  Analysis of MCC950 vs standard FBR treatment dexamethasone in a peripheral nerve damage and microchannel device implantation repair model.


## Contact

Contact *mn367@cam.ac.uk* for bioinformatics related queries.


## References

Anders,S., Pyl,P.T. and Huber,W. (2014) HTSeq--a Python framework to work with high-throughput sequencing data. Bioinformatics, 31, 166–169.

Chen, Z., Huang, A., Sun, J., Jiang, T., Qin, F. X.-F., & Wu, A. (2017). Inference of immune cell composition on the expression profiles of mouse tissue. Scientific Reports, 7, 40508.

Conway, J. R., Lex, A., & Gehlenborg, N. (2017). UpSetR: an R package for the visualization of intersecting sets and their properties. Bioinformatics, 33(18), 2938–2940

Dobin,A., Davis,C.A., Schlesinger,F., Drenkow,J., Zaleski,C., Jha,S., Batut,P., Chaisson,M. and Gingeras,T.R. (2012) STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29, 15–21.

Ewels,P., Magnusson,M., Lundin,S. and Käller,M. (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32, 3047–3048.

Ewels,P., Krueger,F., Käller,M. and Andrews,S. (2017) Cluster Flow: A user-friendly bioinformatics workflow tool. F1000Research, 5, 2824.

García-Alcalde,F., Okonechnikov,K., Carbonell,J., Cruz,L.M., Götz,S., Tarazona,S., Dopazo,J., Meyer,T.F. and Conesa,A. (2012) Qualimap: evaluating next-generation sequencing alignment data. Bioinformatics, 28, 2678–2679.

Gong, T., & Szustakowski, J. D. (2013). DeconRNASeq: a statistical framework for deconvolution of heterogeneous tissue samples based on mRNA-Seq data. Bioinformatics, 29(8), 1083–1085.

Liao,Y., Smyth,G.K. and Shi,W. (2013) featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30, 923–930.

Love,M.I., Huber,W. and Anders,S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15.
