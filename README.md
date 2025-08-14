| File/Folder | Description |
|--|--|
| **Fig1** | **Code and data for constructing the Gephi-based copangraph visualization (Figure 1)** |
| Fig1/code/nodes_of_interest.py: | Constructs a csv classifying nodes as "in METABAT2 bin" or not. Reduces an input gfa to include only nodes <= 5 nodes away from a binned node. |
| Fig1/code/gfa_to_gexf.py: | Converts a gfa file to a gexf file used as input to gephi. Replaces each node of length L by a chain ceil(L/2500) nodes. |	
| Fig1/code/metabat.sh: | Bins the MEGAHIT contigs with MetaBAT2. | 
| Fig1/code/gtdbtk.sh: | Taxonomically classifies the MetaBAT2 bins. |
| Fig1/data/raw_graph.gfa.gz: | Copangraph of the gut metagenomic samples |
| Fig1/data/raw_graph.ncolor.csv.gz: | csv of the node x sample occurrence matrix output by copangraph. |
| Fig1/data/sample_map.csv: | List of samples in the copangraph graph. |
| Fig1/data/persistence_table.csv: | Whether the sample has a persistence or clearance label. |
| Fig1/data/metabat2bins.tar: | MetaBAT2 bins in format <sample_name>.bin.fasta.gz |
| Fig1/data/binned_nodes.csv.gz: | csv produced by nodes_of_interest.py listing whether each node contains intervals from binned or only unbinned contigs. |
| Fig1/data/reduced.gfa.gz: | Reduced gfa produced by nodes_of_interest.py. |
| Fig1/data/gtdbtk.bac120.summary.tsv: | taxonomic classification of each contig by GTDB-Tk. |
| Fig1/data/Figure2.gephi: | Gephi file specifying the layout of the final graph visualized in Figure 1. |
|--|--|
| **Fig3bd** | **Code and data to construct connectivity and coverage F-score plots on simulated samples (Figures 3a,b).** |
| `Fig3bd/run_simulated_eval.py/` | Constructs coverage and connectivity TP, FP, FN counts for a set of input sequence graphs against a reference metagenome. Input graphs and reference metagenome are specified by a yaml in `./data/Fig3bd/eval_yamls`. Writes the results to a csv, `<sample>_megahit_<depth>_quality.csv` (found in evaluation), where `<sample>` is the sample name and `<depth>` is the number of read pairs (in millions) used to construct the graphs. |
| `Fig3bd/plot_simulation_results.py/` | Takes the `<sample>_megahit_<depth>_quality.csv`'s as input, and calcualtes connectivity/coverage precision, recall, and F-score, and plots them by read depth depth. |
| `AssemblyParser.py, parse_seq.py, evaluate_assembly.py/` | Auxiliary files used by run_simulated_eval.py. |
| `./data/Fig3bd/MEGAHIT_contigs/` | MEGAHIT assemblies (final.contigs.fa) of the simulated samples at each depth. Format `<sample>_<depth>M_final.contigs.fa`, where `<sample>` is the sample name and `<depth>` is the numbef of read pairs (in millions). |
| `./data/Fig3bd/MEGAHIT_graph/` | megahit_toolkit constructed assembly graphs (k141.fastg) of the simulated samples at each depth. Format `<sample>_<depth>M_k141.fastg`. |
| `./data/Fig3db/copangraph/` | copangraphs constructed from the simulated samples at each depth. Format `<sample>_<depth>.gfa`.|
| `./data/Fig3db/MetaCarvel/` | MetaCarvel constructed scaffold graphs of the simulated date at each depth. Format `./<sample>_<depth>/scaffold_graph.seq.gfa`. |
| `./data/camisim/<sample>_out/` | Each <sample>_out directory contains CAMISIM3-simulated simulated metagenomic data: sequencing reads, reference genomes, abundance profile, etc. Relative to `<sample>_out`, the sequencing reads are located at `<sample>_out/data/reads`, and the reference metagenome used as a ground truth by run_simulated_eval.py to calculate TP/FP/FN values is located at `<sample>_out/genomes/<sample>_reference_metagenome.fa`. |
| `./data/camisim/biom_files/` | Contains `<sample>.biom` files constructed from MetaPhlAn4 files `<sample>.species.mpa` of real metagenomes. `<sample>.biom` files are used to by camisim to construct the simulated metagenomes. |
| `./data/camisim/default_config.ini` | configuration file used by CAMISIM to construct simulated metagenomes; modified default config to set the number of simulated reads to ~ 80 million (by setting Gbp sequenced to 24). |
| `./evaluation/` | Contains the `<sample>_megahit_<depth>_quality.csv` files output by `run_simulated_eval.py`. Contains pdfs used in manuscript (`megahit_cov_F-score_labels.pdf` is Figure 3b, `megahit_cnx_F-score_labels.pdf` is Figure 3d). Also contains one-sided wilcoxon paired tests, comparing copangraph against other tools (`megahit_cov_F-score_wcox.csv, megahit_cnx_F-score_wcox.csv`, used for Figure 3b,d respectively). |
|--|--|
| **Fig3efgh** | **Code and data to construct F-score plots on real samples (Figures 3e-h).** |
| `Fig3efgh/run_real_eval.py` | Constructs coverage and connectivity TP, FP, FN counts of a set of input sequence graphs against a reference metagenome. Input graphs and reference metagenomes are specified by a yaml. For Figures 3e,f, which compare multi-sample graphs, these yamls are in `data/Fig3efgh/coassembly_eval_yamls`. For Figures 3g,h these yamls are in `data/Fig3efgh/assembly_eval_yamls`. For Figures 3e,f, outputs csvs `<n_samples>_sample_coassembly_analysis_quality.csv` (found in `coassembly_evaluation`) where `<n_samples` is either 3, 6, or 9 depending on the coassembly size. For Figures 3g,h, outputs csvs `<sample>_quality.csv` (found in `assembly_evaluation`) where `<sample>` is the metagenomic sample name. |
| `Fig3efgh/plot_real_data_results.py/` | Takes the `*_quality.csv` in `coassembly_evaluation` or `assembly_evaluation` and a flag `assembly` or `coassembly` and calculates connectivity/coverage precision, recall, and F-score plots (F-score plots used in Figures 3efgh). |
| `AssemblyParser.py, parse_seq.py, evaluate_assembly.py` | Auxiliary files used by run_real_eval.py. |
| `coasssembly_sample_names/` | List of samples names for the 3, 6, and 9 sample coassemblies. |
| `coassembly_MEGAHIT_contigs/` | MEGAHIT coassemblies of the 3, 6, and 9 pooled samples. |
| `coassembly_MEGAHIT_graph/` | megahit_toolkit constructed coassembly graphs of the 3, 6, 9, pooled samples. |
| `coassembly_MetaCarvel/` | MetaCarvel scaffold graphs of the 3, 6 pooled samples (based on the 6, 9 MEGAHIT coassemblies). |
| `coassembly_copangraph/` | Multi-sample copangraphs of the 3, 6, and 9 samples. |
| `long_read_mag_construction/flye_long_read_assemblies` | Flye assemblies of the long reads from each of the 9 Illumina HiFi metagenome samples. |
| `long_read_mag_construction/contigs_gt_1M` | Contigs from long-read assemblies with length > 1 million base pairs, each in a separate fasta file. |
| `long_read_mag_construction/checkm_out/low_contam_bins` | Contigs > 1 million base pairs (each considered long-read bin), with contamination below 5% as calculated by CheckM. |
| `long_read_mag_construction/dreped_bins_99/dereplicated_genomes` | Long-read bins dereplicated at 99% ANI. Each is considered a MAG, and used to construct the real data reference metagenomes. Format is `GHERIG_<sample>_contig_<n>.fa`, where `<sample>` is the name of the sample, and `<n>` is the contig name. |
| `coassembly_reference_metagenomes/` | Reference metagenomes for the 3, 6, 9 sample analysis (Figures 3e,f). Constructed by concatenating the MAGs from `long_read_mag_construction/dreped_bins_99/dereplicated_genomes` of the relevant samples. |
| `coassembly_eval_yamls/` | yamls input to `run_real_eval.py` specifying the graphs and reference metagenome for each coverage / connectivity evaluation used in Figures 3e,f. |
| `coassembly_evaluation/` | `<n_sample>_sample_coassembly_analysis_quality.csv` files output by `run_real_eval.py`. Contains pdfs used in manuscript (`coassembly_analysis_co_cov_F-score_labels.pdf`, and `coassembly_analysis_co_cnx_F-score_labels.pdf`, Figure 3e, 3f respectively). |
| `assembly_MEGAHIT_contigs/` | MEGAHIT assemblies (i.e, constructed from a single metagenome) for each of the 9 short-read metagenomic samples. |
| `assembly_MEGAHIT_graphs/` | megahit_toolkit assembly graphs for each of the 9 real short-read metagenomic samples. |
| `single_sample_MetaCarvel/` | MetaCarvel scafold graphs for each of the 9 short-read metagenomc samples. |
| `extracted_copangraph/` | Copangraphs, one for each of the 9 short-read samples, extracted from the nine-sample copangraph (`coassembly_copangraph/9_sample_graph.gfa`). |
| `single_sample_reference_metagenomes/` | The reference metagenome from each sample, constructed by concatenating MAGs in ` long_read_mag_construction/dreped_bins_99/dereplicated_genomes` of the same sample. |
| `assembly_eval_yamls/` | yamls input to `run_real_eval.py` specifying the graphs and reference metagenome for each coverage / connnectivity evaluation used in Figures 3g,h. |
| `assembly_evaluation/` | `<sample>_quality.csv` files output by `run_real_eval.py`. Contains pdfs used in manuscript (`ss_cov_F-score_labels.pdf`, and `ss_cnx_F-score_labels.pdf`, Figure 3g, 3h, respectively). |
|--|--|
| **Fig4a** | **Code and data to plot number of copangraph nodes as a function of the homology threshold and mean ANI of a population of related genomes. ** |
|--|--|
| **Fig4b** | **Code and data to plot taxonomic richness of nodes as a function of homology threshold** |
| `construct_taxonomy_curve.py` | Code to plot the (log) of the taxonomic richness as a function of homology threshold for the simulated samples. Takes as input: 1) the `<sample>_out` CAMISIM simulated directories, which contain taxonomic information required to calculate richness; 2) BLAST records of all contigs agains all reference genomes for each metagenomic samples, which is requires to calculate which contigs derive from which genomes; 3) Copangraph fasta files, describing which contig intervals belong to which nodes. |
| `<sample>_20M.blast.txt` | Tab-formated blast records of the MEGAHIT contigs from each `<sample>` against `coassembly_reference.fa`, which is the non-redundant set of all reference genomes from all reference metagenomes. |
| `graph_sd<threshold>.fasta` | Copangraph .fasta files output from constructing a multi-sample copangraph on the ten simulated metagenomic samples (used in Figure 3b,d), when setting the homology threshold to `<threshold>`. A particular threshold value `<threshold>` equates to collapsing sequences with an identity >= `100 - 0.<threshold>`%. |
| `all_alns.csv` | All BLAST records passing the filtering thresholds of >= 98% identity and alignment length being 98%-102% of the contig length. |
| `alpha_divs.csv` | Output file used to plot Figure 4b. Lists the taxonomic richeness of each node at each homology threshold. |
| `coassembly_reference.fa` | (non-redundant) concatenation of all reference genomes used in all ten simulated metagenomes (used to BLAST contigs against) |
| `taxonomic_curve.pdf` | Figure output by `construct_taxonomic_curve.py`; Figure 4b in manuscript. |
|--|--|
| **Fig4cdef** | **Code and data to show the suitability of copangraphs (vs coassembly graphs) for comparative metagenomic analysis**. |
| `mixing_cals.py` | Constructs intermediate data tables for plotting. Takes as input nine-sample coassembly reference metagenome (`coassembly_reference_metagenomes/9_sample_reference.fasta`) and a sequence graph, calculates n50, graph complexity (num. edges), coverage F-score and multi-genomic node proportion. |
| `gherig_analysis.ipynb` | Constructs plots based intermediate data tables used for Figures 4c-f. |
| `MEGAHIT_graphs/k<kmer>.fastg ` | MEGAHIT coassembly graphs constructed with megahit_toolkit on the nine real samples, where `<kmer>` specifies the kmer size used by MEGAHIT to construct the contigs. |
| `copangraphs/gherig_copan_sd_0.<threshold>.gfa` | multi-sample copangraphs constucted on the nine real samples, each with a different homology threshold, `<threshold>`. A particular `<threshold>` equates to collapsing sequences with an identity >= `100 - 0.<threshold>`%. 
|--|--|
| **Fig4ghij** | **Code and data to compare the ability of copangraphs and coassembly graphs to resolve between-sample strain variability.** |
| `mixing_cals_bacmeta.py` | Constructs intermediate data tables for plotting. Takes as input the population of genomes simulated under an evolutionary model (as reference) anda graph. |
| `bacmeta_analysis.ipynb` | Construct plots based on intermediate data tables for Figures 4ghij. |
| `copangraphs/`| Copangraphs constructed from the population of 25 genomes at different homology thresholds. Filename structure is `<rep>_5000_0.<threshold>_copangraph.gfa` where a particular `<threshold>` equates to collapsing sequences with an identity >= `100 - 0.<threshold>`%, `<rep>` is the replicate. |
| `MEGAHIT_graphs/` | megahit_toolkit coassembly graphs constructed from the population of 25 genomes at different kmer sizes. Filename structure is `<rep>_5000_k<kmer>.fastg` where `<rep>` is the replicate and `<kmer>` is the kmer size used to construct the graph. |
| `fastq/` | Simulated sequencing reads constructed from the genome populations. |
| `genomes/` | Fasta files, `<rep>_Gen5000_reference.fasta`, of the 25 genomes in each population, where `<rep>` is a replicate population. |
|--|--|
| **Fig5** | **Code and data to quantify resource utilization** |
| `plot_resource.ipynb` | Plotting code to construct Figure 5a,b. Takes as input `resource_dat.csv`. |
| `resource_dat.csv` | Memory and runtime to process each sample with MEGAHIT, metaSPAdes, or copangraph |
|--|--|
|**Fig6b** | **Code and data for VRE persistence prediction**|
|--|--|
|**Fig6cegh** | **Code and data for comparative metagenomics with copangraph**|
|**Fig6df**| **Code and data for VRE-isolate - copangraph bubble analysis** |


