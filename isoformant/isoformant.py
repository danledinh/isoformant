import os, glob, random, subprocess, warnings, collections
import numpy as np
import pandas as pd
from plotnine import *
import plotnine
import khmer
import scanpy as sc
import anndata as ad
from tqdm import tqdm
import pysam
import genomeview
from IPython.display import display
from scipy import sparse
from scipy.stats import chisquare
from spoa import poa

# Set a seed value
seed_value= 100
os.environ['PYTHONHASHSEED']=str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)

def shuffle_csv(incsv):
    """Shuffle CSV records.
    Performs task in-place.

    :param incsv: CSV file path
    :type incsv: str
    """

    outfa = incsv + '.tmp'
    with open(outfa, 'w') as outf:
        with open(incsv) as f:
            data = f.readlines()
            data_idx = [x for x in range(0, len(data), 1)]
            random.Random(seed_value).shuffle(data_idx)
            for i in data_idx:
                outf.write(data[i])
    os.rename(outfa, incsv)

def process_bam(bamfile, chrom, start, end, 
                suffix, outcsv, outbam, 
                qual_cutoff = 10, len_cutoff = 300):
    """Performs a series of tasks on BAM file.
    Filter BAM reads: (1) not secondary/supplemental read, (2) min mean base quality, (3) min read length.
    Write to disk in BAM format.
    Write to disk in CSV format. Randomly shuffled.

    :param bamfile: Path to sorted BAM file (.bai index file required in same directory)
    :type bamfile: str
    :param chrom: Reference sequence name of ROI
    :type chrom: str
    :param start: 1-based start coordinate of ROI
    :type start: int
    :param end: 1-based end coordinate of ROI
    :type end: int
    :param suffix: Identifier to be appended to read names
    :type suffix: str
    :param outcsv: Path to CSV output file
    :type outcsv: str
    :param outbam: Path to BAM output file
    :type outbam: str
    :param qual_cutoff: Minimum mean base quality, defaults to 10
    :type qual_cutoff: float, optional
    :param len_cutoff: Minimum read length, defaults to 300
    :type len_cutoff: int, optional
    
    :return: Number of records passing filters
    :rtype: int
    """
    print('>>Task: Filter BAM file')
   
    # initiate counter
    keep_n = 0
    
    # parse input BAM
    with pysam.AlignmentFile(bamfile, 'rb') as samfile:
        # open output BAM
        with pysam.AlignmentFile(outbam, 'wb', template=samfile) as writebam:
            # open output CSV
            with open(outcsv, 'w') as outf:
                # fetch ROI records
                for idx, s in tqdm(enumerate(samfile.fetch(chrom, start, end))):
                    # filter out secondary/supplemental reads
                    if not s.is_secondary and not s.is_supplementary:
                        if np.mean(s.query_qualities) >= qual_cutoff and s.query_length >= len_cutoff:
                            # write to disk 
                            ## CSV
                            outf.write(f'{suffix}\t{s.query_name}\t{s.query_sequence}\n')
                            ## BAM
                            writebam.write(s)
                            # update count
                            keep_n += 1
                            
    # shuffle FA                        
    shuffle_csv(outcsv)
    
    # print info
    print('\t' + 'BAM = ' + bamfile)
    print('\t' + '# of records in ROI:', idx+1)
    print('\t' + '# passing:', keep_n)
    
    return keep_n

def preprocess_reads(bams_list, chrom, start, end, dest_dir, max_reads=1000, ncore=1, qual_cutoff = 10, len_cutoff = 300):
    """Module to process a list of BAM files.
    Coordinates 'process_bam' task: extract passing reads.
    Creates merged BAM file and writes to disk.
    Creates merged min-depth balanced CSV file and writes to disk.

    :param bams_list: List of sorted BAM file paths (.bai required in same directory)
    :type bams_list: list
    :param chrom: Reference sequence name of ROI
    :type chrom: str
    :param start: 1-based start coordinate of ROI
    :type start: int
    :param end: 1-based end coordinate of ROI
    :type end: int
    :param dest_dir: Path to directory for output files
    :type dest_dir: str
    :param max_reads: Maximum possible reads considered, defaults to 1000
    :type max_reads: int, optional
    :param ncore: Number of cores, defaults to 1
    :type ncore: int, optional
    :param qual_cutoff: Minimum mean base quality, defaults to 10
    :type qual_cutoff: float, optional
    :param len_cutoff: Minimum read length, defaults to 300
    :type len_cutoff: int, optional
    
    :return: Path to merged CSV file
    :rtype: str
    :return: Path to sorted merged BAM file
    :rtype: str     
    """
    
    print('>Module: Process BAM files')
    
    # check for dest dir. create if not extant
    if os.path.exists(dest_dir) == False:
        os.mkdir(dest_dir)
    
    # Parse list of BAMs
    sample_dict = {}
    sample_n_l = []
    for sortedbamfile in bams_list:
        sortedbamindex = sortedbamfile + '.bai'
        
        # extract filename
        prefix = sortedbamfile.split('/')[-1]
        
        # generate output filenames
        outcsv = f'{dest_dir}/{prefix}.{chrom}:{start}-{end}.csv'
        outbam = f'{dest_dir}/{prefix}.{chrom}:{start}-{end}.bam'
        sortedoutbam = f'{dest_dir}/{prefix}.{chrom}:{start}-{end}.sorted.bam'
        
        # process_bam task + sort + index
        keep_n = process_bam(sortedbamfile, chrom, start, end, prefix, outcsv, outbam, qual_cutoff, len_cutoff)
        sample_n_l.append(keep_n)
        sample_dict[prefix] = [sortedoutbam, outcsv]
        pysam.sort('-@', str(ncore), '-o', sortedoutbam, outbam)
        pysam.index(sortedoutbam)
        
    # merge roi bam files + sort + index
    roi_list = [f[0] for key, f in sample_dict.items()]
    combined_bam = f'{dest_dir}/combined.bam'
    sortedcombined_bam = f'{dest_dir}/combined.sorted.bam'
    merge_parameters = ['-f', combined_bam] + roi_list
    pysam.merge(*merge_parameters)
    pysam.sort('-@', str(ncore), '-o', sortedcombined_bam, combined_bam)
    pysam.index(sortedcombined_bam)

    # merge csv files (balanced depth via random downsample)
    min_reads = np.array(sample_n_l).min()
    ## enforce max number of reads
    if min_reads > max_reads:
        min_reads = max_reads
    combined_csv = f'{dest_dir}/combined.csv'
    with open(combined_csv, 'w') as outfile:
        for key, f in sample_dict.items():
            with open(f[1], 'r') as infile:
                for n,line in enumerate(infile):
                    if n < min_reads:
                        outfile.write(line)
                    else:
                        break
                        
    return combined_csv, sortedcombined_bam
    
def kmerize(combined_csv, ksize=5):
    """Ingests CSV file and produces k-mer frequency anndata object.

    :param combined_csv: Path to CSV file
    :type combined_csv: str
    :param ksize: k-mer size, defaults to 5
    :type ksize: int, optional
    
    :return: k-mer frequencies as anndata object
    :rtype: anndata object
    """
    print('>Module: Read k-merization')
    print('>>Task: k-merize CSV reads')      
    print('\t' + f'k = {ksize}')
    
    # determine width of k-mer frequency table
    nkmers = 4**ksize
    
    # count number of seqs in CSV
    with open(combined_csv, 'r') as incsv:
        for nrow, row in enumerate(incsv):
            pass
    entry_count = nrow+1
    
    # create placeholders for IDs and counts
    kmer_arr = np.zeros((nkmers, entry_count), dtype = 'int')
    arr_idx = 0
    
    # parse CSV
    sample_id_l = []
    read_name_l = []
    with open(combined_csv, 'r') as incsv:
        for row in tqdm(incsv):
            row_elems = row.strip().split('\t')
            sample_id = row_elems[0]
            read_name = row_elems[1]
            seq = row_elems[2]
            
            # Initialize countgraph
            ktable = khmer.Countgraph(ksize, nkmers + 10, 1)

            # count all k-mers in the given string
            ktable.consume(seq)

            # capture full kmer counts
            k_n_list = [ktable.get(i) for i in range(nkmers)]

            # update kmer count arr
            kmer_arr[:,arr_idx] = k_n_list 

            # update arr pointer
            arr_idx = arr_idx + 1
            
            # update lists
            sample_id_l.append(sample_id)
            read_name_l.append(read_name)
    
    # extract kmer strings
    kmer_vars = [ktable.reverse_hash(i) for i in range(nkmers)]
    
    # normalize to column sum
    kmer_arr = kmer_arr/kmer_arr.sum(0)
                    
    # create ad obj
    print('>>Task: Convert k-mer pandas dataframe to anndata object')
    adata_reads = ad.AnnData(X=kmer_arr).T
    adata_reads.X = sparse.csr_matrix(adata_reads.X)
    adata_reads.var_names = kmer_vars
    adata_reads.obs_names = [f'{x}' for x in range(len(adata_reads.obs))]
    
    # extract features
    adata_reads.obs['sample_id'] = sample_id_l
    adata_reads.obs['read_name'] = read_name_l
    
    return adata_reads

def pca_pipeline(adata_reads, groupby_list, components=['1,2'], n_comps=100, plot=False):
    """PCA pipeline performed in-place.
    Scale data and perform PCA.
    Option to print PC scatter plot.

    :param adata_reads: Read x k-mer data
    :type adata_reads: Anndata object
    :param groupby_list: List of feature names to color-code in plot
    :type groupby_list: list
    :param components: List of PC pairs to plot (e.g. ['1,2', ...]), defaults to ['1,2']
    :type components: list, optional
    :param n_comps: Number of PCs to compute, defaults to 100
    :type n_comps: int, optional
    :param plot: 'True' to plot, 'False' to pass, defaults to 'False'
    :type plot: bool, optional
    """
    print('>>Task: PCA')
    print('\t' + f'# of PCs = {n_comps}')
    
    sc.pp.scale(adata_reads)
    sc.tl.pca(adata_reads,
              n_comps=n_comps,
             )
    if plot == True:
        sc.pl.pca(adata_reads,
                  components=components,
                  color = groupby_list,
                  ncols=1
                 )
            
def umap_pipeline(input_adata, 
                  groupby_list,
                  n_pcs = 3,
                  n_neighbors = 15,
                  min_dist = 0.5,
                  resolution = 1,
                  components=['1,2'],
                  plot = False,):
    """UMAP pipeline performed in-place.
    Build KNN graph, compute UMAP, and cluster reads using Leiden algorithm. 
    Option to print UMAP scatter plot.

    :param input_adata: Read x k-mer data (processed by `pca_pipeline`)
    :type input_adata: Anndata object
    :param groupby_list: List of feature names to color-code in plot
    :type groupby_list: list
    :param n_pcs: Number of PCs for KNN graph, defaults to 3
    :type n_pcs: int, optional
    :param n_neighbors: Number of neighbors per neighborhood, defaults to 15
    :type n_neighbors: int, optional
    :param min_dist: Minimum distance to neighbor, defaults to 0.5
    :type min_dist: float, optional
    :param resolution: Leiden clustering resolution, defulats to 1
    :type resolution: float, optional
    :param components: List of UMAP component pairs to plot (e.g. ['1,2', ...]), defaults to ['1,2']
    :type components: list, optional
    :param plot: 'True' to plot, 'False' to pass, defaults to 'False'
    :type plot: bool, optional
    """
    print('>>Task: UMAP')
    print('\t' + f'# of PCs = {n_pcs}')
    print('\t' + f'# of neighbors = {n_neighbors}')
    print('\t' + f'Minimum distance = {min_dist}')
    print('\t' + f'Resolution = {resolution}')    

    sc.pp.neighbors(input_adata, n_pcs=n_pcs, n_neighbors = n_neighbors, random_state=seed_value) 
    sc.tl.umap(input_adata, min_dist = min_dist, random_state=seed_value)
    sc.tl.leiden(input_adata, resolution=resolution)
    input_adata.obs['leiden'] = [x.zfill(3) for x in input_adata.obs['leiden']]
        
    if plot == True:
        # plot not reproducible. must restart kernal and run full script: https://github.com/theislab/scanpy/issues/2014
        sc.pl.umap(input_adata,
                  components=components,
                  color = groupby_list,
                  ncols=1
                 )

def create_viz_bams(input_adata, combined_bam, dest_dir, cluster_label='leiden', bam_n=10, ncore=1):
    """Create BAM files for visualizations.
    Returns a dictionary mapping cluster ID to BAM file.

    :param input_adata: Read x k-mer data (processed by `umap_pipeline`)
    :type input_adata: Anndata object
    :param combined_bam: Path to BAM file associated with processed reads
    :type combined_bam: str
    :param dest_dir: Path to directory for output files
    :type dest_dir: str
    :param cluster_label: Cluster label, defaults to 'leiden'
    :type cluster_label: str, optional
    :param bam_n: Maximum number of records per BAM file, defaults to 10
    :type bam_n: int, optional
    :param ncore: Number of cores, defulats to 1
    :type ncore: int, optional
    
    :return: Dictionary mapping cluster to BAM path
    :rtype: dict
    :return: Dictionary mapping cluster to read sequences
    :rtype: dict
    """
    print('>>Task: Create visualization BAM files')
    print('\t' + f'Maximum number of records per BAM = {bam_n}')
    
    # prepare dir for bams
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
    else:
        rm_list = glob.glob(f'{dest_dir}/*')
        for f in rm_list:
            os.remove(f)

    # get bam_n number of random read indices per cluster
    bam_df = input_adata.obs
    bam_df_slice = bam_df.sample(frac=1, random_state=seed_value).groupby(cluster_label).tail(bam_n)
    bam_write_dict = {}
    for cluster, qname in zip(bam_df_slice[cluster_label],bam_df_slice['read_name']):
        bam_write_dict[qname]=cluster

    # write bam records if matching index
    seq_dict = {}
    with pysam.AlignmentFile(combined_bam, 'rb') as samfile:
        sorted_clusters = sorted(list(set(bam_write_dict.values())))
        files = [pysam.AlignmentFile(f'{dest_dir}/clust{i}.bam', 'wb', template=samfile) for i in sorted_clusters]

        for s in samfile:
            cluster_id = bam_write_dict.get(s.query_name, None)
            if cluster_id is not None:
                files[sorted_clusters.index(cluster_id)].write(s)
                if seq_dict.get(cluster_id, None) is not None:
                    seq_dict[cluster_id] = seq_dict.get(cluster_id) + [s.query_sequence]
                else:
                    seq_dict[cluster_id] = [s.query_sequence]

        for f in files:
            f.close()

    #sort and index bam
    bam_dict = {}
    for cluster_id in sorted_clusters:
        outbam = f'{dest_dir}/clust{cluster_id}.bam'
        sortedbam = f'{dest_dir}/clust{cluster_id}.sorted.bam'
        bam_dict[cluster_id] = sortedbam
        pysam.sort('-o',sortedbam,outbam)            
        pysam.index(sortedbam)
        
    return bam_dict, seq_dict

def roi_fa(ref_fa, chrom, out_fa):
    """Trim FASTA to reference sequence of interest.
    
    :type ref_fa: str
    :param cnsfa: Path to input query FASTA file
    :param chrom: Reference sequence name of ROI
    :type chrom: str
    :param out_fa: Path to output file
    :type out_fa: str
    """
    
    with pysam.Fastafile(ref_fa) as infa:
        roi = infa.fetch(chrom)
        with open(out_fa, 'w') as outf:
            outf.write(f'>{chrom}\n')
            outf.write(roi + '\n')
    
def minimap2_launcher(minimap2_path, ref_fa, cnsfa, cnsbam, sirv=False, messages=False):
    """Launch minimap2 alignment.

    :param minimap2_path: Path to minimap2 executable
    :type minimap2_path: str
    :param ref_fa: Path to reference FASTA file
    :type ref_fa: str
    :param cnsfa: Path to input query FASTA file
    :type cnsfa: str
    :param cnsbam: Path to output query BAM file
    :type cnsbam: str
    :param sirv: True to benchmark SIRVs, defaults to False
    :type sirv: bool, optional
    :param messages: True to print STDOUT/STDERR messages, defaults to False
    :type messages: bool, optional
    """
    print('>>Task: Align sequences')
    
    # minimap2 path can be found inside conda env with "which minimap2" command
    command_l = [minimap2_path, 
                '-ax', 'splice',
                '-o', cnsbam,
                ref_fa,
                cnsfa
               ]
    
    # option for benchmarking SIRVs
    if sirv == True:
        command_l.insert(3, '--splice-flank=no')
   
    process = subprocess.Popen(command_l,
                               stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE,
                               universal_newlines=True,
                              )
    stdout, stderr = process.communicate()

    if messages == True:
        print(command_l)
        print('stdout:', stdout)
        print('stderr:', stderr)
    
def cns(minimap2_path, seq_dict, dest_dir, ref_fa, cluster_label = 'leiden', ncore=1, sirv=False, messages=False):
    """Consensus calling module.

    :param minimap2_path: Path to minimap2 executable
    :type minimap2_path: str
    :param seq_dict: Dictionary mapping cluster to read sequences
    :type seq_dict: dict
    :param dest_dir: Path to directory for output files
    :type dest_dir: str     
    :param ref_fa: Path to reference FASTA file
    :type ref_fa: str
    :param cluster_label: Cluster label, defaults to 'leiden'
    :type cluster_label: str, optional
    :param ncore: Number of cores, defaults to 1
    :type ncore: int, optional
    :param sirv: True to benchmark SIRVs, defaults to False
    :type sirv: bool, optional
    :param messages: True to print STDOUT/STDERR messages, defaults to False
    :type messages: bool, optional
    
    :return: Path to consensus sequence set BAM file
    :rtype: str
    :return: Dictionary mapping cluster-specific consensus sequence to BAM path
    :rtype: dict
    """
    print('>Module: Consensus calling')
    
    # prepare directory for outputs
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
    else:
        rm_list = glob.glob(f'{dest_dir}/*')
        for f in rm_list:
            os.remove(f)
    
    # POA step
    print('>>Task: POA step')
    cns_l = []
    clust_l = []
    for cluster_id in tqdm(seq_dict.keys()):
        seq_l = seq_dict.get(cluster_id)
        consensus, _ = poa(seq_l)
        cns_l.append(consensus)
        clust_l.append(cluster_id)
    
    # write to FASTA
    cnsfa = f'{dest_dir}/cns.fa'
    with open(cnsfa, 'w') as faf:
        for qname, qseq in zip(clust_l, cns_l):
            faf.write(f'>{qname}' + '\n')
            faf.write(f'{qseq}' + '\n')   
            
    # map CNS to ref
#     print('>>Task: Align consensus sequences')
    cnsbam = cnsfa + '.bam'
    sorted_cnsbam = cnsfa + '.sorted.bam'
    minimap2_launcher(minimap2_path, ref_fa, cnsfa, cnsbam, sirv=sirv, messages=messages)
    pysam.sort('-@', str(ncore), '-o', sorted_cnsbam, cnsbam)
    pysam.index(sorted_cnsbam)
    
    # Split CNS set to individual cluster-specific BAMs
    cns_dict = {}
    with pysam.AlignmentFile(sorted_cnsbam, 'rb') as samfile:
        for s in samfile:
            # split by qname and write new bam
            qname = s.qname
            outbam = f'{dest_dir}/split_cns_{qname}.bam'
            cns_dict[qname] = outbam
            with pysam.AlignmentFile(outbam, 'wb', template=samfile) as writebam:
                writebam.write(s)
            pysam.index(outbam)
    
    return sorted_cnsbam, cns_dict

def plot_highlight_COI(input_adata, COI, cluster_label = 'leiden'):
    """Highlight cluster of interest in UMAP plot.

    :param input_adata: Read x k-mer data (processed by `umap_pipeline`)
    :type input_adata: Anndata object
    :param COI: Cluster label
    :type COI: str
    :param cluster_label: Cluster label, defaults to 'leiden'
    :type cluster_label: str, optional
    """
    
    umap_df = pd.DataFrame(input_adata.obsm['X_umap'])
    umap_df.columns = ['xvar','yvar']
    umap_df[cluster_label] = input_adata.obs[cluster_label].values
    umap_df['COI'] = [x == COI for x in umap_df[cluster_label]]
        
    plotnine.options.figure_size = (4,4)
    plot = (ggplot(umap_df)
            + theme_bw()
            + theme(aspect_ratio = 9/16,
                    axis_text = element_blank()
                   )
            + geom_point(aes('xvar','yvar',color='COI'))
            + ggtitle(f'cluster: {COI}')
            + labs(x='',y='')
           )
    print(plot)
    
def balanced_chisq_test(adata, cond_list=[], cond_label='sample_id', cluster_label='leiden'):
    """Chi-squared goodness-of-fit test.
    Assumes class balance. 

    :param input_adata: Read x k-mer data (processed by `umap_pipeline`)
    :type input_adata: Anndata object
    :param cond_list: List of condition labels to be considered, defaults to []
    :type cond_list: list, optional
    :param cond_label: Condition label, defaults to 'sample_id'
    :type cond_label: str, optional
    :param cluster_label: Cluster label, defaults to 'leiden'
    :type cluster_label: str, optional
    
    :return: Test results
    :rtype: pandas dataframe
    """
    print('>>Task: Chi-sqaured test')
    if len(cond_list) == 0:
        cond_list = sorted(list(set(adata.obs[cond_label])))

    print('\t' + f'Considering {cond_list} from class = "{cond_label}"')

    cluster_l = []
    pval_l = []
    for cluster in set(adata.obs[cluster_label]):
        cluster_slice = adata.obs[adata.obs[cluster_label] == cluster]
        cluster_slice = cluster_slice[[x in cond_list for x in cluster_slice[cond_label]]]
        obs = np.array(cluster_slice[cond_label].value_counts())
        if np.any(obs<5):
            warnings.warn('Chi-squared test overrejects when any observed frequency less than 5')
        _, pval = chisquare(obs)
        cluster_l.append(cluster)
        pval_l.append(pval)

    results_df = pd.DataFrame({cluster_label:cluster_l,
                               'pval':pval_l, 
                              })
    results_df['pval_adj'] = results_df['pval']*len(results_df)
    results_df['pval_adj'] = [1 if x>=1 else x for x in results_df['pval_adj']]
    results_df = results_df.sort_values('pval_adj')
    
    return results_df

def plot_cluster_occupancy(adata, cond_list=[], cond_label='sample_id', cluster_label='leiden', fig_size=(6,6), subplot_hjust=0.3, save=None):
    """Plots cluster occupancy by condition label.

    :param input_adata: Read x k-mer data (processed by `umap_pipeline`)
    :type input_adata: Anndata object
    :param cond_list: List of condition labels to be considered, defaults to []
    :type cond_list: list, optional
    :param cond_label: Condition label, defaults to 'sample_id'
    :type cond_label: str, optional
    :param cluster_label: Cluster label, defaults to 'leiden'
    :type cluster_label: str, optional
    :param fig_size: Figure size, defaults to (10,6)
    :type fig_size: tuple, optional
    :param subplot_hjust: Subplot horizontal spacing, defaults to 0.3
    :type fig_size: float, optional
    :param save: Path to save plot, defaults to None
    :type save: None, optional
    
    :return: Frequency table
    :rtype: pandas dataframe
    """
    if len(cond_list) == 0:
        cond_list = sorted(list(set(adata.obs[cond_label])))
    
    adata_slice = adata.obs[[x in cond_list for x in adata.obs[cond_label]]]
    prop_df = (pd.DataFrame(adata_slice
                           .groupby(cluster_label)[cond_label]
                           .value_counts(normalize=False)
                          )
               .rename(columns={cond_label:'frequency'})
               .reset_index()
               .sort_values('frequency', ascending=False)
              )
    prop_df['group_rank'] = (prop_df
                             .groupby(cluster_label)
                             .rank(method='first')
                             .astype(int)
                             .astype('category')
                            )

    plotnine.options.figure_size=fig_size
    plot = (ggplot(prop_df)
            + theme_bw()
            + theme(subplots_adjust={'wspace': subplot_hjust},
                    axis_text_x=element_blank()
                   )
            + geom_bar(aes('group_rank', 'frequency', fill=cond_label), stat='identity', position='dodge')
            + facet_wrap(f'~{cluster_label}', scales = 'free_y')
            + labs(x='',y='frequency')
           )
    print(plot)
    
    if save:
        ggsave(plot, save)
    return prop_df.sort_values([cluster_label, cond_label])

def plot_tracks(track_dict, chrom, start, end, ref_fa, left_margin = 0, save=None):
    """Given dictionary, plot reads tracks.

    :param track_dict: Dictionary of track labels mapped to BAM files 
    :type input_adata: dict
    :param chrom: Reference sequence name of ROI
    :type chrom: str
    :param start: 1-based start coordinate of ROI
    :type start: int
    :param end: 1-based end coordinate of ROI
    :type end: int
    :param ref_fa: Path to reference FASTA file
    :type ref_fa: str
    :param left_margin: Left margin in genomic coodinates, defaults to 0
    :type left_margin: int, optional
    :param save: Path to save plot, defaults to None
    :type save: None, optional
    """
    
    doc = genomeview.visualize_data(track_dict, chrom, int(start)-int(left_margin), int(end), ref_fa)
    for key in track_dict.keys():
        doc.get_tracks(key)[0].color_fn = lambda x: "darkgray"
        doc.get_tracks(key)[0].nuc_colors = {"A":"blue", "C":"blue", "G":"blue", "T":"blue", "N":"yellow"}
        doc.get_tracks(key)[0].insertion_color = "red"
        doc.get_tracks(key)[0].clipping_color = "cyan"
        doc.get_tracks(key)[0].deletion_color = "green"
        doc.get_tracks(key)[0].insertion_color = "red"
        doc.get_tracks(key)[0].min_indel_size = 5
    display(doc)
    
    if save:
        genomeview.save(doc, save)