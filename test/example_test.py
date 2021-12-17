#!/usr/bin/env python3
"""
Run example and check output validity.
"""

def main():
    
    import subprocess, stat, sys, os
    subprocess = subprocess.Popen("which minimap2",
                              shell=True, 
                              stdout=subprocess.PIPE)
    minimap2_path = subprocess.stdout.read().decode("utf-8").strip()
    print(minimap2_path)
    st = os.stat(minimap2_path)
    os.chmod(minimap2_path, st.st_mode | stat.S_IEXEC)
    
    sys.path.append(os.path.abspath('./isoformant/'))
    from isoformant import preprocess_reads, kmerize, pca_pipeline, umap_pipeline, create_viz_bams, cns, balanced_chisq_test, plot_cluster_occupancy

    ncore=1 # CPUs 
    ref_fa = './example/SIRV_isoforms_multi-fasta_170612a.fasta' # reference genome 
    gene_track = './example/SIRV_isoforms_multi-fasta-annotation_C_170612a_corrected_norm.bed' # Optional, isoform annotations
    dest_dir = './example/outputs' # output directory 

    # ingest list of BAM files
    bams_list = [
        './example/SIRV1_cDNA_sorted.bam',
        './example/SIRV1_RNA_sorted.bam',
    ]

    # assign region-of-interest (ROI)
    chrom, start, end = ('SIRV1', 1, 12643)

    # process BAM
    combined_csv, combined_bam = preprocess_reads(bams_list, 
                                                  chrom, # ROI
                                                  start, # ROI
                                                  end, # ROI
                                                  dest_dir, 
                                                  max_reads=1000, # limits reads to 'max_reads' per BAM
                                                  ncore=ncore, # CPUs
                                                  qual_cutoff=15, # min. mean base quality
                                                  len_cutoff=500) # min. read length
    
    # k-merize sequences and create adata object
    adata = kmerize(combined_csv, 
                    ksize = 7 # k-mer size
                   )

    # pca pipeline: dimensionality reduction
    pca_pipeline(adata, groupby_list=['sample_id'], n_comps=100, plot=True)

    # umap pipeline: clustering
    umap_pipeline(adata, groupby_list=['leiden', 'sample_id'], 
              n_pcs = 20, # number of PCs
              n_neighbors = 10, 
              min_dist = 0.5, 
              resolution = 1, # clustering resolution
              plot = True
             )

    # Create BAMs for visualization
    bam_dict, seq_dict = create_viz_bams(adata, combined_bam, dest_dir + '/viz', cluster_label='leiden', bam_n=10, ncore=ncore)
    
    # Consensus calling: POA consensus and mapping with minimap2
    sorted_cnsbam, cns_dict = cns(minimap2_path,
                                    seq_dict, 
                                    dest_dir + '/cns', 
                                    ref_fa,
                                    cluster_label = 'leiden', 
                                    ncore = ncore,
                                    sirv=True
                                   )
    
    # Chi-squared goodness-of-fit test
    results_df = balanced_chisq_test(adata)
    
    # Plot cluster occupancy by condition label
    occupancy_df = plot_cluster_occupancy(adata, fig_size=(8,6))
    
if __name__ == "__main__":
    main()