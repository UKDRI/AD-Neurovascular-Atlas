#!/usr/bin/env python3

import pandas as pd
import numpy as np
from arboreto.algo import grnboost2
from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='Run GRNBoost2 for gene regulatory network inference')
    parser.add_argument('--input', required=True, help='Path to input expression matrix')
    parser.add_argument('--tf-list', required=True, help='Path to TF list file')
    parser.add_argument('--output', required=True, help='Path to output network file')
    parser.add_argument('--cores', type=int, default=1, help='Number of cores to use')
    args = parser.parse_args()
    
    print(f"Loading expression matrix from {args.input}")
    expr_matrix = pd.read_csv(args.input, index_col=0)
    
    # Transpose if needed (samples should be rows, genes columns)
    if expr_matrix.shape[0] < expr_matrix.shape[1]:
        print("Transposing expression matrix...")
        expr_matrix = expr_matrix.T
        
    print(f"Expression matrix shape: {expr_matrix.shape} (samples × genes)")
    
    # Load TFs
    print(f"Loading TF list from {args.tf_list}")
    with open(args.tf_list) as f:
        tf_names = [line.strip() for line in f]
    print(f"Loaded {len(tf_names)} transcription factors")
    
    # Set up parallel processing
    print(f"Setting up Dask cluster with {args.cores} workers")
    cluster = LocalCluster(n_workers=args.cores, threads_per_worker=1)
    client = Client(cluster)
    print(f"Dashboard link: {client.dashboard_link}")
    
    # Run GRNBoost2
    print("Running GRNBoost2 network inference...")
    with ProgressBar():
        network = grnboost2(expression_data=expr_matrix, 
                           tf_names=tf_names,
                           verbose=True)
    
    # Save the network
    print(f"Saving network to {args.output}")
    network.to_csv(args.output, sep='\t', index=False)
    print("GRNBoost2 network inference completed.")
    
    # Close the client
    client.close()
    cluster.close()
    
if __name__ == "__main__":
    main()
