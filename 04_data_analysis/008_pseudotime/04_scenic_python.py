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
    parser.add_argument('--input', required=True, help='Path to input expression matrix (cells as rows, genes as columns)')
    parser.add_argument('--tf-list', required=True, help='Path to TF list file')
    parser.add_argument('--output', required=True, help='Path to output network file')
    parser.add_argument('--cores', type=int, default=1, help='Number of cores to use')
    args = parser.parse_args()
    
    print(f"\n📖 Loading expression matrix from: {args.input}")
    expr_matrix = pd.read_csv(args.input, index_col=0)
    print(f"✅ Expression matrix loaded with shape: {expr_matrix.shape} (samples × genes)")
    print(f"📝 First 5 gene names (columns): {expr_matrix.columns[:5].tolist()}\n")

    # Load TFs
    print(f"📖 Loading TF list from: {args.tf_list}")
    with open(args.tf_list) as f:
        tf_names = [line.strip() for line in f]
    print(f"✅ Loaded {len(tf_names)} transcription factors")

    # Check how many TFs are actually in the expression matrix
    matched_tfs = [tf for tf in tf_names if tf in expr_matrix.columns]
    print(f"🔍 Found {len(matched_tfs)} TFs present in the expression matrix out of {len(tf_names)}")
    if len(matched_tfs) == 0:
        print("⚠️  No matching TFs found in the expression matrix — check your inputs.")
        return

    # Set up parallel processing
    print(f"\n⚙️  Setting up Dask cluster with {args.cores} workers")
    cluster = LocalCluster(n_workers=args.cores, threads_per_worker=1)
    client = Client(cluster)
    print(f"✅ Dask dashboard link: {client.dashboard_link}\n")

    # Run GRNBoost2
    print("🚀 Running GRNBoost2 network inference...")
    with ProgressBar():
        network = grnboost2(expression_data=expr_matrix, 
                            tf_names=matched_tfs,
                            verbose=True)

    print(f"✅ Inferred network with {network.shape[0]} edges")
    
    # Save the network
    print(f"\n💾 Saving network to: {args.output}")
    network.to_csv(args.output, sep='\t', index=False)
    print("🎉 GRNBoost2 network inference completed successfully.")

    # Close the client
    client.close()
    cluster.close()
    print("🛑 Dask cluster shut down.")

if __name__ == "__main__":
    main()
