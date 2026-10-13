import argparse
import os
import pandas as pd
from src.parser import read_fasta, write_fasta
from src.detector import OutlierClassifier

def main():
    # 1. Setup the command line argument parser
    parser = argparse.ArgumentParser(
        description="Filters divergent homologous sequences (outliers) in a multiple FASTA alignment using RKHS theory."
    )
    
    # Required arguments
    parser.add_argument("-i", "--input", required=True, help="Path to the input FASTA alignment file.")
    parser.add_argument("-o", "--output", required=True, help="Path for the output CSV/FASTA file or directory.")
    
    # Algorithm parameters (with default values)
    parser.add_argument("--kernel", type=str, default="id", choices=["id", "exp", "poly"], help="Type of kernel to use (id, exp, poly).")
    parser.add_argument("--gamma", type=float, default=1.0, help="Gamma parameter for exp and poly kernels (default: 1.0).")
    parser.add_argument("--delta", type=float, default=0.05, help="Delta risk parameter for the RKHS threshold (default: 0.05).")
    
    # Additional options
    parser.add_argument("--clean-fasta", action="store_true", help="If set, generates a new FASTA file without outliers in addition to the CSV.")
    parser.add_argument("--verbose", action="store_true", help="Prints detection details for each sequence.")

    # Parse arguments provided by the user in the terminal
    args = parser.parse_args()

    # 2. Execute the pipeline with provided arguments
    print(f"Loading alignment: {args.input}")
    if not os.path.exists(args.input):
        print(f"Error: The file {args.input} does not exist.")
        return

    # A. Parse the FASTA file
    seq_names, sequences = read_fasta(args.input) 
    print(f"{len(sequences)} sequences loaded. Computing the '{args.kernel}' kernel matrix...")

    # B. Initialize the classifier
    try:
        classifier = OutlierClassifier(sequences, seq_names, kernel=args.kernel, gamma=args.gamma)
    except ValueError as e:
        print(f"Mathematical Error: {e}")
        return

    # C. Run outlier detection
    print(f"Detecting outliers (delta={args.delta})...")
    df_results = classifier.fit_predict(delta=args.delta, verbose=args.verbose)
    
    # D. Save the results
    # Create the output directory if it does not exist
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    
    df_results.to_csv(args.output, index=False)
    total_outliers = df_results["is_outlier"].sum()
    print(f"\nDone! {total_outliers} outliers detected.")
    print(f"Detailed results saved to: {args.output}")
    
    # E. (Optional) Generate the clean FASTA
    if args.clean_fasta:
        # Extract names of sequences that are NOT outliers
        good_seq_names = df_results[~df_results["is_outlier"]]["protein_name"].tolist()
        
        # Generate the new filename and write the clean FASTA
        clean_fasta_path = args.output.replace(".csv", "_clean.fasta")
        write_fasta(clean_fasta_path, good_seq_names, sequences)
        print(f"Cleaned alignment saved to: {clean_fasta_path}")

if __name__ == "__main__":
    main()