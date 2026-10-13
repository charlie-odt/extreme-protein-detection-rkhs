import os

def read_fasta(file_path):
    """
    Reads a FASTA file and extracts sequence names and strings.
    Handles multi-line sequences correctly.
    
    Args:
        file_path (str): Path to the .fasta or .aln file.
        
    Returns:
        tuple: (list of names, list of sequences)
    """
    names = []
    sequences = []
    
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            current_name = ""
            current_seq = []
            
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith(">"):
                    if current_name:
                        names.append(current_name)
                        sequences.append("".join(current_seq))
                    
                    current_name = line[1:] 
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Don't forget the last sequence in the file!
            if current_name:
                names.append(current_name)
                sequences.append("".join(current_seq))
                
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: The file {file_path} was not found.")
        
    # Validation step: in a multiple sequence alignment, 
    # all sequences must have the exact same length.
    if sequences:
        first_len = len(sequences[0])
        if any(len(s) != first_len for s in sequences):
            raise ValueError(f"Invalid alignment in {file_path}: Sequences have different lengths.")
            
    return names, sequences

def write_fasta(file_path, good_seq_names, all_sequences, all_names=None):
    """
    Writes a new FASTA file containing only the non-outlier sequences.
    
    Args:
        file_path (str): Output path for the new FASTA file.
        good_seq_names (list): List of sequence names to keep (the 'normals').
        all_sequences (list): List of all original sequences.
        all_names (list, optional): List of all original sequence names.
                                    If None, assumes good_seq_names and all_sequences match.
    """
    # Map names to their sequences for safe extraction
    if all_names:
        seq_dict = dict(zip(all_names, all_sequences))
        sequences_to_write = [(name, seq_dict[name]) for name in good_seq_names if name in seq_dict]
    else:
        # Fallback if just zipped directly (less safe if order changed during processing)
        sequences_to_write = zip(good_seq_names, all_sequences)

    with open(file_path, "w", encoding="utf-8") as f:
        for name, seq in sequences_to_write:
            # Standard FASTA format: sequence name starts with '>'
            f.write(f">{name}\n")
            
            # Break sequence into chunks of 80 characters for standard FASTA compliance
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")