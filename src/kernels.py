import numpy as np

def compute_identity_matrix(sequences):
    seq_array = np.array([list(seq) for seq in sequences])
    N = len(sequences)
    mat = np.zeros((N, N))
    
    for i in range(N):
        mask_not_gap = ~((seq_array[i] == '-') & (seq_array[i:] == '-'))
        matches = ((seq_array[i] == seq_array[i:]) & mask_not_gap).sum(axis=1)
        lengths = mask_not_gap.sum(axis=1)
        
        safe_lengths = np.where(lengths > 0, lengths, 1)
        scores = matches / safe_lengths
        
        mat[i, i:] = scores
        mat[i:, i] = scores
        
    return mat

def create_similarity_matrix(sequences, kernel="id", gamma=1.0):
    base_mat = compute_identity_matrix(sequences)
    
    if kernel == "id":
        return base_mat
    elif kernel == "exp":
        return np.exp(-gamma * (1.0 - base_mat))
    elif kernel == "poly":
        return base_mat ** gamma
    else:
        raise ValueError(f"Unknown kernel '{kernel}'. Choose 'id', 'exp' or 'poly'.")