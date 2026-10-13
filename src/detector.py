import numpy as np
from src.kernels import create_similarity_matrix 

class OutlierClassifier:
    def __init__(self, sequences, seq_names, kernel="id", gamma=1.0):
        """
        Initializes the classifier directly with sequences, avoiding slow CSV reads.
        """
        self._n = len(sequences)
        self._seq_names = np.array(seq_names)
        
        # Uses the optimized vectorized NumPy matrix generator
        self._K = create_similarity_matrix(sequences, kernel=kernel, gamma=gamma)

        # Verifies that the similarity matrix is semi-positive definite
        min_eig = np.min(np.linalg.eigvals(self._K))
        if min_eig < -1e-10:
            raise ValueError("The similarity matrix is not semi-positive definite.")

    def fit_predict(self, delta=0.05, verbose=False):
        results = []
        for i in range(self._n):
            is_outlier, f_val, thresh = self.fit_predict_one(i, delta=delta, verbose=verbose)
            results.append({
                "protein_name": self._seq_names[i],
                "is_outlier": is_outlier,
                "score": round(f_val - thresh, 4),
                "f_outlier": round(f_val, 4),
                "threshold": round(thresh, 4),
            })
            
        # Optional: You can still return a Pandas DataFrame at the very end
        # because returning a neat table of results is what users expect.
        import pandas as pd
        return pd.DataFrame(results).sort_values(by="score", ascending=False)
    
    def fit_predict_one(self, index_outlier, delta=0.05, verbose=False):
        # Remove the outlier from the dataset and create the similarity matrix
        K_train = np.delete(
            np.delete(self._K, index_outlier, axis=0), index_outlier, axis=1
        )
        m = self._n - 1
        
        # Compute the mean of the distance matrix (without outlier)
        moy = np.mean(K_train)

        # Compute the values of f for the sequences in the training set
        # Vectorized version for speed:
        f_train = np.diag(K_train) - 2 * (1/m) * np.sum(K_train, axis=1) + moy

        # Compute the threshold for outlier detection
        c = 1.0 / np.sqrt(m) 
        rad = 2.0 * c
        threshold = np.mean(f_train) + rad + np.sqrt(np.log(1 / delta) / (2 * m))

        # Identify outliers
        f_outlier = 1 - 2 * (1/m) * (np.sum(self._K[index_outlier, :]) - 1) + moy

        if verbose:
            print(f"Sequence: {self._seq_names[index_outlier]}")
            print(f"f_outlier: {f_outlier:.4f} | threshold: {threshold:.4f}")
            print("OUTLIER" if f_outlier > threshold else "NORMAL")
        
        return f_outlier > threshold, f_outlier, threshold