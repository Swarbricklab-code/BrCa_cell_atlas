# -*- coding: utf-8 -*-
"""
Convert sparse matrix to dense 

AUTHOR: Daniel Roden
CREATED: 15/05/2020
MODIFIED: 14/07/2020

NOTES: depends on snakemake input and output file info
"""
# =============================
print("CONFIG", flush=True)

from scipy import sparse, io
import pandas as pd
import numpy as np
from datetime import datetime 

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("----------------------------------", flush=True)
print("Start Time =", current_time, flush=True)
print("----------------------------------", flush=True)

# ----------------------------------------
# PARAMS

# INPUT
file_sparse_matrix = snakemake.input["sparse_matrix"]
print("file_sparse_matrix: " + file_sparse_matrix, flush=True)

file_sparse_matrix_columns = snakemake.input["cells"]
print("file_sparse_matrix_columns: " + file_sparse_matrix_columns, flush=True)

file_sparse_matrix_rows = snakemake.input["genes"]
print("file_sparse_matrix_rows: " + file_sparse_matrix_rows, flush=True)

# OUTPUT
file_dense_matrix = snakemake.output["dense_matrix"]
print("file_dense_matrix: " + file_dense_matrix, flush=True)

# WILDCARDS


# --------------------------------------------------------
print("READING sparse matrix")
counts_sparsematrix = io.mmread(file_sparse_matrix)
counts_densematrix = counts_sparsematrix.toarray()

row_names = np.genfromtxt(file_sparse_matrix_rows, dtype=str)
col_names = np.genfromtxt(file_sparse_matrix_columns, dtype=str)

print("WRITING dense matrix")
df = pd.DataFrame(counts_densematrix, columns=col_names, index=row_names)
df.to_csv(file_dense_matrix, sep='\t', header=True, index=True, index_label='Somelabel')

# --------------------------------------------------------
print("COMPLETED", flush=True)
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("----------------------------------", flush=True)
print("End Time =", current_time, flush=True)
print("----------------------------------", flush=True)
