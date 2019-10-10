import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3
sc.logging.print_versions()

results_file = '../data/processed/test.csv'

adata = sc.read_csv('data/raw/sc_10x_5cl.count.csv', first_column_names=True)


print(adata)

sc.pl.highest_expr_genes(adata, n_top=20)