# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
# %run setup.ipynb

# %%
chrom = '2L'

# %%
callset_pass = zarr.open_group('../phase2.AR1/variation/main/zarr2/ag1000g.phase2.ar1.pass.biallelic', mode='a')
callset_pass

# %%
#snp positions
pos = allel.SortedIndex(callset_pass_biallelic[chrom]['variants/POS'])
pos

# %%
