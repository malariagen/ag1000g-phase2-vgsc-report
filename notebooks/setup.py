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

# %% [markdown]
# # Setup for ag1000g phase 2 analysis

# %%
# %%HTML
<style type="text/css">
.container {
    width: 100%;
}
</style>

# %%
# python standard library
import sys
import os
import operator
import itertools
import collections
import functools
import glob
import csv
import datetime
import bisect
import sqlite3
import subprocess
import random
import gc
import shutil
import shelve
import contextlib
import tempfile
import math
import warnings

# %%
# plotting setup
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
import seaborn as sns
sns.set_context('paper')
sns.set_style('ticks')
# use seaborn defaults
rcParams = plt.rcParams
rcParams['savefig.jpeg_quality'] = 100

# %%
# %matplotlib inline
# %config InlineBackend.figure_formats = {'retina', 'png'}

# %%
# general purpose third party packages
import numpy as np
nnz = np.count_nonzero
import scipy
import scipy.stats
import scipy.spatial.distance
import numexpr
import h5py
import tables
import bcolz
import dask
import dask.array as da
import pandas as pd
import IPython
from IPython.display import clear_output, display, HTML
import sklearn
import sklearn.decomposition
import sklearn.manifold
import petl as etl
etl.config.display_index_header = True
import humanize
from humanize import naturalsize, intcomma, intword
import zarr
from scipy.stats import entropy
import lmfit

# %%
#analysis packages
import allel

# %%
sys.path.insert(0, '../agam-report-base/src/python')
from util import *

# %%
from ag1k import phase2_ar1

# %%
# This is a symlink in your root directory
# eg: ln -s /kwiat/vector/ag1000g/release/phase2.AR1 .
phase2_ar1.init("../phase2.AR1")

# %%
region_vgsc = SeqFeature('2L', 2358158, 2431617, label='Vgsc')

# %%
import veff

# %%
