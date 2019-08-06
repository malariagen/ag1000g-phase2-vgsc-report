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
# # Extract mutations in VGSC
#
# This notebook extracts data on all mutations in the VGSC gene.

# %% [markdown]
# ## Setup

# %%
# %run setup.ipynb

# %%
# download gene annotations from vectorbase
# !wget \
#     --no-clobber \
#     -O ../data/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.4.gff3.gz \
#     https://www.vectorbase.org/download/anopheles-gambiae-pestbasefeaturesagamp44gff3gz


# %%
# download the Davies et al. (2007) gene models
# !wget \
#     --no-clobber \
#     -O ../data/davies_vgsc_model_20170125.gff3 \
#     http://alimanfoo.github.io/assets/davies_vgsc_model_20170125.gff3


# %%
# load the vectorbase geneset
geneset_agamp44 = allel.FeatureTable.from_gff3('../data/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.4.gff3.gz',
                                               attributes=['ID', 'Parent'])
geneset_agamp44 = geneset_to_pandas(geneset_agamp44)
geneset_agamp44.head()

# %%
# subset to VGSC
geneset_agamp44_vgsc = geneset_agamp44.query(region_vgsc.query_str).copy()
# replace CDS IDs as not informative
geneset_agamp44_vgsc['ID'].values[(geneset_agamp44_vgsc.type == 'CDS').values] = ''
geneset_agamp44_vgsc.type.value_counts()

# %%
# load the Davies geneset
geneset_davies = allel.FeatureTable.from_gff3('../data/davies_vgsc_model_20170125.gff3',
                                              attributes=['ID', 'Parent'])
geneset_davies = geneset_to_pandas(geneset_davies)
geneset_davies.head()

# %%
# make a combined geneset
geneset_vgsc_combined = pandas.concat([geneset_agamp44_vgsc, geneset_davies])
geneset_vgsc_combined.query("type == 'mRNA'")

# %%
# setup a variant annotator
annotator = veff.Annotator(
    fasta_path='../phase2.AR1/genome/agamP3/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa', 
    gff3_path=['../data/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.4.gff3.gz',
               '../data/davies_vgsc_model_20170125.gff3'],
    seqid='2L'
)

# %%
# identify VGSC transcripts
transcript_ids = [f.feature_id for f in annotator.get_children('AGAP004707')]
transcript_ids

# %%
# tabulate Davies exons
tbl_davies_exons = (
    etl
    .fromdataframe(geneset_davies)
    .eq('type', 'CDS')
    .cutout('Parent', 'source', 'type', 'score', 'strand', 'phase')
    .merge(key=('start', 'end'))
    .rename('seqid', 'exon_seqid')
    .rename('ID', 'exon')
    .rename('start', 'exon_start')
    .rename('end', 'exon_end')
    .movefield('exon_seqid', 0)
)
tbl_davies_exons.displayall()

# %% [markdown]
# ## Extract table of variants

# %%
callset = phase2_ar1.callset
callset

# %%
# what fields are available?
print(', '.join(callset['2L/variants']))

# %%
#get SNPEFF annotations from HDF5 - also find out why these aren't in zarr format
snpeff_h5_fn = '../phase2.AR1/variation/main/hdf5/all_snpeff/ag1000g.phase2.ar1.snpeff.AgamP4.2.2L.h5'
snpf = h5py.File(snpeff_h5_fn, mode='r')

# %%
# what SNPEFF fields are available?
print(', '.join(snpf['2L/variants/ANN'].dtype.names))

# %%
samples = phase2_ar1.df_samples
samples.head()

# %%
#samples needs to have a numeric index - this should probably be tested in future versions
samples = samples.reset_index()


# %% [markdown]
# ### breakdown table code - not working

# %%
# variants = callset[seqid]['variants']
# ann = snpf[seqid]['variants']['ANN']
# pos = allel.SortedIndex(variants['POS'])

# %%
# start = region_vgsc.start
# end = region_vgsc.end

# %%
# loc = pos.locate_range(start, end)
# genotype = allel.GenotypeArray(callset[seqid]['calldata/genotype'][loc])

# %%
# acs = genotype.count_alleles_subpops(max_allele=3, subpops=subpops)

# %%
def tabulate_variants(callset, snpeff, seqid, start, end, pop_ids, subpops):
    """Build a table of variants for a given callset and genome region."""
    
    variants = callset[seqid]['variants']
    ann = snpeff[seqid]['variants']['ANN']
    pos = allel.SortedIndex(variants['POS'])
    loc = pos.locate_range(start, end)
    genotype = allel.GenotypeArray(callset[seqid]['calldata/genotype'][loc])
    acs = genotype.count_alleles_subpops(max_allele=3, subpops=subpops)
    
    # extract columns
    variants_fields = [
        'CHROM',
        'POS',
        'num_alleles',
        'REF',
        'ALT',
        'AC',
        'FILTER_PASS',
        'NoCoverage',
        'LowCoverage',
        'HighCoverage',
        'LowMQ',
        'HighMQ0',
        'RepeatDUST',
        'RepeatMasker',
        'RepeatTRF',
        'FS',
        'HRun',
        'QD',
        'ReadPosRankSum',
    ]
    ann_fields = ['Allele', 'Annotation', 'HGVS_c', 'HGVS_p', 'Feature_ID', 'CDS_pos']
    cols = (
        [variants[f][loc] for f in variants_fields] + 
        [ann[loc][f] for f in ann_fields] + 
        [acs[p].to_frequencies() for p in pop_ids]
    )

    def split_alleles(row):
        for i in range(row.num_alleles - 1):
            # break down alleles
            out = [
                row['CHROM'], 
                row['POS'], 
                row['num_alleles'], 
                row['REF'], 
                row['ALT'][i], 
                row['AC'][i], 
                i, 
            ]
            # add in remaining variant annotations
            out += [row[f] for f in variants_fields[6:]]
            # SNPEFF annotation only applies to first allele
            if i == 0:
                out += [row[f] for f in ann_fields]
            else:
                out += [None for f in ann_fields]
            # add in population allele frequencies
            out += [row[p][i+1] for p in pop_ids]
            yield out
        
    tbl = (
        etl
        .fromcolumns(cols, header=variants_fields + ann_fields + list(pop_ids))
        .rowmapmany(split_alleles, header=variants_fields[:6] + ['ALTIX'] + variants_fields[6:] + ann_fields + list(pop_ids), failonerror=True)
        .convert('CHROM REF ALT Allele Annotation HGVS_c HGVS_p Feature_ID'.split(), lambda v: str(v, 'ascii'))
        .rename({f: 'SNPEFF_' + f for f in ann_fields})
        .rename({p: 'AF_%s' % p for p in pop_ids})
        .addfield('check_allele', lambda row: row['SNPEFF_Allele'] is None or row['SNPEFF_Allele'] == row['ALT'])
    )
    
    return tbl


# %%
pop_ids = phase2_ar1.pop_ids
print(', '.join(pop_ids))

# %%
subpops = {p: samples[samples.population == p].index.values.tolist() for p in pop_ids}

# %%
# build a table of variants from phase 1
tbl_variants_phase2 = tabulate_variants(callset, snpf, 
                                        seqid=region_vgsc.seqid, start=region_vgsc.start, end=region_vgsc.end, 
                                        pop_ids=pop_ids, subpops=subpops)
tbl_variants_phase2

# %%
#let's have a look shall we
tbl_variants_phase2.tocsv('phase2_variants.csv')

# %% [markdown]
# ## Annotate effects for all transcripts

# %%
cds_effects = [
    'NON_SYNONYMOUS_CODING', 
    'SYNONYMOUS_CODING',    
]
intron_effects = [
    'INTRONIC', 
    'SPLICE_CORE',
    'SPLICE_REGION',        
]
selected_effects = cds_effects + intron_effects


# %%
def lpop(l, default=None):
    """Pop the first item from a list if not empty."""
    try:
        return l[0]
    except IndexError:
        return default



# %%
def transcript_effect(transcript_id):
    def f(row):
        e = lpop([e for e in row.VEFF if e.transcript_id == transcript_id])
        if e and e.effect in cds_effects:
            return (e.effect, e.aa_change)
        elif e and e.effect in intron_effects:
            return (e.effect, e.intron_cds_5prime, e.intron_5prime_dist, e.intron_cds_3prime, e.intron_3prime_dist)
        else:
            return None
    return f



# %%
tbl_variants_phase2_eff = (
    tbl_variants_phase2
    # join in Davies exon information
    .intervalleftjoin(
        # don't include shorter exon alternatives
        tbl_davies_exons.select('exon', lambda v: v[-1] != '-'),
        lkey='CHROM', rkey='exon_seqid', lstart='POS', rstart='exon_start', lstop='POS', rstop='exon_end', include_stop=True)
    .cutout('exon_seqid')
    .addfield('VEFF', lambda row: [e for e in annotator.get_effects(chrom=row.CHROM, pos=row.POS, ref=row.REF, alt=row.ALT) 
                                   if e.effect in selected_effects])
    .addfield(transcript_ids[0], transcript_effect(transcript_ids[0]))
    .addfield(transcript_ids[1], transcript_effect(transcript_ids[1]))
    .addfield(transcript_ids[2], transcript_effect(transcript_ids[2]))
    .addfield(transcript_ids[3], transcript_effect(transcript_ids[3]))
    .addfield(transcript_ids[4], transcript_effect(transcript_ids[4]))
    .addfield(transcript_ids[5], transcript_effect(transcript_ids[5]))
    .addfield(transcript_ids[6], transcript_effect(transcript_ids[6]))
    .addfield(transcript_ids[7], transcript_effect(transcript_ids[7]))
    .addfield(transcript_ids[8], transcript_effect(transcript_ids[8]))
    .addfield(transcript_ids[9], transcript_effect(transcript_ids[9]))
    .addfield(transcript_ids[10], transcript_effect(transcript_ids[10]))
    .addfield(transcript_ids[11], transcript_effect(transcript_ids[11]))
    .addfield(transcript_ids[12], transcript_effect(transcript_ids[12]))
    .cutout('VEFF')
    .replaceall('.', None)
    .replaceall('', None)
    .cache()
)

# %%
tbl_variants_phase2_eff


# %% [markdown]
# ## Inspect missense variants

# %%
def simplify_missense_effect(v):
    if v and v[0] == 'NON_SYNONYMOUS_CODING':
        return v[1]
    else:
        return ''

    
td_styles = {
    'FILTER_PASS': lambda v: 'background-color: red' if not v else '',
    'NoCoverage': lambda v: 'background-color: red' if v > 1 else '',
    'LowCoverage': lambda v: 'background-color: red' if v > 76 else '',
    'HighCoverage': lambda v: 'background-color: red' if v > 15 else '',
    'LowMQ': lambda v: 'background-color: red' if v > 76 else '',
    'HighMQ0': lambda v: 'background-color: red' if v > 1 else '',
    'RepeatDUST': lambda v: 'background-color: red' if v else '',
    'FS': lambda v: 'background-color: red' if v > 60 else '',
    'QD': lambda v: 'background-color: red' if v < 5 else '',
    'ReadPosRankSum': lambda v: 'background-color: red' if v < -8 else '',
    'HRun': lambda v: 'background-color: red' if v > 4 else '',
    'num_alleles': lambda v: 'background-color: orange' if v > 2 else '',
}


def tr_style(row):
    """Colour row by alternate allele count."""
    return 'background-color:rgba(0, 255, 0, %.3f)' % (min(1, row['AC']/100))


tbl_variants_phase1_missense = (
    tbl_variants_phase1_eff
    .select(lambda row: any(row[t] and row[t][0] == 'NON_SYNONYMOUS_CODING' for t in transcript_ids))
    .convert(transcript_ids, simplify_missense_effect)
)
tbl_variants_phase1_missense.displayall(td_styles=td_styles, tr_style=tr_style)


# %% [markdown]
# ## Inspect splice site variants

# %%
def simplify_intron_effect(v):
    if v and v[0] in ['SPLICE_REGION', 'SPLICE_CORE']:
        if math.fabs(v[2]) < math.fabs(v[4]):
            return v[1], v[2]
        else:
            return v[3], v[4]
    else:
        return ''

    
td_styles = {
    'FILTER_PASS': lambda v: 'background-color: red' if not v else '',
    'NoCoverage': lambda v: 'background-color: red' if v > 1 else '',
    'LowCoverage': lambda v: 'background-color: red' if v > 76 else '',
    'HighCoverage': lambda v: 'background-color: red' if v > 15 else '',
    'LowMQ': lambda v: 'background-color: red' if v > 76 else '',
    'HighMQ0': lambda v: 'background-color: red' if v > 1 else '',
    'RepeatDUST': lambda v: 'background-color: red' if v else '',
    'FS': lambda v: 'background-color: red' if v > 60 else '',
    'QD': lambda v: 'background-color: red' if v < 5 else '',
    'ReadPosRankSum': lambda v: 'background-color: red' if v < -8 else '',
    'HRun': lambda v: 'background-color: red' if v > 4 else '',
    'num_alleles': lambda v: 'background-color: orange' if v > 2 else '',
}


def tr_style(row):
    """Colour row by alternate allele count."""
    return 'background-color:rgba(0, 255, 0, %.3f)' % (min(1, row['AC']/100))


tbl_variants_phase1_splice = (
    tbl_variants_phase1_eff
    .select(lambda row: any(row[t] and row[t][0] in ['SPLICE_REGION', 'SPLICE_CORE'] for t in transcript_ids))
    .convert(transcript_ids, simplify_intron_effect)
)
tbl_variants_phase1_splice.displayall(td_styles=td_styles, tr_style=tr_style)

# %% [markdown]
# ## Write out variants to file

# %%
(tbl_variants_phase1_eff
 .teepickle('../data/tbl_variants_phase1.pkl')
 .convert(transcript_ids, lambda v: ':'.join(map(str, v)))
 .replaceall(None, 'NA')
 .tocsv('../data/tbl_variants_phase1.csv')
)

# %%
# check OK
etl.frompickle('../data/tbl_variants_phase1.pkl')

# %%
etl.fromcsv('../data/tbl_variants_phase1.csv')

# %%
