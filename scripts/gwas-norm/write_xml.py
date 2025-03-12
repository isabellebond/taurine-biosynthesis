#!/usr/bin/env python
"""
XML writer script for gwas norm
"""

from gwas_norm import common
from gwas_norm.elements import study_obj, gwas_data_obj, file_obj, \
    analysis_obj, test_obj, phenotype_obj as ph, cohort_obj
import urllib.parse
import urllib.request
import urllib.error
import os
import csv
import gzip
import sys
import argparse
import glob
import decimal
from pathlib import Path
import pandas as pd
import shutil

# The column mappings for a data file
INDIR = "/lustre/projects/mol_cardio/resources/gwas-format/data/raw/surendran_metabolome_36357675/formatted"
STUDY_NAME = "surendran_metabolome"
EFFECT_TYPE = 'beta' #results are output from bolt-llm
COHORT = 'Multiple'
target_genome_assemblies = ['b38']
# CONSORTIUM = "DECODE"
COMPRESSION = 'gzip'
SOURCE_ASSEMBLY = 'b38'
ANALYSIS_TYPE = 'trait'
PUBMED_ID = '36357675'
INFO = ""
PVALUE_LOGGED = False
URL = 'https://www.nature.com/articles/s41591-022-02046-0'
DELIMITER = "\t"


colmap = dict(
    var_id='MarkerName',
    chr_name='chromosome',
    start_pos='start_pos',
    end_pos='end_pos',
    pvalue='P-value',
    effect_size='Effect',
    standard_error='StdErr',
    effect_allele="effect_allele",
    other_allele="other_allele",
    number_of_samples="N",
    effect_allele_freq = 'Freq1')

# The dest dir is completely relative to the ENV dest dir
DEST_DIR = "/home/rmgpibo/taurine-biosynthesis/data/gwas-norm/{0}_{1}".format(STUDY_NAME, PUBMED_ID)

infiles = os.listdir(INDIR)

# The central GwasData object
gd = gwas_data_obj.GwasData()

# data has files at analysis level (single file per analysis), not the study level
st = study_obj.Study(
    STUDY_NAME,
    INDIR,
    DEST_DIR,
    SOURCE_ASSEMBLY,
    pubmed_id=PUBMED_ID,
    target_genome_assemblies=target_genome_assemblies,
    url=URL
)


for i in infiles:
    trait = i.split('_')[-1]
    trait = trait.split('.')[0]    # Create analysis (analysis name is SeqID)

    a = analysis_obj.AnalysisFile(
        trait, ANALYSIS_TYPE, EFFECT_TYPE
        )

    # create phenotypes
    pheno_def = ph.Definition(trait, type = "text")
    a.phenotype = ph.Phenotype(pheno_def)

    # add files to the analysis
    file_element = file_obj.GwasFile(
        i, colmap,
        pvalue_logged=PVALUE_LOGGED,
        compression=COMPRESSION, delimiter=DELIMITER,
        file_check=True
        )

    a.add_file(file_element)
    st.add_analysis(a)

gd.add_study(st)

# write to xml
gd.write(f'data/gwas-norm/{STUDY_NAME}.xml')