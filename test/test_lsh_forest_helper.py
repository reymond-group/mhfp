import pytest
import numpy as np
from scipy.spatial.distance import jaccard
from rdkit.Chem import AllChem
from mhfp.encoder import MHFPEncoder
from mhfp.lsh_forest import LSHForestHelper

# Keeping tests barebone and simple

mhfp_encoder = MHFPEncoder()
lfh = LSHForestHelper()

drugbank = []

with open('test/drugbank.smi') as f:
    for line in f.readlines():
        mol = AllChem.MolFromSmiles(line.strip().split()[0])
        if mol:
            drugbank.append(mhfp_encoder.encode_mol(mol))

for i, fp in enumerate(drugbank):
    lfh.add(i, fp)

lfh.index()

def test_setup():
    assert len(drugbank) == 226

def test_add():
    assert len(lfh.lsh_forest.keys) == 226

def test_query():
    assert lfh.query(drugbank[0], 6, drugbank) == [0, 3, 4, 7, 5, 9]
    assert lfh.query(drugbank[0], 6, drugbank, 100) == [0, 3, 4, 7, 5, 9]
    print(drugbank[0])



