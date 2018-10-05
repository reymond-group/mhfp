import pytest
import numpy as np
from scipy.spatial.distance import jaccard
from rdkit.Chem import AllChem
from mhfp.encoder import MHFPEncoder
from mhfp.lsh_forest import LSHForest

# Keeping tests barebone and simple

mhfp_encoder = MHFPEncoder()
lf = LSHForest()

drugbank = []

with open('test/drugbank.smi') as f:
    for line in f.readlines():
        mol = AllChem.MolFromSmiles(line.strip().split()[0])
        if mol:
            drugbank.append(mhfp_encoder.encode_mol(mol))

for i, fp in enumerate(drugbank):
    lf.add(i, fp)

lf.index()

def test_setup():
    assert len(drugbank) == 226

def test_add():
    assert len(lf.keys) == 226

def test_query():
    assert lf.query(drugbank[0], 6) == [0, 3, 4, 197, 9, 149]



