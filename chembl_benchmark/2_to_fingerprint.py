import time
import pandas as pd
from rdkit.Chem import AllChem
from mhfp.encoder import MHFPEncoder

# Let's get sparkly
from pyspark import SparkContext, SparkConf
sparkConf = SparkConf().setAppName('ChEMBL Benchmark: Fingerprint Calculation')
sc = SparkContext(conf=sparkConf)


def convert(subset):
    target = '/cluster/chembl/chembl.' + str(subset) + '.smi'
    actives = pd.read_csv(target, sep=' ', usecols=[0], header=None)
    
    mh = MHFPEncoder()

    with open('/cluster/chembl/chembl.' + str(subset) + '.mhfp6', 'w+') as f:
        for _, row in actives.iterrows():
            mol = AllChem.MolFromSmiles(row[0])
            if mol:
                fp_vals = ','.join(map(str, mh.encode_mol(mol)))

                f.write(fp_vals + '\n')

    with open('/cluster/chembl/chembl.' + str(subset) + '.mhecfp4', 'w+') as f:
        for _, row in actives.iterrows():
            mol = AllChem.MolFromSmiles(row[0])
            if mol:
                fp_vals = ','.join(map(str, mh.from_sparse_array([*AllChem.GetMorganFingerprint(mol, 2).GetNonzeroElements()])))

                f.write(fp_vals + '\n')

    with open('/cluster/chembl/chembl.' + str(subset) + '.ecfp4', 'w+') as f:
        for _, row in actives.iterrows():
            mol = AllChem.MolFromSmiles(row[0])
            if mol:
                fp_vals = ','.join(map(str, AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)))

                f.write(fp_vals + '\n')
                
time.sleep(60)
sc.parallelize(range(20)).foreach(convert)