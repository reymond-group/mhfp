import time
import pandas as pd
import numpy as np
from operator import itemgetter
from rdkit.Chem import AllChem
from scipy.spatial.distance import jaccard

# Let's get sparkly
from pyspark import SparkContext, SparkConf
sparkConf = SparkConf().setAppName('ChEMBL Benchmark: Create Sorted Lists (Linear Scan)')
sc = SparkContext(conf=sparkConf)

def process(subset):
    mhfps = []
    ecfps = []
    mhecfps = []

    np.random.seed(42)
    samples = np.random.choice(10000, 20, replace=False)

    with open('/cluster/chembl/chembl.' + str(subset) + '.mhfp6') as f:
        for line in f.readlines():
            mhfps.append(np.array(list(map(int, line.split(','))), dtype=np.uint32))

    with open('/cluster/chembl/chembl.' + str(subset) + '.ecfp4') as f:
        for line in f.readlines():
            ecfps.append(np.array(list(map(int, line.split(','))), dtype=np.uint8))

    with open('/cluster/chembl/chembl.' + str(subset) + '.mhecfp4') as f:
        for line in f.readlines():
            mhecfps.append(np.array(list(map(int, line.split(','))), dtype=np.uint32))

    for sample in samples:
        query = mhecfps[sample]
        dists = []

        for i, mhecfp in enumerate(mhecfps):
            dists.append((i, jaccard(query, mhecfp)))
        
        dists = sorted(dists, key=itemgetter(1))

        with open('/cluster/chembl/dists_mhecfp_' + str(subset) + '_' + str(sample) + '.txt', 'w+') as f:
            for dist in dists:
                f.write(str(dist[0]) + ',' + str(dist[1]) + '\n')

        query = mhfps[sample]
        dists = []

        for i, mhfp in enumerate(mhfps):
            dists.append((i, jaccard(query, mhfp)))
        
        dists = sorted(dists, key=itemgetter(1))

        with open('/cluster/chembl/dists_mhfp_' + str(subset) + '_' + str(sample) + '.txt', 'w+') as f:
            for dist in dists:
                f.write(str(dist[0]) + ',' + str(dist[1]) + '\n')


        query = ecfps[sample]
        dists = []

        for i, ecfp in enumerate(ecfps):
            dists.append((i, jaccard(query, ecfp)))
        
        dists = sorted(dists, key=itemgetter(1))

        with open('/cluster/chembl/dists_ecfp_' + str(subset) + '_' + str(sample) + '.txt', 'w+') as f:
            for dist in dists:
                f.write(str(dist[0]) + ',' + str(dist[1]) + '\n')


time.sleep(90)
sc.parallelize(range(20)).foreach(process)
