import time
import random
import pandas as pd
import numpy as np
from operator import itemgetter
from rdkit.Chem import AllChem
from scipy.spatial.distance import jaccard

# Let's get sparkly
from pyspark import SparkContext, SparkConf
sparkConf = SparkConf().setAppName('ChEMBL Benchmark: Subsetting')
sc = SparkContext(conf=sparkConf)


def process(seed):
    smis = []

    with open('/cluster/chembl/chembl.smi') as f:
        for line in f.readlines():
            smis.append(line)

    random.seed(seed)
    subset = random.sample(smis, 10000)

    with open('/cluster/chembl/chembl.' + str(seed) + '.smi', 'w+') as f:
        for line in subset:
            f.write(line)

    


time.sleep(120)
sc.parallelize(range(20)).foreach(process)
