import time
import datetime
import pandas as pd
import numpy as np
from rdkit.Chem import AllChem
from mhfp.encoder import MHFPEncoder
from mhfp.lsh_forest import LSHForestHelper
from sklearn.neighbors import NearestNeighbors
from operator import itemgetter

# Let's get sparkly
from pyspark import SparkContext, SparkConf
sparkConf = SparkConf().setAppName('ChEMBL Benchmark: Recovery')
sc = SparkContext(conf=sparkConf)


def process(subset):
    np.random.seed(42)
    samples = np.random.choice(10000, 20, replace=False)

    for fp in ['mhfp6', 'ecfp', 'mhecfp4']:
        fps = []

        if fp == 'mhfp6':
            with open('/cluster/chembl/chembl.' + str(subset) + '.mhfp6') as f:
                for line in f.readlines():
                    fps.append(np.array(list(map(int, line.split(','))), dtype=np.uint32))
        elif fp == 'ecfp4':
            with open('/cluster/chembl/chembl.' + str(subset) + '.ecfp4') as f:
                for line in f.readlines():
                    fps.append(np.array(list(map(int, line.split(','))), dtype=np.uint8))
        else:
            with open('/cluster/chembl/chembl.' + str(subset) + '.mhecfp4') as f:
                for line in f.readlines():
                    fps.append(np.array(list(map(int, line.split(','))), dtype=np.uint32))

        nn_linear = {}

        if fp == 'mhfp6':
            for sample in samples:
                nn_linear[sample] = []
                with open('/cluster/chembl/dists_mhfp_' + str(subset) + '_' + str(sample) + '.txt') as f:
                    for line in f.readlines():
                        nn_linear[sample].append(int(line.split(',')[0]))
        elif fp == 'ecfp4':
            for sample in samples:
                nn_linear[sample] = []
                with open('/cluster/chembl/dists_ecfp_' + str(subset) + '_' + str(sample) + '.txt') as f:
                    for line in f.readlines():
                        nn_linear[sample].append(int(line.split(',')[0]))
        else:
            for sample in samples:
                nn_linear[sample] = []
                with open('/cluster/chembl/dists_mhecfp_' + str(subset) + '_' + str(sample) + '.txt') as f:
                    for line in f.readlines():
                        nn_linear[sample].append(int(line.split(',')[0]))

        
        f_perf = open('/cluster/chembl/perf_' + fp + '_' + str(subset) + '.txt', 'w+')

        # MHFP
        for kc in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            if fp == 'mhfp6':
                for npf in [8, 16, 32, 64, 128, 256, 512]:
                    # LHS Forest
                    t_from = datetime.datetime.now()
                    lf = LSHForestHelper(2048, n_prefix_trees=npf)
                    for i, e in enumerate(fps):
                        lf.add(i, e)
                    lf.index()
                    t_built = datetime.datetime.now()
                    build_time = (t_built - t_from).total_seconds() * 1000.0

                    for k in [5, 10, 50, 100]:
                        for sample in samples:
                            f_perf.write(str(sample) + ',' + str(len(fps)) + ',' + str(k) + ',' + str(kc) + ',' + 'lf' + ',' + str(npf) + ',')
                            f_perf.write(str(build_time) + ',')

                            t_from = datetime.datetime.now()
                            # Find k * 10 nearest neighbours to account for approximation
                            # then get the exact distances and take top k
                            nn = lf.query(fps[sample], k, fps, kc)
                            
                            t_queried = datetime.datetime.now()
                            f_perf.write(str((t_queried - t_from).total_seconds() * 1000.0) + ',')

                            nn = np.array(nn)
                            f_perf.write(str(len(np.intersect1d(nn_linear[sample][:k], nn))) + '\n')

            # MHECFP
            if fp == 'mhecfp4':
                for npf in [8, 16, 32, 64, 128, 256, 512]:
                    # LHS Forest
                    t_from = datetime.datetime.now()
                    lf = LSHForestHelper(2048, n_prefix_trees=npf)
                    for i, e in enumerate(fps):
                        lf.add(i, e)
                    lf.index()
                    t_built = datetime.datetime.now()
                    build_time = (t_built - t_from).total_seconds() * 1000.0
                    
                    for k in [5, 10, 50, 100]:
                        for sample in samples:
                            f_perf.write(str(sample) + ',' + str(len(fps)) + ',' + str(k) + ',' + str(kc) + ',' + 'lf_2' + ',' + str(npf) + ',')
                            f_perf.write(str(build_time) + ',')

                            t_from = datetime.datetime.now()
                            # Find k * 10 nearest neighbours to account for approximation
                            # then get the exact distances and take top k
                            nn = lf.query(fps[sample], k, fps, kc)
                            
                            t_queried = datetime.datetime.now()
                            f_perf.write(str((t_queried - t_from).total_seconds() * 1000.0) + ',')

                            nn = np.array(nn)
                            f_perf.write(str(len(np.intersect1d(nn_linear[sample][:k], nn))) + '\n')

        # kNN
        if fp == 'ecfp4':
            for k in [5, 10, 50, 100]:
                t_from = datetime.datetime.now()
                nbrs = NearestNeighbors(n_neighbors=k, algorithm='balltree').fit(fps)
                t_built = datetime.datetime.now()
                build_time = (t_built - t_from).total_seconds() * 1000.0

                for sample in samples:
                    f_perf.write(str(sample) + ',' + str(len(fps)) + ',' + str(k) + ',' + 'knn' + ',' + '0' + ',')
                    f_perf.write(str(build_time) + ',')

                    _, indices = nbrs.kneighbors([fps[sample]])
                    nn = indices[0]
                    
                    t_queried = datetime.datetime.now()
                    f_perf.write(str((t_queried - t_from).total_seconds() * 1000.0) + ',')

                    nn = np.array(nn)
                    f_perf.write(str(len(np.intersect1d(nn_linear[sample][:k], nn))) + '\n')


        f_perf.close()


time.sleep(90)

sc.parallelize(range(20)).foreach(process)