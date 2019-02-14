import os
import sys
import random
import pickle
import numpy as np
from timeit import default_timer as timer
from mhfp.encoder import MHFPEncoder
from mhfp.lsh_forest import LSHForest
from rdkit.Chem import AllChem

# config = mstmap.LayoutConfiguration()
# # config.merger = mstmap.Merger.Solar
# # print(config)

# # TODO: Fails for disconnected components!
# u = mstmap.VectorUint([0, 1, 2, 3, 4])
# v = mstmap.VectorUint([1, 2, 0, 4, 3])
# w = mstmap.VectorFloat([1.0, 1.0, 1.0, 2.0, 6.0])
# x, y = mstmap.layout(5, u, v, config, w)

# print(x)
# print(y)

enc = MHFPEncoder(512)



fps = []

if not os.path.isfile('fps.dat'):
    with open('drugbank.smi', 'r') as f:
        i = 0
        for line in f:
            smiles = line.split()[0].strip()
            mol = AllChem.MolFromSmiles(smiles)
            if mol:
                fps.append(enc.encode_mol(mol))
            i += 1
            if i > 2000: break
    pickle.dump(fps, open('fps.dat', 'wb'))
else:
    fps = pickle.load(open('fps.dat', 'rb'))


m = enc.encode("N=C(N)NCCC[C@H](NC(=O)[C@@H]1CCCN1C(=O)[C@H](N)Cc1ccccc1)C(=O)N1CCC[C@H]1C(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)N[C@@H](CC(=O)N)C(=O)NCC(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CC(C)C)C(=O)O")
n = enc.encode("O=C(N1[C@@H](CCC1)C(=O)NNC(=O)N)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H]1NC(=O)CC1)Cc1[nH]cnc1)Cc1c2c([nH]c1)cccc2)CO)Cc1ccc(O)cc1)COC(C)(C)C)CC(C)C)CCCN=C(N)N")
q = enc.encode("[C@@H](C(=O)NCC(=O)N[C@@H](C(=O)N[C@@H](C)C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)NCCO)Cc1c2c(cccc2)[nH]c1)CC(C)C)Cc1c[nH]c2c1cccc2)CC(C)C)Cc1c[nH]c2c1cccc2)CC(C)C)Cc1c[nH]c2c1cccc2)C(C)C)C(C)C)C(C)C)CC(C)C)(C(C)C)NC=O")
r = enc.encode('CNCNCNCNC')

lf_classic = LSHForest(512, 64)
for i, e in enumerate(fps):
    lf_classic.add(i, e)
lf_classic.index()

start = timer()
result = lf_classic.query(r, 5)
end = timer()

print(end - start)

print(result)















# import os
# import sys
# import random
# import pickle
# from timeit import default_timer as timer
# from mhfp.encoder import MHFPEncoder
# from mhfp.lsh_forest import LSHForest
# from mstmap import mstmap
# from rdkit.Chem import AllChem

# # config = mstmap.LayoutConfiguration()
# # # config.merger = mstmap.Merger.Solar
# # # print(config)

# # # TODO: Fails for disconnected components!
# # u = mstmap.VectorUint([0, 1, 2, 3, 4])
# # v = mstmap.VectorUint([1, 2, 0, 4, 3])
# # w = mstmap.VectorFloat([1.0, 1.0, 1.0, 2.0, 6.0])
# # x, y = mstmap.layout(5, u, v, config, w)

# # print(x)
# # print(y)


# fps = []
# enc = MHFPEncoder(128)

# if not os.path.isfile('fps.dat'):
#     with open('drugbank.smi', 'r') as f:
#         i = 0
#         for line in f:
#             smiles = line.split()[0].strip()
#             mol = AllChem.MolFromSmiles(smiles)
#             if mol:
#                 fps.append(enc.encode_mol(mol))
#             i += 1
#             if i > 3: break
#     pickle.dump(fps, open('fps.dat', 'wb'))
# else:
#     fps = pickle.load(open('fps.dat', 'rb'))


# m = enc.encode("CNCNCNCNCNCNC")
# n = enc.encode("CCCCCCCCCCCCC")
# q = enc.encode("CCCCCCCCCCCCO")

# start = timer()
# lf_classic = LSHForest(128, 8)
# for i, fp in enumerate(fps):
#     lf_classic.add(i, m)
# lf_classic.index()
# end = timer()
# print(end - start)

# start = timer()
# for _ in range(1000):
#     result = lf_classic.query(q, 5)
# end = timer()
# print(end - start)

# print(result)

# start = timer()
# lf = mstmap.LSHForest()
# for i, fp in enumerate(fps):
#     lf.add(i, mstmap.VectorUint(m))
# lf.index()
# end = timer()
# print(end - start)

# start = timer()
# for _ in range(1000):
#     result = lf.query(mstmap.VectorUint(q), 5)
# end = timer()

# print(end - start)
# print(result)