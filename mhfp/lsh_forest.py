import warnings
from operator import itemgetter
from collections import defaultdict
from mhfp.encoder import MHFPEncoder

# Adapted from:
# https://github.com/ekzhu/datasketch/blob/master/datasketch/lshforest.py  

class LSHForestHelper():
  """A class simplyfing the use of the LSH Forest implementation. Adding paramter kc to imporve search results.
  """
  def __init__(self, dims = 2048, n_prefix_trees = 64):
    """The parameter dims has to be set to the same number with which the MHFP fingerprints were encoded.

    Keyword arguments:
        dims {int} -- The number of dimensions the MHFP fingerprints to be indexed were encoded in (default: {2048})
        n_prefix_trees {int} -- The number of prefix trees to use for indexing (default: {64})

    """

    self.dims = dims
    self.n_prefix_trees = n_prefix_trees

    self.lsh_forest = LSHForest(self.dims, n_prefix_trees=self.n_prefix_trees)

  def add(self, key, mhfp):
    """Adds a MHFP stored in a numpy.ndarray to the LSH Forest. The index() method has to be run before
    entries can be queried.

    Keyword arguments:
        key -- The key for the MHFP entry.
        mhfp {numpy.ndarray} -- A MHFP array containing 32-bit hashes
    """

    self.lsh_forest.add(key, mhfp)

  def index(self):
    """Indexes the added entities. This method has to be run after all MHFPs have been added or when
    new ones were added. Has to be run before querying.
    """

    self.lsh_forest.index()

  def query(self, query_mhfp, k, data, kc=10):
    """ Query the LSH Forest. k is multiplied by kc when querying the LSH Forest. The k nearest neighbors
    ist then selected from the list of k * kc nearest neighbors using linear scan.

    Keyword arguments:
        query_mhfp {numpy.ndarray} -- The query MHFP fingerprint.
        k {int} -- The number of nearest neighbors to be returned from the approximate nearest neighbors
        data {dict} -- The MHFP values indexed with the same key supplied to add()
        kc {int} -- The factor in kc * k for additional approximate nearest neighbors to be searched for (default: {10})
    Returns:
      numpy.ndarray -- An array containing the merged hash values.
    """

    results = self.lsh_forest.query(query_mhfp, k * kc)
    return LSHForestHelper._get_knn(query_mhfp, results, k, data)

  @staticmethod
  def _get_knn(query_mhfp, ann, k, data):
    """ Brute-force search for selecting k nearest neighbors from k * kc  approximate nearest neighbors.

    Keyword arguments:
        query_mhfp {numpy.ndarray} -- The query MHFP fingerprint.
        ann {list} -- A list of indices of approximate nearest neighbors of size k * kc to be brute-force searched
        k {int} -- The number of nearest neighbors to be returned from the approximate nearest neighbors
        data {dict} -- The MHFP values indexed with the same key supplied to add()
    """

    dists = []

    for index in ann:
        dists.append((index, 1.0 - MHFPEncoder.distance(query_mhfp, data[index])))
    
    dists.sort(key=itemgetter(1), reverse=True)
    return [x[0] for x in dists[:k]]


class LSHForest():
  """A class for indexing and querying MHFP fingerprints using the LSH Forest algorithm.
  """

  def __init__(self, dims = 2048, n_prefix_trees = 64):
    """The parameter dims has to be set to the same number with which the MHFP fingerprints were encoded.

    Keyword arguments:
        dims {int} -- The number of dimensions the MHFP fingerprints to be indexed were encoded in (default: {2048})
        n_prefix_trees {int} -- The number of prefix trees to use for indexing (default: {64})
    """

    self.dims = dims
    self.n_prefix_trees = n_prefix_trees

    self.max_depth = int(self.dims / self.n_prefix_trees)
    self.hashtables = [defaultdict(list) for _ in range(self.n_prefix_trees)]
    self.hashranges = [(i * self.max_depth, (i + 1) * self.max_depth) for i in range(self.n_prefix_trees)]
    self.keys = dict()
    self.clean = False

    self.hashtables_sorted = [[] for _ in range(self.n_prefix_trees)]

  def add(self, key, mhfp):
    """Adds a MHFP stored in a numpy.ndarray to the LSH Forest. The index() method has to be run before
    entries can be queried.

    Keyword arguments:
        key -- The key for the MHFP entry.
        mhfp {numpy.ndarray} -- A MHFP array containing 32-bit hashes
    """

    # self.keys[key] = [self._swap(mhfp[start:end]) for start, end in self.hashranges]
    self.keys[key] = 0
    
    k = [self._swap(mhfp[start:end]) for start, end in self.hashranges]
    for h, hashtable in zip(k, self.hashtables):
      hashtable[h].append(key)
    
    self.clean = False

  def index(self):
    """Indexes the added entities. This method has to be run after all MHFPs have been added or when
    new ones were added. Has to be run before querying.
    """
    for i, hashtable in enumerate(self.hashtables):
      self.hashtables_sorted[i] = [k for k in hashtable.keys()]
      self.hashtables_sorted[i].sort()
    
    self.clean = True

  def query(self, mhfp, k):
    """Query the LSH Forest.

    Keyword arguments:
        mhfp -- The query MHFP fingerprint
        k -- The number of nearest neighbors to be searched for
    """
    if not self.clean:
      warnings.warn('Index out of date. Call index() before querying after adding new entities.', UserWarning)

    results = set()
    r = self.max_depth

    while r > 0:
      for key in self._internal_query(mhfp, r):
        results.add(key)
        if len(results) >= k:
          return list(results)
      r -= 1
    return list(results)

  def _internal_query(self, mhfp, r):
    """Internal method, use the query method instead.

    Keyword arguments:
        mhfp -- The query MHFP fingerprint
        r -- The current search depth
    """

    prefixes = [self._swap(mhfp[start:start + r]) for start, _ in self.hashranges]
    len_prefix = len(prefixes[0])

    for hashtable_sorted, hashtable, prefix in zip(self.hashtables_sorted, self.hashtables, prefixes):
      i = self._binary_search(len(hashtable_sorted), lambda x: hashtable_sorted[x][:len_prefix] >= prefix)
      
      if i < len(hashtable_sorted) and hashtable_sorted[i][:len_prefix] == prefix:
        j = i
        while j < len(hashtable_sorted) and hashtable_sorted[j][:len_prefix] == prefix:
          for key in hashtable[hashtable_sorted[j]]:
            yield key
          j += 1

  def _swap(self, hashes):
    """Internal method. Swaps bytes for prefixing.

    Keyword arguments:
        hashes -- The hashes on which to perform a byteswap
    """
    return bytes(hashes.byteswap().data)


  def _binary_search(self, n, func):
    """Internal method. Performs a binary search on a hashtable using func as the comparison function.

    Keyword arguments:
        n -- Upper bound for index
        func -- The comparison fuction
    """

    i = 0
    j = n

    while i < j:
      h = int(i + (j - i) / 2)
      if not func(h):
        i = h + 1
      else:
        j = h
    
    return i