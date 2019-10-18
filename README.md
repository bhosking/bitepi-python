### bitepi

bitepi is a python wrapper around the BitEpi project found at [https://github.com/aehrc/BitEpi](https://github.com/aehrc/BitEpi). It provides a pandas
interface for identification of epistasis interactions.
It exposes a single **Epistasis** class, through which the analysis can be
performed by calling **compute_epistasis**.

Input is two arrays, which can be lists, numpy arrays or pandas dataframes.

sample_array contains mappings of sample names to case (1) or control (0).
Note that the header is ignored for numpy arrays and pandas dataframes, and
should not be present in python lists.

sample | case/control
--- | ---
S1 | 0
S2 | 1
S3 | 1
S4 | 0
S5 | 0

genotype_array contains the genotypes of each sample at each SNP, with a 0, 1
and 2 representing 0|0, 0|1 and 1|1 respectively. Headers are used to match
samples to the sample_array, but the first column's header is ignored.

SNP | S1 | S2 | S3 | S4 | S5
---|---|---|---|---|---
snpA | 0 | 0 | 2 | 1 | 0
snpB | 2 | 1 | 2 | 1 | 2
snpC | 0 | 1 | 1 | 2 | 1
snpD | 0 | 2 | 2 | 2 | 1
snpE | 1 | 1 | 1 | 0 | 1

The sets of samples do not need to match exactly, unless **Epistasis** is
called with *strict_intersect=True*. If the sample sets do not match, analysis is done on the intersect.

The output will be a dictionary, with metric codes e.g. "IG.1" as the keys and
pandas dataframes as the values.

```python
import bitepi
sample_array = [
    ['S1', 0],
    ['S2', 1],
    ['S3', 1],
    ['S4', 0],
    ['S5', 0],
]
genotype_array = [
    ['SNP', 'S1', 'S2', 'S3', 'S4', 'S5'],
    ['snpA', 0, 0, 2, 1, 2],
    ['snpB', 2, 1, 2, 1, 2],
    ['snpC', 0, 0, 1, 1, 1],
    ['snpD', 2, 2, 1, 2, 1],
    ['snpE', 0, 1, 1, 1, 2],
]
epistasis = bitepi.Epistasis(
    genotype_array=genotype_array,
    sample_array=sample_array,
)
interactions = epistasis.compute_epistasis(
    sort=True,
    best_ig=True,
)['best_ig']

print(interactions)
```

This should return:

```
    SNP     SNP_P    PAIR_P  TRIPLET_P  QUADLET_P    SNP_IG   PAIR_IG  TRIPLET_IG  QUADLET_IG  PAIR TRIPLET_1 TRIPLET_2 QUADLET_1 QUADLET_2 QUADLET_3
0  snpA  1.707692  2.000000   2.109091        0.0  1.187692  0.266667    0.109091         0.0  snpE      snpB      snpE      snpA      snpA      snpA
1  snpB  1.642424  1.909091   2.109091        0.0  1.122424  0.266667    0.200000         0.0  snpC      snpC      snpE      snpA      snpA      snpA
2  snpC  1.641026  1.909091   2.109091        0.0  1.121026  0.266667    0.200000         0.0  snpB      snpD      snpE      snpA      snpA      snpA
3  snpD  1.642424  1.909091   2.109091        0.0  1.122424  0.175758    0.200000         0.0  snpE      snpC      snpE      snpA      snpA      snpA
4  snpE  1.733333  2.000000   2.109091        0.0  1.213333  0.266667    0.200000         0.0  snpA      snpC      snpD      snpA      snpA      snpA

```

For higher order interactions (p3, ig3, p4 and ig4) Epistasis may take several
minutes to run, depending on the number of SNPs. If more information is
required when running, the logging level can be increased to `logging.INFO`
or `logging.DEBUG`. `logging.DEBUG` will provide the greatest detail, including
logging from within the binary.

```python
import logging
logging.root.setLevel(logging.DEBUG)
```
