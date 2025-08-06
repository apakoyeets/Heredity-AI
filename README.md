# Heredity

This project implements a probabilistic AI to assess the likelihood that individuals in a family carry a particular gene and exhibit an associated trait — based on inheritance and mutation probabilities modeled by a Bayesian network.

## Overview

Mutations in the **GJB2** gene are one of the primary causes of hearing impairment in newborns. Each person carries two copies of this gene and can have 0, 1, or 2 copies of the mutated form. This project uses information about individuals and their parents to estimate the probability distributions for gene inheritance and expression of the trait.

Given a CSV file with family data, the program infers — for each person — the probabilities that they:
- Have 0, 1, or 2 copies of the gene
- Exhibit or do not exhibit the trait

### Example

```bash
$ python heredity.py data/family0.csv
Harry:
  Gene:
    2: 0.0092
    1: 0.4557
    0: 0.5351
  Trait:
    True: 0.2665
    False: 0.7335
James:
  Gene:
    2: 0.1976
    1: 0.5106
    0: 0.2918
  Trait:
    True: 1.0000
    False: 0.0000
Lily:
  Gene:
    2: 0.0036
    1: 0.0136
    0: 0.9827
  Trait:
    True: 0.0000
    False: 1.0000
