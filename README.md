# fcup-abi
[![CircleCI](https://circleci.com/gh/antonioalmeida/fcup-abi.svg?style=svg)](https://circleci.com/gh/antonioalmeida/fcup-abi)

🍀 **bioseq** - a Python module to handle biological sequences and perform common Bioinformatics operations. Developed as part of 'Algorithms for Bioinformatics', a subject @FCUP.

### Reports
The package's development was divided in 3 stages, each with its own report.

- [Basic Processing of Biological Sequences](https://github.com/antonioalmeida/fcup-abi/blob/master/reports/basic-processing-of-biological-sequences.pdf)
- [Pairwise Sequence Alignment](https://github.com/antonioalmeida/fcup-abi/blob/master/reports/pairwise-sequence-alignment.pdf)
- [Bioinformatics Pipeline](https://github.com/antonioalmeida/fcup-abi/blob/master/reports/bioinformatics-pipeline.pdf)

### Requirements

- Python3

### Install
```shell
$ pip3 install git+https://github.com/antonioalmeida/fcup-abi
```

### Running Demos

```shell
cd demos
$ python3 run_me_<n>.py
```

### Running Tests

```shell
$ python3 -m unittest test/*.py
```

### Usage

Import the package
```python
# sample.py

import bioseq as bs
```

Use the functions as demonstrated on the [demos](https://github.com/antonioalmeida/fcup-abi/tree/master/demos) and the test suite.
