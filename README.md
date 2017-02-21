# A Kernel Independence Test for Geographical Language Variation.

Code for the following paper:

```
D. Nguyen and J. Eisenstein. A Kernel Independence Test for Geographical Language Variation. 
To appear in Computational Linguistics.
```

## Getting Started

The code is in Python 2.7 and makes use of several Python packages:

* descartes
* fiona
* numpy
* pyproj
* scipy
* shapely


### Tests

The following file runs some unit tests

```
python testing.py
```


## Data

The data can be downloaded from http://www.dongnguyen.nl/data/dataset-nguyen-eisenstein-cl2017.zip (274 MB)

* *synthetic_experiments*: The synthetic datasets. Each directory corresponds to one experiment. Each directory contains a results.txt file with the raw results. 
* *shapefiles*: Shapefiles of the Netherlands for plotting the synthetic datasets, aggregating data into bins (for Moran's I), etc.  You'll still need them if you would like to experiment with the synthetic data. 

## Experiments

* *synthetic_data.ipynb*: This notebook plots several selected synthetic datasets and shows how to apply the methods to the different types of data (binary, categorical and frequency data).
* *HSIC.ipynb*: This notebook illustrates HSIC with several synthetic (non-geographical) datasets.
* *plots.r*: Shows how to generate the plots in the paper based on the result files in the *synthetic_experiments* directory.
 
## Command line tool

```
python hsic_wrapper.py 
```

Frequency data (should return 0.00653):

```
python hsic_wrapper.py -d freq -l sample_data/locs1.txt -f sample_data/data1.txt
```

Binary data (should return 0.00241):

```
python hsic_wrapper.py -d bin -l sample_data/locs2.txt -f sample_data/data2.txt
```

Categorical data (should return 0.00111):

```
python hsic_wrapper.py -d cat -l sample_data/locs3.txt -f sample_data/data3.txt
```


## Authors

* **Dong Nguyen** - http://www.dongnguyen.nl
* **Jacob Eisenstein** - http://www.cc.gatech.edu/~jeisenst/

Contact: dong.p.ng@gmail.com