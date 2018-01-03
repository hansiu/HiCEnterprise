# HiCEnterprise
Scripts for prediction of long-range interactions between regions/domains based on Hi-C maps.
Establish contact!

## What is in the package?
HiCEnterprise takes a list of regions/domains and searches for their significant interactions on given Hi-C maps.
Package consists of two types of analyses: regions and domains.

### Regions
This part of the package is an implementation of the method for identification of long-range interacting regions and 
creating interaction profiles based on Hi-C data as introduced in (Won, 2016).  

The significance of the individual pairwise interaction between bins is calculated as the probability of observing a 
stronger contact under the fitted Weibull distribution matched by chromosome and distance. FDR (False Discovery Rate) is
used for correcting for multiple testing.

### Domains
In the domain part of the package the significance of interaction between domains A and B is calculated as the 
probability of observing a stronger contact under the fitted Hypergeometric distribution, where N (population size) is
the sum of all frequencies in the Hi-C map, a is the sum from the column of domain A and b is the sum from the row of
domain B (number of success states in the population and number of draws). 

## Basic usage

`HiCEnterprise type options`

>type could be either regions or domains, options are type dependent.

If you want to see the list of available arguments with short explanations type either `HiCEnterprise regions --help` or `HiCEnterprise domains --help`. Full specifications for each argument can be found in the documentation (docs folder or online).

### Installation

`python setup.py install`

> You may need administrative privileges to do this.

#### Testing the installation

To test if the scripts behave in a way they should be find the main HiCEnterprise directory and run:

```
python setup.py test
```

> pytest package is needed for running the test framework.


### Examples

We provided some example files (used for testing too) with which you can learn to use the program.
Here you can see two example runs:

#### Regions

`HiCEnterprise regions ` <TODO add options>

#### Domains

`HiCEnterprise domains ` <TODO add options>

Additional examples and specifications for each argument can be found in the documentation (docs folder or online).




# 

Hania Kranas, 2017
