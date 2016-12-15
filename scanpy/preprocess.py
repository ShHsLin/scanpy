"""

Reading and Preprocessing Data

From package scanpy (https://github.com/falexwolf/scanpy).

The attribute 'examples' is a dictionary and stores all available examples. Furthermore,
for each example, there is a function that retrieves the data.

Written in Python 3 (compatible with 2) using Numpy, Scipy, Matplotlib
(use Anaconda package for all modules).

Copyright (c) 2016 F. Alexander Wolf (http://falexwolf.de).
   
"""   

import os
# scientific modules
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
from collections import OrderedDict as odict
from string import ascii_lowercase
# scanpy modules
import scanpy.utils as utils
import scanpy.logfile as log

# dictionary storing examples, entries are sorted by categories
# simulated/experimental and within each category alphabetically
examples = {

# simulated data
'krumsiek11' : {
    'doi' : '10.1371/journal.pone.0022649',
    'datafile' : 'data/krumsiek11/sim_000000.txt',
    'dpt' : { 
        'num_branchings' : 2, # detect two branching points (default 1)
        'allow_branching_at_root' : True } # allow branching directly at root         
    },
'toggleswitch' : {
    'datafile' : 'data/toggleswitch/sim_000000.txt',
    },

# experimental data
'moignard15' : {
    'doi' : '10.1038/nbt.3154',
    'datafile' : 'data/moignard15/nbt.3154-S3.xlsx',
    },
'nestorowa16' : {
    'doi' : '10.1182/blood-2016-05-716480',
    'datafile' : 'data/nestorowa16/mouse_bone_marrow_data/normalisedCountsVariableGenes.csv',
    'dpt' : { 'k' : 10 } # increase number of neighbors (default 5)
    },
'paul15' : {
    'doi' : '10.1016/j.cell.2015.11.013',
    'datafile' : 'data/paul15/paul15.h5',
    'dpt' : { 
        'k' : 20, # increase number of neighbors (default 5)
        'fix_nr_neighbors' : True } # set a hard threshold on number of neighbors        
    }
}

def moignard15():
    """ 
    Get data matrix, labels, and xroot for data of Moignard et al. (2015).

    Returns
    -------
    ddata : dict containing
        X : np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        xiroot : np.ndarray or int
            Expression vector or index of root cell.
        rownames : np.ndarray
            Array storing the experimental labels of samples.
        colenames : np.ndarray
            Array storing the names of genes.
    """

    example = examples['moignard15']

    filename = example['datafile']

    df = utils.read_excel_hdf5_to_pandas(filename,sheet='dCt_values.txt')
    X = df.values    

    # the 0th column contains the labels and is omitted.  also omit the 5th
    # column (Eif2b1), the 32nd (Mrpl19), the 37th (Polr2a) and the 46th
    # (last,UBC), as done by Haghverdi et al. (2016)
    genes = np.r_[np.arange(1,5),np.arange(6,32),np.arange(33,37),np.arange(38,46)]
    log.m(0,'selected',len(genes),'genes')
    X = X[:,genes]

    # choose root cell with expression vector xroot
    # this is the choice in Haghverdi et al. (2016)
    # note that in Matlab/R, counting starts at 1
    if True:
        xroot = X[532] 

    # compute xroot as average over all stem cells
    if False:    
        xroot = X[1557:2181].mean(axis=0)

    # labels
    labels = df.iloc[:,0].values.astype(str)
    # gene names
    genenames = np.array(df.columns[genes],dtype=str)
    
    ddata = {
        'X' : X, 'xiroot' : xroot, 'rownames' : labels, 'colnames' : genenames
    }

    return ddata

def _moignard15_colors_labels(labels):
    """
    Return the colors for labels from Moignard et al. (2015).

    Parameters
    ----------
    labels : list of str
        List of strings storing labels.

    Returns
    -------
    colors_labels : list
        List storing the colors associated with labels.
    """
    # coloring according to Moignard 2015 cell groups
    colors_labels = []
    if labels.ndim == 2:
        labels = labels[:,0]        
    for l in labels:
        if 'HF' == l[:2]: colors_labels.append('#D7A83E')
        elif 'NP' == l[:2]: colors_labels.append('#7AAE5D')
        elif 'PS' == l[:2]: colors_labels.append('#497ABC')
        elif '4SG' == l[:3]: colors_labels.append('#AF353A')
        elif '4SFG' == l[:4]: colors_labels.append('#765099')

    colors_labels = np.array(colors_labels)

    return colors_labels

def nestorowa16():
    """ 
    Get data matrix, labels, and xroot for data of Nestorowa et al. (2016).

    Returns
    -------
    See moignard15.
    """

    example = examples['nestorowa16']

    filename = example['datafile']

    df = pd.read_csv(filename)

    # labels in first row (column labels)
    labels = np.array(df.columns[1:],dtype=str) 
    # gene names in first column
    genenames = df.iloc[:,0].values.astype(str)
    # abbreviate genenames
    for ign,gn in enumerate(genenames):
        genenames[ign] = gn.replace('ENSMUSG','')
    # values 
    X = df.values[:,1:].T.astype(float)

    # take logarithm as for Paul et al. (2015)
    if False:
        X = np.log(X+1)

    # some testing output
    if False:
        print(len(labels))
        print(labels)
        print(len(genenames))
        print(genenames)
        print(X.shape)
        print(X)

    xroot = X[920]

    ddata = {
        'X' : X, 'xiroot' : xroot, 'rownames' : labels, 'colnames' : genenames
    }

    return ddata

def paul15():
    """ 
    Get data matrix, labels, and xroot for data of Moignard et al. (2015).

    This largely follows a tutorial by Maren Buttner.

    Returns
    -------
    See moignard15.
    """

    ddata = _paul15_raw()

    X = ddata['X']

    # apply logarithm to plain count data, shifted by one to avoid 
    # negative infinities, as done in Haghverdi et al. (2016)
    X = np.log(X+1)

    # set root cell as in Haghverdi et al. (2016)
    # subtract 1, as in Haghverdi this is a matlab style index
    xroot = X[840]
    
    # write to data dictionary
    ddata['X'] = X
    ddata['xiroot'] = xroot
    
    return ddata

def _paul15_raw():
    """
    Helper function for paul15.
    """

    example = examples['paul15']

    filename = example['datafile']

    ddata = utils.read_hdf5(filename,'data.debatched')
    # data set has automatically been transposed
    X = ddata['X']
    labels = ddata['rownames']
    genenames = ddata['colnames']

    # cluster assocations identified by Paul et al.
    if False:
        clusterid = utils.read_hdf5(filename,'cluster.id')['X']

    # subtract 1, as these are R style idcs
    infogenenames_subidcs = utils.read_hdf5(filename,'info.genes_codes')['X']-1
    infogenenames = utils.read_hdf5(filename,'info.genes_strings')['X']
    log.m(1,'the first 10 informative gene names are')
    log.m(1,infogenenames[:10])

    # just keep the first of the equivalent names for each gene
    genenames = np.array([gn.split(';')[0] for gn in genenames])
    log.m(1,'the first 10 trunkated gene names are')
    log.m(1,genenames[:10])

    # mask array for the informative genes
    infogenes_idcs = np.array([(True if gn in infogenenames else False)
                               for gn in genenames])
    
    # restrict data array to the 3451 informative genes
    X = X[:,infogenes_idcs]
    genenames = genenames[infogenes_idcs]
    log.m(1,'after selecting info genes, the first 10 gene names are')
    log.m(1,genenames[:10])

    ddata = {
        'X' : X, 'rownames' : labels, 'colnames' : genenames
    }

    return ddata

def get_data(exkey):
    """ 
    Retrieve data from file stored in the examples dictionary.

    Parameters
    ----------
    exkey : str
        String that identifies the example in the dictionary examples.

    Returns
    -------
    example : dict
        Dictionary storing meta-information about data.
    ddata : dict containing
        X : np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        xroot : np.ndarray
            Expression vector of root cell.
        rownames : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        colnames : np.ndarray
            Array storing the names of columns (gene names).
    """

    example = examples[exkey]

    # dictionary of all globally defined functions to preprocess
    # data for examples etc.
    functions = globals()
    if exkey in functions:
        ddata = functions[exkey]()
    # read formatted plain text file
    elif example['datafile'][-4:] == '.txt':
        ddata = utils.read_formatted_text(example['datafile'])
    else:
        raise ValueError('do not know how to read data for',exkey,
                         'please generate an entry in the preprocess.examples',
                         'dictionary and define a function for reading data',
                         'with the same name')

    return example, ddata


