"""

User-defined Reading and Preprocessing Data

For package scanpy (https://github.com/falexwolf/scanpy).

Written in Python 3 (compatible with 2) using Numpy, Scipy, Matplotlib
(use Anaconda package for all modules).
   
"""   

import warnings
from sys import path
# scientific modules
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib
with warnings.catch_warnings():
    # import of pyplot raises warning only in python 2
    warnings.simplefilter('ignore')
    import matplotlib.pyplot as pl
# scanpy, first try loading it locally
path.insert(0,'.')
import scanpy.preprocess

# dictionary storing examples
examples = {

# this is the same data as in scanpy, but 
# without a log before transform
'paul15_nolog' : {
    'doi' : '10.1016/j.cell.2015.11.013',
    'datafile' : 'data/paul15/paul15.h5',
    'dpt' : { 'k' : 10 } # increase number of neighbors (default 5) 
    }

}

def paul15_nolog():
    """
    See scanpy.preprocess.paul15 for a version that uses
    logarithmized data as in Haghverdi et al. (2016).
    """

    ddata = scanpy.preprocess._paul15_raw()

    # set root cell as in Haghverdi et al. (2016)
    # subtract 1, as in Haghverdi this is a matlab style index
    xroot = ddata['X'][840]

    # write to data dictionary
    ddata['xiroot'] = xroot
    
    return ddata

def get_data(exkey):
    """ 
    Retrieve data from file stored in the examples dictionary.

    This is almost the same function as scanpy.preprocess.get_data

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
            Array storing the names of columns (genes).

    """

    example = examples[exkey]

    # dictionary of all globally defined functions to preprocess
    # data for examples etc.
    functions = globals()
    if exkey in functions:
        ddata = functions[exkey]()
    else:
        raise ValueError('do not know how to read data for',exkey,
                         'please generate an entry in the preprocess.examples',
                         'dictionary and define a function for reading data',
                         'with the same name')

    return example, ddata

