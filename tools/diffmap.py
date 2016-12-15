"""

Just a wrapper that calls scanpy.diffmap

From package scanpy (https://github.com/falexwolf/scanpy).

Written in Python 3 (compatible with 2) using Numpy, Scipy, Matplotlib
(use Anaconda package for all modules).

Copyright (c) 2016 F. Alexander Wolf (http://falexwolf.de).
   
"""   

from sys import path
# my examples
import preprocess
# scanpy, first try loading it locally
path.insert(0,'.')
import scanpy as sc
import scanpy.dpt
import scanpy.utils
import scanpy.preprocess

if __name__ == '__main__':

    # concatenate my examples and the existing scanpy examples
    all_examples = sc.utils.merge_dicts(preprocess.examples,
                                        sc.preprocess.examples)

    args = sc.utils.read_args(sc.dpt.__doc__,all_examples)

    # add an argument indicating that the example is from the tools directory
    if args['exkey'] in preprocess.examples:
        args['preprocess'] = preprocess

    sc.dpt.main_diffmap(args)

    


