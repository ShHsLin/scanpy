"""

Logfile, global variables, timing, verbosity

From package scanpy (https://github.com/falexwolf/scanpy).

Sets various global variables and provides convenient log message operations.

Copyright (c) 2016 F. Alexander Wolf (http://falexwolf.de).

Note
----

This is partly based on

http://stackoverflow.com/questions/1557571/how-to-get-time-of-a-python-program-execution
   
"""   

import atexit
from time import clock
from functools import reduce


verbosity = 1
""" Global verbosity level, choose from {0,...,6}. """

suffix = ''
""" Global suffix, which is appended to basename of output and figure files.
"""

extd = 'h5'
""" Global file extension format for data storage. 

Allowed are 'h5' (hdf5), 'xlsx' (Excel) or 'csv' (comma separated value
file).
"""

extf = 'pdf'
""" Global file extension for saving figures.

Recommended are 'png' and 'pdf'. Many other formats work as well.
"""

recompute = False
""" Don't use the results of previous calculations.

Recompute and overwrite previous files.  
"""

savefigs = False
""" Save plots/figures as files in directory 'figs'.

Do not show plots/figures interactively.
"""

def filename(filename=''):
    """ 
    Define filename of logfile. 

    If not defined, log output will be to the standard output.

    Parameters
    ----------
    filename : str
        Filename of 
    """
    global logfilename
    logfilename = filename
    # if providing a logfile name, automatically set verbosity to a very high level
    verbosity(5)

def m(v=0,*msg):
    """ 
    Write message to log output, depending on verbosity level.

    Parameters
    ----------
    v : int
        Verbosity level of message.
    *msg : 
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    if verbosity > v:
        mi(*msg)

def mi(*msg):
    """ 
    Write message to log output, ignoring the verbosity level.

    Parameters
    ----------
    *msg : 
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    if logfilename == '':
        # in python 3, the following works
        # print(*msg)
        # due to compatibility with the print statement in python 2 we choose
        print(' '.join([str(m) for m in msg]))
    else:
        out = ''
        for s in msg:
            out += str(s) + ' '
        with file(filename) as f:
            f.write(out + '\n')

def mt(v=0,*msg):
    """ 
    Write message to log output and show computation time.

    Both depending on chosen verbosity level.

    Parameters
    ----------
    v : int
        Verbosity level of message.
    *msg : str
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    if verbosity > v:
        global intermediate
        now = clock()
        elapsed_since_start = now - start
        elapsed = now - intermediate
        intermediate = now
        mi(_sec_to_str(elapsed),'-',*msg)
        # mi(separator)

def aa(p):
    """ 
    Add global arguments that serve as global variables.

    Parameters
    ----------
    p : ArgumentParser
        From the argparse module.

    Returns
    -------
    args : dict
        Dictionary of non-global variables.
    """
    aa = p.add_argument
    aa('--suffix',
       type=str,default='',metavar='SUF',
       help='Specify suffix to append to file basenames'
            ' (default: %(default)s).')
    aa('--subsample',
       type=int,default=1,metavar='SS',
       help='Specify subsample > 1 if you want to use a fraction of 1/subsample'
            ' of the data'
            ' (default: %(default)d).')
    aa('--recompute',
       action='store_const', default=False, const=True,
       help='Do not use the results of previous calculations.' 
            ' Recompute and overwrite previous files'
            ' (default: %(default)d).')
    aa('--savefigs',
       type=str,default='',const='pdf',nargs='?',metavar='EXTF',
       help='Save figures to files, do not show them as interactive plots.'
            ' You can specify the file by providing the extension, e.g. "pdf"'
            ' or "png". Just specifying --savefigs will save to "pdf"'
            # here, the standard initialization of default is replaced
            # as we need to check whether the user has actually supplied 
            # something below
            ' (default: do not save figures).')
    aa('--extd',
       type=str,default='h5',
       help='specify file extension and by that'
            ' choose format for saving results, either "h5", "xlsx" or "csv"'
            ' h5 is much faster and needs much less storage than xlsx and csv'
            ' (default: %(default)s).')
    aa('--log',
       action='store_const', default=False, const=True,
       help='write logfile to out/example..._log.txt'
            ' (default: %(default)d)')
    aa('-v','--verbosity',
       type=int,default=1,metavar='V',
       help='specify 6 >= VERBOSITY > 1 for more output, and 0 for no output'
            ' (default: %(default)d).')

    args = vars(p.parse_args())

    global suffix
    suffix = args['suffix']
    args.pop('suffix')

    global recompute
    recompute = args['recompute']
    args.pop('recompute')

    global verbosity
    verbosity = args['verbosity']
    args.pop('verbosity')

    global savefigs
    global extf
    if args['savefigs'] == '':
        savefigs = False
    else:   
        savefigs = True
        extf = args['savefigs']
    args.pop('savefigs')

    global extd
    extd = args['extd']
    args.pop('extd')


    return args

def _sec_to_str(t):
    """ 
    Format time in seconds.

    Parameters
    ----------
    t : int
        Time in seconds.
    """
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll,b : divmod(ll[0],b) + ll[1:],
            [(t*1000,),1000,60,60])

def _terminate():
    """ 
    Function called when program terminates.
    
    Similar to mt, but writes total runtime.
    """
    if verbosity > 0:
        now = clock()
        elapsed_since_start = now - start
        mi(separator)
        mi(_sec_to_str(elapsed_since_start),'- total runtime')
        mi(separator)


# further global variables
start = clock() # clock() is deprecated since version python version 3.3
intermediate = start
logfilename = ''
separator = 80*"-"

atexit.register(_terminate)
