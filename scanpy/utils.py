"""

Utility functions and classes

From package scanpy (single-cell toolkit in Python).

Written in Python 3 (compatible with 2) using Numpy, Scipy, Matplotlib, Pandas
(use Anaconda package for all modules).

Copyright (c) 2016 F. Alexander Wolf (http://falexwolf.de).
   
"""   

# standard modules
import os
import argparse
import h5py
import sys
import warnings
# scientific modules
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
import scipy.cluster
import matplotlib
with warnings.catch_warnings():
    # import of pyplot raises strange warning only in python 2
    warnings.simplefilter('ignore')
    import matplotlib.pyplot as pl
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from collections import OrderedDict as odict
# local modules
import scanpy.plot as mpl_plot
import scanpy.logfile as log

sppars = matplotlib.figure.SubplotParams

def read_args(description,examples_dict):
    """
    Read and process arguments.

    Further arguments are read in scanpy.logfile.aa.
    """

    epilog = '\n'
    for k,v in examples_dict.items():
        epilog += '    ' + k + '\n'
    p = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=('shorthands for examples: possible values'+epilog))

    aa = p.add_argument
    aa('exkey',
       type=str, default='toggleswitch',
       help='Specify the "example key" (just a shorthand), which is used'
            ' to look up files and parameters. See possible values below'
            ' (default: %(default)s).')
    args = log.aa(p)

    if args['exkey'] not in examples_dict.keys():
        raise ValueError('Use one of possible values for EXAMPLE:' + epilog)
    
    if args['log']:
        prefix = args['exkey']
        sig = '_ss{:02}'.format(subsample)
        filename = prefix + sig + '_log.txt'
        log.filename(filename)

    return args

def merge_dicts(*dicts):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    Note
    ----
    http://stackoverflow.com/questions/38987/how-to-merge-two-python-dictionaries-in-a-single-expression
    """
    result = {}
    for d in dicts:
        result.update(d)
    return result

def scale_to_zero_one(x):
    """
    Take some 1d data and scale it so that min matches 0 and max 1.
    """
    
    xscaled = x - np.min(x)
    xscaled /= np.max(xscaled)
    
    return xscaled

def plot_zoom(ax,xy='x',factor=1):
    """ 
    Zoom into axis.

    Parameters
    ----------
    """
    limits = ax.get_xlim() if xy == 'x' else ax.get_ylim()
    new_limits = (0.5*(limits[0] + limits[1]) 
                  + 1./factor * np.array((-0.5, 0.5)) * (limits[1] - limits[0]))
    if xy == 'x':
        ax.set_xlim(new_limits)
    else:
        ax.set_ylim(new_limits)

def plot_scatter(ax,X,*args,**kwargs):
    """ 
    Plot scatter plot of data. Just some wrapper of matplotlib.Axis.scatter.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis to plot on.
    X : np.array
        Data array, data to be plotted needs to be in the first two columns.
    """
    if 's' not in kwargs:
        kwargs['s'] = 2 if X.shape[0] > 500 else 10
    if 'edgecolors' not in kwargs:
        kwargs['edgecolors'] = 'face'
    ax.scatter(X[:,0],X[:,1],**kwargs)
    ax.set_xticks([]); ax.set_yticks([])

def plot_scatter_list(Xs,layout='2d',c='blue',highlights=[],title=''):
    """ 
    Plot scatter plot of data.

    Parameters
    ----------
    Xs : np.ndarray or list of np.ndarray
        Single data array or list of data arrays. Rows store observations,
        columns store variables. For example X, or phi or [phi,psi]. Arrays must
        be of dimension ndim=2.
    c : string or np.array or list of np.array
        Colors array, for example [pseudotimes] or [pseudotimes,colors_labels].
    title : str
        Title of plot, default ''. 
    layout : str
        Choose from '2d', '3d' and 'unfolded 3d', default '2d'.

    Returns
    -------
    axs : matplotlib.axis or list of matplotlib.axis
        Depending on whether supplying a single array or a list of arrays,
        return a single axis or a list of axes.
    """
    # if we have a single array, transform it into a list with a single array
    avail_layouts = ['2d','3d','unfolded 3d']
    if layout not in avail_layouts:
        raise ValueError('choose layout from',avail_layouts)
    colors = c
    if type(Xs) == np.ndarray:
        Xs = [Xs]
    if len(colors) == len(Xs[0]) or type(colors) == str:
        colors = [colors]
    figsize = (4*len(colors),4*len(Xs))

    if layout == 'unfolded 3d':
        if len(Xs) != 1:
            raise ValueError('use single 3d array')
        if len(colors) > 1:
            raise ValueError('choose a single color')
        figsize = (4*2,4*2)
        X = Xs[0]
        Xs = [X[:,[1,2]], X[:,[0,1]], X, X[:,[0,2]]]
    fig = pl.figure(figsize=figsize, 
                        subplotpars=sppars(left=0.07,right=0.98,bottom=0.08))
    fig.suptitle(title)
    count = 1
    bool3d = True if layout == '3d' else False
    axs = []
    for X in Xs:
        markersize = 2 if X.shape[0] > 500 else 10
        for icolor,color in enumerate(colors):
            # set up panel
            if layout == 'unfolded 3d' and count != 3:
                ax = fig.add_subplot(2,2,count)
                bool3d = False
            elif layout == 'unfolded 3d' and count == 3:
                ax = fig.add_subplot(2,2,count,
                                     projection='3d')
                bool3d = True
            elif layout == '2d':
                ax = fig.add_subplot(len(Xs),len(colors),count)
            elif layout == '3d':
                ax = fig.add_subplot(len(Xs),len(colors),count,
                                     projection='3d')
            if not bool3d:
                data = X[:,0],X[:,1]
            else:
                data = X[:,0],X[:,1],X[:,2]
            # do the plotting
            ax.scatter(*data,
                       c=color,edgecolors='face',s=markersize,
                       cmap='jet')
            # set the label
            if len(colors) >= 2:
                if icolor == 0:
                    ax.set_title('pseudotime')
                if icolor == 1:
                    ax.set_title('segments')
                if icolor == 2:
                    ax.set_title('experimental labels')
            # output highlighted data points
            for iihighlight,ihighlight in enumerate(highlights):
                data = [X[ihighlight,0]],[X[ihighlight,1]]
                if bool3d:
                    data = [X[ihighlight,0]],[X[ihighlight,1]],[X[ihighlight,2]]
                ax.scatter(*data,c='black',
                           facecolors='black',edgecolors='black', 
                           marker='x', s=50)  
                # the following is a Python 2 compatibility hack
                ax.text(*([d[0] for d in data]+[str(ihighlight)]))
            ax.set_xticks([]); ax.set_yticks([])
            if bool3d:
                ax.set_zticks([]) 
            axs.append(ax)
            count += 1
    # scatter.set_edgecolors = scatter.set_facecolors = lambda *args:None
    return axs[0] if len(axs) == 1 else axs

def plot_diffmap(Xs,layout='2d',axlabels=None,**kwargs):
    """
    See plot_scatter_list for definition of arguments.

    Adds labels to axis.
    """
    axs = plot_scatter_list(Xs,layout=layout,**kwargs)
    if type(axs) is not list:
        axs = [axs]
    bool3d = True if layout == '3d' else False
    # define default axlabels
    if axlabels is None:
        if layout == '2d':
            axlabels = [['DC'+str(i) for i in idcs] 
                         for idcs in [[1,2] for iax in range(len(axs))]]            
        elif layout == '3d':
            axlabels = [['DC'+str(i) for i in idcs] 
                         for idcs in [[1,2,3] for iax in range(len(axs))]]
        elif layout == 'unfolded 3d':
            axlabels = [['DC'+str(i) for i in idcs] 
                         for idcs in [[2,3],[1,2],[1,2,3],[1,3]]]
    # set axlabels
    for iax,ax in enumerate(axs):
        if layout == 'unfolded 3d' and iax != 2:
            bool3d = False
        elif layout == 'unfolded 3d' and iax == 2:
            bool3d = True
        if axlabels is not None:
            ax.set_xlabel(axlabels[iax][0]) 
            ax.set_ylabel(axlabels[iax][1]) 
            if bool3d:
                # shift the label closer to the axis
                ax.set_zlabel(axlabels[iax][2],labelpad=-7)


def plot_arrows_transitions(ax,X,indices,weight=None):
    """ 
    Plot arrows of transitions in data matrix.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    X : np.array
        Data array, any representation wished (X, psi, phi, etc).
    indices : array_like
        Indices storing the transitions.
    """
    step = 1
    width = plot_axis_to_data(ax,0.001)
    if X.shape[0] > 300:
        step = 5
        width = plot_axis_to_data(ax,0.0005)
    if X.shape[0] > 500:
        step = 30
        width = plot_axis_to_data(ax,0.0001)
    head_width = 10*width
    for ix,x in enumerate(X):
        if ix%step == 0:
            X_step = X[indices[ix]] - x
            # don't plot arrow of length 0 
            for itrans in range(X_step.shape[0]):
                alphai = 1
                widthi = width
                head_widthi = head_width
                if weight is not None:
                    alphai *= weight[ix,itrans]
                    widthi *= weight[ix,itrans]
#                     head_widthi *= weight[ix,itrans]
                if np.any(X_step[itrans,:1]):
                    ax.arrow(x[0], x[1],
                             X_step[itrans,0], X_step[itrans,1],
                             length_includes_head=True,
                             width=widthi, 
                             head_width=head_widthi,
                             alpha=alphai,
                             color='grey')

def plot_get_ax_size(ax,fig):
    """ 
    Get axis size

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    fig : matplotlib.Figure
        Figure.
    """
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    print(width,height)
    width *= fig.dpi
    height *= fig.dpi
    print(width,height)

def plot_axis_to_data(ax,width):
    """ 
    For a width in axis coordinates, return the corresponding in data
    coordinates.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    width : float
        Width in xaxis coordinates.
    """
    xlim = ax.get_xlim()
    widthx = width*(xlim[1] - xlim[0])
    ylim = ax.get_ylim()
    widthy = width*(ylim[1] - ylim[0])
    return 0.5*(widthx + widthy)

def plot_axis_to_data_points(ax,points_axis):
    """ 
    Map points in axis coordinates to data coordinates.

    Uses matplotlib.transform.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    points_axis : np.array
        Points in axis coordinates.
    """
    axis_to_data = ax.transAxes + ax.transData.inverted()
    return axis_to_data.transform(points_axis)

def plot_data_to_axis_points(ax,points_data):
    """ 
    Map points in data coordinates to axis coordinates.

    Uses matplotlib.transform.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    points_axis : np.array
        Points in data coordinates.
    """
    data_to_axis = axis_to_data.inverted()
    return data_to_axis(points_data)

def plot_timeseries(*args,**kwargs):
    """ 
    Plot X. See plot_timeseries_subplot.
    """
    pl.figure(figsize=(8,4),
              subplotpars=sppars(left=0.12,right=0.98,bottom=0.13))
    plot_timeseries_subplot(*args,**kwargs)

def plot_timeseries_subplot(X,varnames=[],highlightsX=[],
                            c = None,
                            xlabel='segments / pseudotime order',
                            ylabel='gene expression',
                            yticks=None,
                            xlim=None):
    """ 
    Plot X.
    """
    for i in range(X.shape[1]):
        pl.scatter(np.arange(X.shape[0]),X[:,i],
                   marker = '.', edgecolor = 'face',
                   c=mpl_plot.cl[i] if c is None else c,
                   label=(varnames[i] if len(varnames) > 0 else ''))
    ylim = pl.ylim()
    for ih,h in enumerate(highlightsX):
        pl.plot([h,h],[ylim[0],ylim[1]],
                '--',color='black')
    pl.ylim(ylim)
    if xlim is not None:
        pl.xlim(xlim)
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)
    if yticks is not None:
        pl.yticks(yticks)
    if len(varnames) > 0:
        pl.legend(frameon=False)

def plot_timeseries_as_heatmap(X,varnames=[]):
    """ 
    Plot timeseries as heatmap.

    Parameters
    ----------
    X : np.ndarray
        Data array.
    varnames : array_like
        Array of strings naming variables stored in columns of X.
    """
    if len(varnames) == 0:
        varnames = np.arange(X.shape[1])
    if varnames.ndim == 2:
        varnames = varnames[:,0]
    fig = pl.figure(figsize=(6,8))
    pl.imshow(np.array(X.T,dtype=np.float_),aspect='auto',
              interpolation='nearest')
    pl.colorbar(shrink=0.5)
    pl.yticks(range(X.shape[1]),varnames)

def check_datafile(filename):
    """
    Check whether the file exceeds 500 Bytes, otherwise download.

    If the file is very small, it likely is just a placeholder for a link to
    github. Download from github in this case.
    """
    # if file is smaller than 500 Bytes = 0.5 KB
    threshold = 500
    if os.path.getsize(filename) < threshold:
        # download the file
        # note that this has 'raw' in the address
        github_baseurl = r'https://github.com/falexwolf/scanpy/raw/master/'
        fileurl = github_baseurl + filename
#         log.m(0,'size of file',filename,'is below',threshold/500,' kilobytes') 
#         log.m(0,'    presumably it is a placeholder for a git-lfs file')
#         log.m(0,'    append _small to the filename')
#         log.m(0,'    and download large file from github using url')
#         log.m(0,'    ' + fileurl)
#         log.m(0,'    if you installed git-lfs, you can use \'git lfs checkout\'')
#         log.m(0,'    to update all files with their actual content') 
        log.m(0,'downloading data file from github using url')
        log.m(0,fileurl)
        log.m(0,'...this may take a while but only happens once')
        # make a backup of the small file
        from shutil import move
        basename, ext = os.path.splitext(filename)
        move(filename,basename+'_small'+ext)
        # download the file
        if sys.version_info >= (3,0): # Python 3
            from urllib.request import urlretrieve
        else: # Python 2
            from urllib import urlretrieve
        try:
            urlretrieve(fileurl,filename)
        except Exception as e:
            log.m(0,e)
            log.m(0,'when calling urlretrieve() in module scanpy.utils.utils')
            log.m(0,'--> is the github repo/url private?')
            log.m(0,'--> if you have access, manually download the file')
            log.m(0,fileurl)
            log.m(0,'replace the small file',filename,'with the downloaded file')
            quit()

def loadtxt_header(filename):
    """ 
    Return data as array and the header as string.
    """
    check_datafile(filename)
    header = ''
    for line in open(filename):
        if line.startswith("#"):
            header += line
        else:
            break
    return np.loadtxt(filename), header

def read_formatted_text(filename):
    """ 
    Read Excel or hdf5 file and return pandas data frame.

    Might have problems with hdf5 files that are not created with pandas.

    Reading from hdf5 is preferred, as it's faster. Creates an
    hdf5 file if it's not present yet.
    """

    X, header = loadtxt_header(filename)
    genenames = np.arange(X.shape[1]).astype(str)
    if len(header) > 0:
        potentialnames = header.split('\n')[-2].strip('#').split()
        if len(potentialnames) == X.shape[1]:
            genenames = np.array(potentialnames)

    # skip first column, it stores a time label
    X = X[:,1:]
    genenames = genenames[1:]
    # the root cell is the first in the datafile
    # this is the convention used for simulations
    xroot = X[0]
    labels = np.arange(X.shape[0]).astype(str)

    ddata = {
        'X' : X, 'xiroot' : xroot, 'rownames' : labels, 'colnames' : genenames
    }

    return ddata

def read_excel_hdf5_to_pandas(filename,sheet=''):
    """ 
    Read Excel or hdf5 file and return pandas data frame.

    Might have problems with hdf5 files that are not created with pandas.

    Reading from hdf5 is preferred, as it's faster. Creates an
    hdf5 file if it's not present yet.

    Parameters
    ----------
    filename : str
        Filename of data file.
    sheet : str
        Name of sheet in Excel or hdf5 file.

    Returns
    -------
    df : pd.DataFrame
        Pandas data frame.
    """
    filename_hdf5 = filename.replace('.xlsx','.h5')
    if os.path.exists(filename_hdf5):
        check_datafile(filename_hdf5)
    # if there is no hdf5 file, read from excel file
    # and write the hdf5 file
    if not os.path.exists(filename_hdf5):
        log.m(0,'reading file',filename)
        check_datafile(filename)
        try:
            df = pd.read_excel(filename,sheet)        
        except Exception as e:
            # in case this raises an error using Python 2 install xlrd via
            # sudo pip install xlrd
            print('try installing xlrd using "sudo pip install xlrd"') 
            raise e
        # TODO!!
        # write as hdf5 file, in Python 3 using numpy the following works
        df.to_hdf(filename_hdf5,sheet)
        # but for Python 2, we should use
        # with h5py.File(filename_hdf5, 'w') as f:
        #     f.create_dataset(filename_hdf5,data=df.values)
        # with a work around for genenames and labels
    # read from hdf5 file, because it's faster
    # pandas has restrictions here, so you might use read_hdf5 
    # if reading from a hdf5 file that was not created using pandas
    # see, for example, the bug report 
    # http://stackoverflow.com/questions/33641246/pandas-cant-read-hdf5-file-created-with-h5py
    else:
        log.m(0,'reading file',filename_hdf5)
        df = pd.read_hdf(filename_hdf5,sheet)
    return df

def read_hdf5(filename,key=''):
    """ 
    Read dataset from hdf5 file.

    Parameters
    ----------
    filename : str
        Filename of data file.
    key : str
        Name of dataset in the file.

    Returns
    -------
    d : dict
        Dictionary consisting of np.array, rownames and colnames.
    """
    # log.m(0,'reading file',filename)
    check_datafile(filename)
    f = h5py.File(filename, 'r')
    # the following is necessary in Python 3, because only
    # a view and not a list is returned
    avail_keys = [k for k in f.keys()]
    if key != '':
        # the slicing operation in the following line transforms
        # the hdf5 object into a np.array
        # alternatively one might use "()"
        # that's very strange, but somehow a transposition is made
        # when transforming to numpy
        X = f[key][:]
        if X.dtype.kind == 'S':
            X = X.astype(str)
        rownames = np.array([])
        colnames = np.array([])
        if key+'_rownames' in avail_keys:
            colnames = f[key+'_rownames'][:].astype(str)
        if key+'_colnames' in avail_keys:
            rownames = f[key+'_colnames'][:].astype(str)
    else:
        print('The file',filename,'stores the following keys:\n',avail_keys,
              '\n Call the function with one of them.')
    return {'X' : X, 'rownames' : rownames, 'colnames' : colnames}


def read_file_to_dict(filename,ext='h5'):
    """ 
    Read Excel file and return dict with sheet names as keys.

    Parameters
    ----------
    filename : str
        Filename of data file.

    Returns
    -------
    result : dict
        Returns OrderedDict.
    """
    d = odict([])
    if ext == 'h5':
        with h5py.File(filename, 'r') as f:
            for key in f.keys():
                # the '()' means 'read everything'
                # alternatively ':' works if one does not read a scalar type
                value = f[key][()]
                # print(key)
                # print(type(value),value.dtype)
                if value.dtype.kind == 'S':
                    d[key] = value.astype(str)
                else:
                    d[key] = value
    elif ext == 'xlsx':
        xl = pd.ExcelFile(filename)
        for sheet in xl.sheet_names:
            d[sheet] = xl.parse(sheet).values
    return d

def write_dict_to_file(filename,d,ext='h5'):
    """ 
    Write content of dictionary to Excel file.

    Parameters
    ----------
    filename : str
        Filename of data file.
    d : dict
        Dictionary storing array_like data.
    ext : string
        Determines file type, allowed are 'h5' (hdf5),
        'xlsx' (Excel) or 'csv' (comma separated value file).
    """
    directory = os.path.dirname(filename)
    if not os.path.exists(directory):
        os.makedirs(directory)
    if ext == 'h5':
        with h5py.File(filename, 'w') as f:
            for key, value in d.items():
                if type(value) != np.ndarray:
                    value = np.array(value)
#                 print(key)
#                 print(type(value),value.dtype,value.dtype.kind,value.shape)
                if value.dtype.kind == 'U':
                    f.create_dataset(key,data=value.astype(np.string_))
                else:
                    try:
                        f.create_dataset(key,data=value)
                    except Exception as e:
                        log.m(0,'error creating dataset for key =', key)
                        raise e
        # using pandas, should be slower than using h5py
        # but would be fully compatible for reading with pandas
        # with pd.HDFStore('store.h5') as store:
        #     for key, value in d.items():
        #         store[key] = pd.DataFrame(value)
    elif ext == 'xlsx':
        with pd.ExcelWriter(filename,engine='openpyxl') as writer:
            for key, value in d.items():
                pd.DataFrame(value).to_excel(writer,key)

def read_parameters(filename):
    """ 
    Read parameter dictionary from data file.

    Assumes that parameters are specified in the format:  
    # key = value | type
    where 'type' needs to be a valid python type.
    
    Parameters
    ----------
    filename : str
        Filename of data file.

    Returns
    -------
    pars : dict
        Dictionary that stores parameters.
    """
    pars = {}
    for line in open(filename):
        if line.startswith('#') and '=' in line:
            line = line[1:]
            keyval, type = line.split('|')[:2]
            key, val = keyval.split('=')
            key = key.strip(); val = val.strip(); type = type.strip('\n\t ')
            if type != 'str':
                val = eval(type+'('+val+')')
            pars[key] = val
    return pars

def subsample(X,subsample=1,seed=0):
    """ 
    Subsample a fraction of 1/subsample samples from the rows of X.

    Parameters
    ----------
    X : np.ndarray
        Data array.
    subsample : int
        1/subsample is the fraction of data sampled, n = X.shape[0]/subsample.
    seed : int
        Seed for sampling.

    Returns
    -------
    Xsampled : np.ndarray
        Subsampled X.
    rows : np.ndarray 
        Indices of rows that are stored in Xsampled.
    """
    if subsample == 1:
        return X, np.arange(X.shape[0],dtype=int)
    n = int(X.shape[0]/subsample)
    log.m(0,'subsampling to',n,'of',X.shape[0],'data points')
    if seed == 0:
        # this sequence is defined simply by skipping rows
        # is faster than sampling
        rows = np.arange(0,X.shape[0],subsample,dtype=int)
        Xsampled = np.array(X[rows])
    if seed > 0:
        np.random.seed(seed)
        Xsampled, rows = subsample_n(X,n=X.shape[0]/subsample)
    return Xsampled, rows

def subsample_n(X,n=0,seed=0):
    """ 
    Subsample n samples from rows of array.

    Parameters
    ----------
    X : np.ndarray
        Data array.
    seed : int
        Seed for sampling.

    Returns
    -------
    Xsampled : np.ndarray
        Subsampled X.
    rows : np.ndarray 
        Indices of rows that are stored in Xsampled.
    """
    if n < 0:
        raise ValueError('n must be greater 0')
    np.random.seed(seed)
    n = X.shape[0] if (n == 0 or n > X.shape[0]) else n
    rows = np.random.choice(X.shape[0],size=n,replace=False)
    Xsampled = np.array(X[rows])
    return Xsampled, rows

def comp_distance(X,metric='euclidean'):
    """ 
    Compute distance matrix for data array X
    
    Parameters
    ----------
    X : np.ndarray
        Data array (rows store samples, columns store variables).
    metric : string
        For example 'euclidean', 'sqeuclidean', see sp.spatial.distance.pdist.

    Returns
    -------
    D : np.ndarray
        Distance matrix.
    """
    D = sp.spatial.distance.pdist(X,metric=metric)
    D = sp.spatial.distance.squareform(D)
    log.mt(0,'computed euclidian distance matrix')
    log.m(5,D)
    if False:
        pl.matshow(D)
        pl.colorbar()
    return D

def hierarch_cluster(M):
    """ 
    Cluster matrix using hierarchical clustering.

    Parameters
    ----------
    M : np.ndarray
        Matrix, for example, distance matrix.

    Returns
    -------
    Mclus : np.ndarray 
        Clustered matrix.
    indices : np.ndarray
        Indices used to cluster the matrix.
    """
    link = sp.cluster.hierarchy.linkage(M)
    indices = sp.cluster.hierarchy.leaves_list(link)
    Mclus = np.array(M[:,indices])
    Mclus = Mclus[indices,:]
    log.mt(0,'clustered matrix')
    if False:
        pl.matshow(Mclus)
        pl.colorbar()
    return Mclus, indices

def plot_neighborhood(gp,g,i,X_hist,nbrs):
    fig = pl.figure(figsize=(4,4))
    ax = fig.add_subplot(111,projection='3d')
    ax.set_title('x-y: G'+str(gp)+' $\stackrel{?}{\leftarrow}$ z: G'+str(g))
    ax.scatter(X_hist[:,gp,0],X_hist[:,gp,1],X_hist[:,g,0],
               c=np.arange(len(X_hist)),edgecolors='face',s=2,alpha=0.5)
    ax.scatter(X_hist[nbrs[2][i],gp,0],X_hist[nbrs[2][i],gp,1],X_hist[nbrs[2][i],g,0],
               marker='x',c='black',edgecolors='face',s=10,alpha=1)
    pl.show(block=False)

def functest(X,X_hist,gs,sol_chains=[],indices=[],indices_and=[],verbosity=0,
             g_control=-1,drift=False):
    """ Statistical test for functional dependency between time series.
    """
    for gp in gs:
        for g in gs:
            if gp != g:
                indices_and[gp,g] = indices[:,gp] * indices[:,g]
    #
    trange = np.arange(X.shape[0])
    X_domain = np.array(X_hist)
    # account for g_control variable, combinatorial regulation
    if g_control != -1:
        X_dummy = np.array([X_hist[:,g_control,0] for g in range(X.shape[1])]).T
        X_domain = np.array(np.concatenate((X_hist,X_dummy[:,:,np.newaxis]),axis=2))
    if drift:
        # most naive drift computation ever
        # this assumes that the drift is weighted with a factor 10
        X_drift = 10*(X_hist[:,:,1] - X_hist[:,:,0])
        X_domain = np.array(np.concatenate((X_hist,X_drift[:,:,np.newaxis]),axis=2))
    # the neighbors array stores for each time point and each variable
    # how many neighbors fall into a ball of the local rate of change
    neighbors_number_dscale = np.ones((X_domain.shape[0],X_domain.shape[1],2))
    neighbors_dscale_cod = np.ones((X_domain.shape[0],X_domain.shape[1],1))
    neighbors_indices = [[[] for j in range(X_domain.shape[1])] 
                              for i in range(X_domain.shape[0])]
    neighbors_dist_tree = [ None for j in range(X.shape[1])]
    # one-dimensional version
    neighbors_dist_matrix_cod = [ None for j in range(X.shape[1])]
    # should depend on the distance norm chosen!
    noise = 1e-2
    # only consider those genes that have some non-trivial data
    igs_delete = []
    for ig,g in enumerate(gs):
        if not np.any(indices[:,g]):
            if verbosity > 1:
                print('no dynamic information for variable '+str(g))
            igs_delete.append(ig)
    if len(igs_delete) > 0:
        gs = np.delete(gs,np.array(igs_delete))
    # compute distances for each variable
#     print(len(X_domain))
#     print(X_domain[indices[:,1],1][:6])
#     print(X_domain[indices[:,1],1][-5:])
#     print(X[indices[:,1],0][:8])
#     print(sp.spatial.distance.squareform(
#                                   sp.spatial.distance.pdist(X_domain[indices[:,1],1],
#                                                              metric='euclidean')))
#     print(len(trange[indices[:,1]]))
#     print(trange[indices[:,1]][:8])
    for g in gs:
        # construct a tree of distances
        neighbors_dist_tree[g] = scipy.spatial.KDTree(X_domain[indices[:,g],g])
        # take all indices here, that is, use the full trange
        neighbors_dist_matrix_cod[g] = sp.spatial.distance.squareform(
                                        sp.spatial.distance.pdist(X[:,g,np.newaxis],
                                                               metric='euclidean'))
        for t in trange[indices[:,g]]:
            d = sp.spatial.distance.euclidean(X_domain[t,g],X_domain[t-1,g])
            d += sp.spatial.distance.euclidean(X_domain[t+1,g],X_domain[t,g])
            d_cod = sp.spatial.distance.euclidean(X[t,g],X[t-1,g])
            d_cod += sp.spatial.distance.euclidean(X[t+1,g],X[t,g])
            d *= 0.5
            d_cod *= 0.5
            # for the "convergence" check, we need the noise assumption,
            # i.e., checking whether a variable converges earlier than another
            # for the "monotonicity" check, it is not necessary
            # it's rather harmful as we are then sensitive to converging solutions,
            # although the latter ok
#             if not diag:
#                 d += noise
#                 d_cod += noise
            # search all points that are equivalent to this point X_domain[t,g]
            # within the local scale of variation d
#            indices_nl = neighbors_dist_tree[g].query_ball_point(X_domain[t,g],d)
            indices_n = neighbors_dist_tree[g].query(X_domain[t,g],3)[1]; indices_nl = indices_n.tolist()
            # we don't remove the t index itself
            # remove t index itself, here is the only difference to the
            # class CCM, this is done in a quite strange way, here
#            print(g,t,trange[indices[:,g]][indices_nl])
            # remove the time t from the neighbor list
#             try:
#                 it = trange[indices[:,g]][indices_nl].tolist().index(t)
#                 indices_nl.pop(it)
#             except ValueError:
#                 print(g,t,'t itself not neighbor list')
#             print(g,t,trange[indices[:,g]][indices_nl])
            indices_n = np.array(indices_nl)
#             if g == 1 and t <= 9:
#                 print(indices_n)
#             if g == 0 and t <= 9:
#                 print(t,X[t,g],d_cod)
            neighbors_indices[t][g] = indices_n
            # store how many neighbors there are
            neighbors_number_dscale[t,g,0] = len(indices_n)
            # store the local scale of variation
            neighbors_number_dscale[t,g,1] = d
            neighbors_dscale_cod[t,g,0] = d_cod

#     print(neighbors_dscale_cod[indices[:,1],0,0][:8])
    #
    var_codomain = np.zeros((X_domain.shape[1],X_domain.shape[1],X_domain.shape[0]))
    var_dynamic = np.zeros((X_domain.shape[1],X_domain.shape[1],X_domain.shape[0]))
    pred_ccm = np.zeros((X_domain.shape[1],X_domain.shape[1],X_domain.shape[0]))
    # We test whether a functional relation from gp to g exists.
    # The returned test statistic is large if a functional relation is unlikely,
    # using as indicator a high variance in the codomain relative to the dynamic
    # variance in the codomain.
    for gp in gs:
        for g in gs:
            if g != gp:
                for it,t in enumerate(trange[indices[:,gp]]):
                    if indices_and[gp,g,t]:
                        neighbors_indices_and = []
#                         print(neighbors_indices[t][gp])
                        for itp,tp in enumerate(trange[indices[:,gp]][neighbors_indices[t][gp]]):
                            if tp in trange[indices_and[gp,g]]:
                                neighbors_indices_and.append(neighbors_indices[t][gp][itp])
                        neighbors_indices_and = np.array(neighbors_indices_and)
#                         print(neighbors_indices_and)
                        if len(neighbors_indices_and) == 0:
                            var_codomain_t = 0
                            pred_ccm[gp,g,t] = 0
                        else:
                            pred_ccm[gp,g,t] = np.average(X[trange[indices[:,gp]][neighbors_indices_and],g])
                            # compute var_codomain_t as average distance between all neighbors
                            var_codomain_t = 0
                            drift_codomain = 0
                            for tp in trange[indices[:,gp]][neighbors_indices_and]:
                                for tpp in trange[indices[:,gp]][neighbors_indices_and]:
                                    if tp != tpp:
                                        var_codomain_t += neighbors_dist_matrix_cod[g][tp,tpp]
                            if len(neighbors_indices_and) > 1:
                                var_codomain_t = var_codomain_t/len(neighbors_indices_and)/(len(neighbors_indices_and)-1)
                            # compute var_codomain_t as average distance to current point
                            var_codomain_t = np.average(
                                  neighbors_dist_matrix_cod[g][t,trange[indices[:,gp]][neighbors_indices_and]])
                        var_codomain[gp,g,t] = var_codomain_t
                        var_dynamic[gp,g,t] = neighbors_dscale_cod[t,g,0]
#                         print(X_domain[neighbors_indices[t][gp],gp])
#                         print(gp,g,t,
#                               neighbors_dist_matrix[gp,t,neighbors_indices[t][gp]],
#                               neighbors_indices[t][gp])
#                        print(X[neighbors_indices[t][gp],g])

#                        print(gp,g,t,pred_ccm[gp,g,t])

#                 if gp == 1:
#                     print(len(var_dynamic[gp,g,indices[:,gp]]))
#                     print(var_dynamic[gp,g,indices[:,gp]][:5])
#                     print(var_codomain[gp,g,indices[:,gp]][:5])
#                     print(np.average(var_dynamic[gp,g,indices[:,gp]]) - np.average(var_codomain[gp,g,indices[:,gp]]))
#                     print(np.average(var_dynamic[gp,g,indices_and[gp,g]]) - np.average(var_codomain[gp,g,indices_and[gp,g]]))

    return var_codomain, var_dynamic, pred_ccm, neighbors_indices

def smooth_data_hist(X_reconst,sol_chains,k=5,dimE=2,future=True):
    """ Smooth data in intervals defined by sol_chains, return X_hist.
        
        future: allows to switch between derivative and future values
        k : how much neighboring points for running average
        dimE : embedding dimension
    """
    sol_reconst = np.zeros(len(sol_chains))
    for i in range(len(sol_chains)):
        sol_reconst[i] = sol_chains[i][-1]
    X_reconst_smooth = np.zeros(X_reconst.shape)
    X_reconst_slope = np.zeros(X_reconst.shape)
    for ig,Xg in enumerate(X_reconst.T):
        sol_old = 0
        for isol,sol in enumerate(sol_reconst):
            if future:
                Xg_sol_smooth, Xg_sol_slope = run_ave_future(Xg[sol_old:sol],k)
            else:
                Xg_sol_smooth, Xg_sol_slope = run_ave_slope(Xg[sol_old:sol],k)
            X_reconst_smooth[sol_old:sol,ig] = Xg_sol_smooth
            X_reconst_slope[sol_old:sol,ig] = Xg_sol_slope
            sol_old = sol

    # Build lagged vectors / shadow manifold.
    # from point values and estimate of slope
    X_hist = np.zeros((X_reconst.shape[0],X_reconst.shape[1],dimE))
    # zero dimension
    X_hist[:,:,0] = X_reconst_smooth
    # first dimension: estimate derivative at each point
    if dimE == 2:
        X_hist[:,:,1] = X_reconst_slope
    #
    return X_hist

def check_convergence(y,sol_chains=[],k=5,p_thres=0.01,noise=None):
    """ use running average to determine convergence
    """
    y_ave = np.zeros(y.shape)
    y_change = np.zeros(y.shape)
    y_std = np.zeros(y.shape)
    indices = np.zeros(y.shape,dtype=bool)
    N = 0
    for isol,sol in enumerate(sol_chains):
        y_ave[sol[0]:sol[1]], y_change[sol[0]:sol[1]], y_std[sol[0]:sol[1]] = run_ave_change_std(y[sol[0]:sol[1]],k)
        noise_estim = y_std[sol[1]-1] if noise is None else noise
        z = np.abs(y_change[sol[0]:sol[1]]/(noise_estim/np.sqrt(k)))
        p = 1-sp.stats.norm.cdf(z)
#         print(y_std[sol[1]-1]/np.sqrt(k))
#         print(y_change[sol[0]:sol[1]])
#         print(z)
#         print(p)
#         l = sol[1]-sol[0]
#         sol1 = (l-next((i for i,pi in enumerate(p[::-1]) if pi < p_thres),l))
#         for i in range(sol1):
#             indices[sol[0]+i] = True
#         N += sol1
        indices[sol[0]:sol[1]] = p < p_thres
#        print(isol,y[sol[0]:sol[0]+20])
#        print(indices[sol[0]:sol[0]+20])
    # only subsequent values indicate progression, otherwise, this might just be an outlier
    # set it to False in this case
    for t in range(len(indices)-1,1,-1):
        if indices[t]:
            if not indices[t-1]:
                indices[t] = False
    # for each interval, append a few time steps to capture transient behavior
    for t in range(len(indices)-2*k,0,-1):
        if indices[t]:
            if not indices[t+1]:
                for delta in range(1,10*k):
                    if t+delta >= indices.shape[0]:
                        break
                    indices[t+delta] = True
    N = sum(indices)
    return y_ave, y_change, y_std, indices, N
 
def run_ave_change_std(y,k=5):
    """ Running average, change rate, standard deviation for 1d array.
    """
    y_ave = np.zeros(y.shape)
    y_change = np.zeros(y.shape)
    y_std = np.zeros(y.shape)
    for i in range(y.shape[0]):
        i_start = max(0,i-k/2)
        i_end = min(y.shape[0]-1,i+k/2)
        y_ave[i] = np.average(y[i_start:i_end+1])
        y_std[i] = np.std(y[i_start:i_end+1])
    y_change[:-1] = y_ave[1:] - y_ave[:-1]
    return y_ave, y_change, y_std

def run_ave_std(y,k=5):
    """ Running average and standard deviation for 1d array.
    """
    y_ave = np.zeros(y.shape)
    y_std = np.zeros(y.shape)
    for i in range(y.shape[0]):
        i_start = max(0,i-k/2)
        i_end = min(y.shape[0],i+k/2)
        y_ave[i] = np.average(y[i_start:i_end+1])
        y_std[i] = np.std(y[i_start:i_end+1])
    return y_ave, y_std

def run_ave(y,k=5):
    """ Running average for 1d array.
    """
    y_ave = np.zeros(y.shape)
    for i in range(y.shape[0]):
        i_start = max(0,i-k/2)
        i_end = min(y.shape[0],i+k/2)
        y_ave[i] = np.average(y[i_start:i_end+1])
    return y_ave

def run_ave_future(y,k=5):
    """ Running average for 1d array. Returned are also the shifted values,
        shifted by one time step into the future. A time step is defined by k.
    """
    y_ave = np.zeros(y.shape)
    y_fut = np.zeros(y.shape)
    for i in range(y.shape[0]):
        i_start = max(0,i-k/2)
        i_end = min(y.shape[0]-1,i+k/2)
        y_ave[i] = np.average(y[i_start:i_end+1])
    #
    y_fut[:-k] = y_ave[k:] 
    y_fut[-k:] = y_ave[-1]*np.ones(k)
    return y_ave, y_fut

def run_ave_slope(y,k=5):
    """ Running average and slope for 1d array.
    """
    y_ave = np.zeros(y.shape)
    y_slope = np.zeros(y.shape)
    for i in range(y.shape[0]):
        i_start = max(0,i-k/2)
        i_end = min(y.shape[0]-1,i+k/2)
        y_ave[i] = np.average(y[i_start:i_end+1])
    # here we use an even higher value for k
    for i in range(y.shape[0]):
        i_start = max(0,i-k/2)
        i_end = min(y.shape[0]-1,i+k/2)
        y_slope[i] = 10*(y_ave[i_end] - y_ave[i_start])/(i_end-i_start)
    return y_ave, y_slope

def cheb_fit_slope(y,k=5):
    """ Fit Chebyshev polynomials and return slope.
    """
    x = np.linspace(-1,1,len(y))
    coef = np.polynomial.chebyshev.chebfit(x, y, k)
    y_cheb = np.polynomial.chebyshev.chebval(x, coef)
    #
    y_slope = np.array(y_cheb[1:] - y_cheb[:-1])
    y_slope = 10*np.r_[y_slope[:1],y_slope]
    return y_cheb, y_slope

def corr_partial(df):
    """ Compute partial correlation for data frame
        parameters:
          df: data frame
        returns:
          Corr: numpy array of partial correlations
    """
    Corr = df.corr(method='pearson')
    Corr = np.array(np.linalg.inv(Corr))
    Diag = np.array(np.diag(Corr))
    Corr = -Corr
    for i in range(Corr.shape[0]):
         for j in range(Corr.shape[1]):
               Corr[i,j] /= np.sqrt(Diag[i]*Diag[j])
    return Corr

def mutual_info(x, y, bins=10):
    """ compute mutual information for two 1d arrays
    """
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = sklearn.metrics.mutual_info_score(None, None, contingency=c_xy)
    return mi

def mutual_info_A(A, bins=10):
    """ compute mutual information of columns of an array A 
    """
    matMI = np.zeros((A.shape[1],A.shape[1]))
    for ix in np.arange(A.shape[1]):
        for jx in np.arange(ix,A.shape[1]):
            matMI[ix,jx] = mutual_info(A[:,ix], A[:,jx], bins)
            matMI[jx,ix] = matMI[ix,jx]    
    return matMI

class CCM:
    """ We wonder whether gene g regulates gene g'. The only input data that we use
        to answer is, is an array of gene expression values with time labels.
    """

    def __init__(self,X,verbosity=0,
                 spacings=[200, 100, 50, 25, 12],
                 dimE = 0):
        """ Parameters
            X : data array, sorted according to increasing time indices
            verbosity : integer 
            spacings : list of tau values, lags in time series
            dimE : dimension of Euclidian state space, if 0 defaults to the number
                   of columns in X
        """
        self.verbosity = verbosity
        self.X = np.array(X)
        self.dimE = self.X.shape[1] if dimE == 0 else dimE
        self.dimE1 = self.dimE-1
        self.indexrange = list(range(self.X.shape[0]))
        self.trangeorig = self.indexrange[self.dimE1:]
        # spacings, that is \tau values 
        self.spacerange = spacings
        # futuretimes stores the number of time steps that we predict
        # into the future
        self.futuretimes = [0]
        # this is going to be compute at runtime
        self.trange = [] 
        self.trangelengths = []
        # output
        if self.verbosity > 0:
            print('... length data array',self.X.shape[0])
            print('... spacings considered:',self.spacerange)
        
    def rho_TwoColsFuture(self,g,gp):
        """ Computes the correlation of the predictor for
            g given the data for gp, and the actual data for g.

            This is done for a sequence of spacings (i.e. sequence 
            of time series lengths) given in the constructor.
            The prediction is made for future time values.

            Parameters
            g : column index of predicted variable (predictee)
            gp : column index of the predictor

            Returns
            rho : an array, first dim: correlation values with increasing time series
                  length, second dim: number of future time steps
        """
        if self.verbosity > 0: 
            print('... compute correlation for g and gp:',g,gp)
        # return an array that stores correlation of predictor with
        # and actually measured values for different numbers of 
        # future time steps 
        rho = np.zeros((len(self.spacerange),len(self.futuretimes)))
        # inits these class variables to have zero lengths
        self.trange = []
        self.trangelengths = []
        # loop over spacings
        if self.verbosity > 1: 
            print('L','futuretime','rho')
        # loop over different future time steps (fts) values
        for ispace,space in enumerate(self.spacerange):
            # we consider a coarsened grid of time points
            self.trange = self.indexrange[::space]
            # we have to cut of the first few time steps for the 
            # the predicition times, as we form lagged coordinate vectors
            offset = space*self.dimE1
            offsetfuture = (len(self.trange)-max(self.futuretimes))*space
            self.trangelag = np.array(self.indexrange[offset:offsetfuture:space])
            self.trangelengths.append(len(self.trangelag))
            # compute correlation coefficient
            if self.verbosity > 2:
                print('reference times')
                print(self.trangelag)
            # if we don't have enough data, exit
            if len(self.trange) < 7:
                print('... choose smaller spacing, not enough data')
                quit()
            # compute prediction
            if self.verbosity == 2:
                print('t+fts*space |','prediction |','actual value |')
            g_hat_list = np.array([ self.compute_g_hat_t(g,gp,t,space,self.futuretimes) 
                                    for t in self.trangelag])
            for ifts,fts in enumerate(self.futuretimes):
                r = np.corrcoef(self.X[self.trangelag+fts*space,g],
                                g_hat_list[:,ifts],rowvar=0)
                rho[ispace,ifts] = r[0,1]
            if self.verbosity > 1: 
                print('correlation')
                for ifts,fts in enumerate(self.futuretimes):
                    print(self.trangelengths[-1],fts,rho[-len(self.futuretimes)+ifts])
        return rho

    def compute_g_hat_t(self,g,gp,t,space,futuretimes):
        """ Compute the predicted value of g for a single time point t+fts*space
            as average over nearest neighbors on the shadow manifold M_gp.

            Here, fts describes the number of time steps into the future. The 
            time step length is measured by space, a space for the time index. 
            fts is expected to be an array or a list.

            We have to find the dimE+1 nearest neighbors in the shadow 
            manifold M_gp.

            Returns prediction (array of the same dimension as fts).
        """
        # compute the array of distances in the shadow manifold 
        # with respect to the current time point t
        dists = [] 
        for tj in self.trangelag:
            dists.append(self.dist(gp,t,tj,space))
        dists = np.array(dists)
        # find the dimE+1 smallest distances in the shadow manifold
        # using argpartition
        # as the following contains the case t==tj where the distance
        # is zero, and which we will throw away, take dimE+2
        # these are the indices for the dimE+1 nearest neighbors
        # the indices have to be used within self.trangelag
        id_tn = np.argpartition(dists,range(self.dimE+2))[:self.dimE+2]
        # get the distance and time point for the nearest neighbor
        d1 = dists[id_tn[1]]
        t1 = self.trangelag[id_tn[1]]
        # control output
        if self.verbosity > 2:
            print('t','time points')
            print(t,np.array(self.trangelag)[id_tn])
            print('t','distances')
            print(t,dists[id_tn])
        # if the nearest neighbor is very close, it get's a high weight
        if d1 < 1e-12: 
            d1 = 1e-6
        # Loop over dimE+1 nearest neighbors.
        # To make prediction for the same time step or fts time steps
        # in the future.
        # Don't consider the case t==tj: start at index 1.
        den = 0 # denominator for normalisation
        prediction = np.zeros(len(futuretimes))
        offset = space*self.dimE1
        #
#         t = 1; tj = 2
#         xlag = self.X[t-offset:t+1:space,gp]
#         xlagj = self.X[tj-offset:tj+1:space,gp]
#         print(xlagj,self.X[tj-offset:tj+1:space,gp])
#         print(gp,g,t,dists[id_tn],id_tn)
#        print(self.X[self.trangelag[id_tn[1:]],g])
#         for idx in id_tn[1:]:
#             tj = self.trangelag[idx]
#             dj = dists[idx]
#             w = np.exp(-dj/d1)
#             # consider how point evolves fts timesteps in the future
#             for i,fts in enumerate(futuretimes):
#                 prediction[i] += self.X[tj+fts*space,g]*w
#             den += w
#         prediction /= den
        prediction = [np.average(self.X[self.trangelag[id_tn[1:]],g])]
#         if t==2 and gp == 2:
#             print(self.trangelag[id_tn])
#             quit()
        #
        np.set_printoptions(precision=2)
#         if self.verbosity > 2:
#             print('t |','fts |','prediction |','actual value |')
#             # the 0th neighbor is the point itself if it is the same manifold
#         if self.verbosity > 1:
#             for i,fts in enumerate(futuretimes):
#                 print(t,fts,
#                       prediction[i],
#                       self.X[t+fts*space,g])
#        print(gp,g,t,prediction[0])
        return prediction
                    
    def dist(self,gp,t,tj,space):
        """ Compute distance of lagged coordinate vectors in the shadow 
            manifold of variable gp.
        
            Parameters
            gp : index of gene
            t, tj : time indices 
            space : index spacing for the current time step considered
        """
        offset = space*self.dimE1
        # lagged-coordinate vector xlag for t
        xlag = self.X[t-offset:t+1:space,gp]
        # lagged-coordinate vector xlag for tj
        xlagj = self.X[tj-offset:tj+1:space,gp]
        return np.linalg.norm(xlag-xlagj)

