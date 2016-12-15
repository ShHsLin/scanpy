"""

Diffusion Maps and Diffusion Pseudotime

From package scanpy (https://github.com/falexwolf/scanpy).

Implements Diffusion Maps and Diffusion Pseudotime for analysis of single-cell
data.

Written in Python 3 (compatible with 2) using Numpy, Scipy, Matplotlib, Pandas
(use Anaconda package for all modules).

Copyright (c) 2016 F. Alexander Wolf (http://falexwolf.de).

References
----------
- Diffusion Maps: Coifman et al., PNAS 102, 7426 (2005).
- Diffusion Pseudotime: Haghverdi et al., Nature Methods 13, 3971 (2016).

See also
--------
- Diffusion Maps applied to single-cell data: Haghverdi et al., Bioinformatics
  31, 2989 (2015).
- Diffusion Maps as a flavour of spectral clustering: von Luxburg,
  arXiv:0711.0189 (2007).

"""   

# standard modules
import os
import warnings
import argparse
from string import ascii_lowercase
from collections import OrderedDict as odict
# scientific modules
import numpy as np
import scipy as sp
import matplotlib
with warnings.catch_warnings():
    # import of pyplot raises strange warning only in python 2
    warnings.simplefilter('ignore')
    import matplotlib.pyplot as pl
# scanpy modules
import scanpy 
import scanpy.preprocess
import scanpy.utils as utils
import scanpy.logfile as log

sppars = matplotlib.figure.SubplotParams

def dpt(X,xiroot,params={}):
    """
    Perform DPT analsysis as of Haghverdi et al. (2016).

    Parameters
    ----------
    X : np.ndarray
        Data array, rows store observations, columns covariates.
    xiroot : np.ndarray or int
        Root of stochastic process on data points (root cell), specified either
        as expression vector of shape X.shape[1] or as index.
    params : dict, optional
        Parameter dictionary for class DPT. See DPT.params attribute
        for description of possible values.

    Returns
    -------
    ddpt : dict
        Result dictionary.
    """

    diffmap, ddpt = _diffmap(X,params)

    # compute M matrix of cumulative transition probabilities,
    # see Haghverdi et al. (2016)
    diffmap.compM()
    # compute DPT distance matrix, called Ddiff
    diffmap.compDdiff()

    # set root point, if it's a gene expression value, first locate the root
    if type(xiroot) == np.ndarray:
        diffmap.find_root(xiroot)
    # if it's an index, directly set the index
    else:
        diffmap.iroot = xiroot
    ddpt['iroot'] = np.array([diffmap.iroot])

    # pseudotimes are distances from root point
    diffmap.set_pseudotimes()
    ddpt['pseudotimes'] = diffmap.pseudotimes

    # detect branchings and partition the data into corresponding segments
    diffmap.branchings_segments()
    ddpt['segs'] = diffmap.segs
    # the name convention for the dictionary aims at making memorizing simple
    # the one for the class member aims at consistence with other class members
    # and local variables used in the class
    ddpt['segtips'] = diffmap.segstips
    ddpt['seglabels'] = diffmap.segslabels
    # the ordering according to segments and pseudotimes
    ddpt['indices'] = diffmap.indices
    ddpt['changepoints'] = diffmap.changepoints

    # TODO: put the following in a place where it's more appropriate
    # this is the color map obtained when calling scatter plots with
    # seglabels as argument for color
    cmap = matplotlib.cm.get_cmap('jet')
    norm = matplotlib.colors.Normalize(vmin=0,vmax=ddpt['segs'].shape[0]-1)
    segcolors = np.array([cmap(norm(c)) for c in range(ddpt['segs'].shape[0])])
    ddpt['segcolors'] = segcolors
    
    return ddpt

def diffmap(X,params={}):
    """
    Compute diffusion map embedding as of Coifman et al. (2005).

    Also implements the modifications to diffusion map introduced by Haghverdi
    et al. (2016).

    Return dictionary that stores the new data representation 'rbasis', which 
    consists of the first few eigenvectors of a kernel matrix of the data, and
    the eigenvalues 'evalues'. 

    If params sets the parameter 'sym' to False, then there is also another
    representation of the data 'lbasis', which corresponds to the left
    eigenbasis of the asymmetric transition matrix.

    Parameters
    ----------
    X : np.ndarray
        Data array, rows store observations, columns covariates.
    params : dict
        Parameter dictionary for class DPT. See DPT.params attribute
        for description of possible values.

    Returns
    -------
    ddmap : dict
        Dictionary that contains 'rbasis' and 'evalues'. It might contain
        'lbasis'.
    """
    diffmap, ddmap = _diffmap(X,params)

    return ddmap

def _diffmap(X,params={}):
    """
    See function diffmap.

    Returns
    -------
    diffmap : DPT
        Instance of DPT object.
    ddmap : dict
        Dictionary as returned by function diffmap.
    """
    # initialization of DPT merely computes diffusion map
    diffmap = DPT(X,params=params)
    # write results to dictionary
    ddmap = odict([])
    ddmap['rbasis'] = diffmap.rbasis
    ddmap['evalues'] = diffmap.evalues
    if 'sym' in params and params['sym'] == False:
        ddmap['lbasis'] = diffmap.lbasis

    return diffmap, ddmap

class DPT:
    """ 
    Diffusion Pseudotime and Diffusion Map.

    Diffusion Pseudotime as of Haghverdi et al. (2016) and Diffusion Map as of
    Coifman et al. (2005).

    Also implements the modifications to diffusion map introduced by Haghverdi
    et al. (2016).
    """

    params = odict([
        ('num_branchings', 1),
        ('method', 'local'),
        ('k', 5),
        ('sigma', 0.2),
        ('fix_nr_neighbors', False),
        ('sym', True),
        ('allow_branching_at_root', False)
    ])
    """ 
    Dictionary storing parameters with default values.

    Is overwritten by __init__.

    Parameters
    ----------
    'num_branchings' : int
        Number of branchings to detect.
    'method' : {'local', 'standard'}
        Allows to choose the kernel for construction of the graph.
    'k' : int
        If method == 'local', you might provide a key 'k' that allows to
        specify the number of nearest neighbors considered to choose the
        kernel width.
    'sigma' : float
        If method == 'standard', you might provide a key 'sigma' that
        allows to choose the width of the Kernel Gaussian.
    'fix_nr_neighbors' : bool
        If set to True, use a hard threshold to restrict the number of
        neighbors to k. Otherwise, use a Gaussian Kernel to assign a low
        weight to distant neighbors.
    'sym' : bool
        Instead of the right eigenbasis of the asymmetric transition
        matrix, consider the basis of the conjugate symmetric kernel
        matrix that shares the same eigenvalues with the transition
        matrix. Note that only the asymmetric version reproduces the
        results by Coifman et al. (2005). Still, the default is the
        symmetric version.
    'allow_branching_at_root' : bool
        Allow to have branching directly at root point.
    """

    def __init__(self,X,params={}):
        """ 
        Initilization of DPT computes a diffusion map.

        The corresponding class members are described below and can be retrieved
        as needed.

        Parameters
        ----------
        X : np.ndarray
            Data array, where the row index distinguishes different
            observations, column index distinguishes different features.
        params : dict, optional
            See attribute params for default keys/values.

        Writes Attributes
        -----------------
        evalues : np.ndarray 
            Eigenvalues of transition matrix
        lbasis : np.ndarray
            Matrix of left eigenvectors (stored in columns).
        rbasis : np.ndarray 
             Matrix of right eigenvectors (stored in columns).
             self.rbasis is projection of data matrix on right eigenvectors, 
             that is, the projection on the diffusion components.
             these are simply the components of the right eigenvectors
             and can directly be used for plotting.
        """
        self.X = X 
        self.N = self.X.shape[0]
        # init parameters 
        for key in params.keys():
            self.params[key] = params[key]
        log.m(0,'computing diffusion map with method',
                '"'+self.params['method']+'"')
        # compute transition matrix
        self.compT()
        # compute spectral embedding
        self.embed()

    def compT(self):
        """ 
        Compute similarity matrix and transition matrix.

        Notes
        -----
        In the code, the following two parameters are set.

        alpha : float
            The density rescaling parameter of Coifman and Lafon (2006). Should
            in all practical applications equal 1: Then only the geometry of the
            data matters, not the sampling density.
        zerodiagonal : bool 
            Set the diagonal of the transition matrix to zero as suggested by
            Haghverdi et al. 2015.

        See also
        --------
        Also Haghverdi et al. (2016, 2015) and Coifman and Lafon (2006) and
        Coifman et al. (2005).
        """
        # compute distance matrix in squared Euclidian norm
        Dsq = utils.comp_distance(self.X,metric='sqeuclidean')
        if self.params['method'] == 'local':
            # choose sigma (width of a Gaussian kernel) according to the
            # distance of the kth nearest neighbor of each point, including the
            # point itself in the count
            k = self.params['k']
            # deterimine the distance of the k nearest neighbors
            indices = np.zeros((Dsq.shape[0],k),dtype=np.int_)
            distances_sq = np.zeros((Dsq.shape[0],k),dtype=np.float_)
            for irow,row in enumerate(Dsq):
                # the last item is already in its sorted position as
                # argpartition puts the (k-1)th element - starting to count from
                # zero - in its sorted position
                idcs = np.argpartition(row,k-1)[:k]
                indices[irow] = idcs
                distances_sq[irow] = np.sort(row[idcs])
            # choose sigma, the heuristic here often makes not much 
            # of a difference, but is used to reproduce the figures
            # of Haghverdi et al. (2016)
            if self.params['fix_nr_neighbors']:
                # as the distances are not sorted except for last element
                # take median
                sigmas_sq = np.median(distances_sq,axis=1)
            else:
                # the last item is already in its sorted position as
                # argpartition puts the (k-1)th element - starting to count from
                # zero - in its sorted position
                sigmas_sq = distances_sq[:,-1]/4
            sigmas = np.sqrt(sigmas_sq)
            log.mt(0,'determined k =',k,
                     'nearest neighbors of each point')
        elif self.params['method'] == 'standard':
            sigmas = self.params['sigma']*np.ones(self.N)
            sigmas_sq = sigmas**2
        # compute the symmetric weight matrix
        Num = 2*np.multiply.outer(sigmas,sigmas)
        Den = np.add.outer(sigmas_sq,sigmas_sq)
        W = np.sqrt(Num/Den)*np.exp(-Dsq/Den)
        # make the weight matrix sparse
        if not self.params['fix_nr_neighbors']:
            W[W < 1e-14] = 0
        else:
            # restrict number of neighbors to k
            Mask = np.zeros(Dsq.shape,dtype=bool)
            for irow,row in enumerate(indices):
                Mask[irow,row] = True
                for j in row:
                    if irow not in indices[j]:
                        Mask[j,irow] = True
            # set all entries that are not nearest neighbors to zero
            W[Mask==False] = 0
        log.mt(0,'computed W (weight matrix) with "fix_nr_neighbors" =',
               self.params['fix_nr_neighbors'])
        if False:
            pl.matshow(W)
            pl.title('$ W$')
            pl.colorbar()
        # zero diagonal as discussed in Haghverdi et al. (2015)
        # then the kernel does not encode a notion of similarity then anymore
        # and is not positive semidefinite anymore
        # in practice, it doesn't matter too much
        zerodiagonal = False
        if zerodiagonal:
            np.fill_diagonal(W, 0)
        # density normalisation 
        # as discussed in Coifman et al. (2005)
        # ensure that kernel matrix is independent of sampling density
        alpha = 1
        if alpha == 0:
            # nothing happens here, simply use the isotropic similarity matrix
            self.K = np.array(W)
        else:
            # q[i] is an estimate for the sampling density at point x_i
            # it's also the degree of the underlying graph
            q = np.sum(W,axis=0)
            # raise to power alpha
            if alpha != 1: 
                q = q**alpha
            Den = np.outer(q,q)
            self.K = W/Den
        log.mt(0,'computed K (anisotropic kernel)')
        if False:
            pl.matshow(self.K)
            pl.title('$ K$')
            pl.colorbar()
        # now compute the row normalization to build the transition matrix T
        # and the adjoint Ktilde: both have the same spectrum
        self.z = np.sum(self.K,axis=0)
        # the following is the transition matrix
        self.T = self.K/self.z[:,np.newaxis]
        # now we need the square root of the density
        self.sqrtz = np.array(np.sqrt(self.z))
        # now compute the density-normalized Kernel
        # it's still symmetric
        szszT = np.outer(self.sqrtz,self.sqrtz)
        self.Ktilde = self.K/szszT
        log.mt(0,'computed Ktilde (normalized anistropic kernel)')
        if False:
            pl.matshow(self.Ktilde)
            pl.title('$ \widetilde K$')
            pl.colorbar()
            pl.show()

    def embed(self,number=10):
        """ 
        Compute eigen decomposition of T.

        Parameters
        ----------
        number : int  
            Number of eigenvalues/vectors to be computed, set number = 0 if
            you need all eigenvectors.
        sym : bool
            Instead of computing the eigendecomposition of the assymetric 
            transition matrix, computed the eigendecomposition of the symmetric
            Ktilde matrix.

        Writes class members
        --------------------
        evalues : np.ndarray 
            Eigenvalues of transition matrix
        lbasis : np.ndarray
            Matrix of left eigenvectors (stored in columns).
        rbasis : np.ndarray 
             Matrix of right eigenvectors (stored in columns).
             self.rbasis is projection of data matrix on right eigenvectors, 
             that is, the projection on the diffusion components.
             these are simply the components of the right eigenvectors
             and can directly be used for plotting.
        """
        self.rbasisBool = True
        # compute the spectrum
        if number == 0:
            w,u = np.linalg.eigh(self.Ktilde)
        else:
            number = min(self.Ktilde.shape[0]-1,number)
            w,u = sp.sparse.linalg.eigsh(self.Ktilde, k=number)
        self.evalues = w[::-1]
        u = u[:,::-1]
        log.mt(0,'computed Ktilde\'s eigenvalues:')
        log.m(0,self.evalues)
        log.m(1,'computed',number,'eigenvalues. if you want more increase the'
                'parameter "number" or set it to zero, to compute all eigenvalues')
        if self.params['sym']:
            self.rbasis = self.lbasis = u
        else:
            # The eigenvectors of T are stored in self.rbasis and self.lbasis 
            # and are simple trafos of the eigenvectors of Ktilde.
            # rbasis and lbasis are right and left eigenvectors, respectively
            self.rbasis = np.array(u/self.sqrtz[:,np.newaxis])
            self.lbasis = np.array(u*self.sqrtz[:,np.newaxis])
            # normalize in L2 norm
            # note that, in contrast to that, a probability distribution
            # on the graph is normalized in L1 norm
            # therefore, the eigenbasis in this normalization does not correspond
            # to a probability distribution on the graph
            self.rbasis /= np.linalg.norm(self.rbasis,axis=0,ord=2)
            self.lbasis /= np.linalg.norm(self.lbasis,axis=0,ord=2)

    def _embed(self):
        """
        Checks and tests for embed.
        """
        if False:
            pl.semilogy(w,'x',label=r'$ \widetilde K$')
            pl.show()
        if log.verbosity > 2:
            # output of spectrum of K for comparison
            w,v = np.linalg.eigh(self.K)
            log.mi('spectrum of K (kernel)')
        if log.verbosity > 3:
            # direct computation of spectrum of T
            w,vl,vr = sp.linalg.eig(self.T,left=True)
            log.mi('spectrum of transition matrix (should be same as of Ktilde)')

    def compM(self):
        """ 
        The M matrix is the matrix that results from summing over all powers of
        T in the subspace without the first eigenspace.

        See Haghverdi et al. (2016).
        """
        # the projected inverse therefore is
        self.M = sum([self.evalues[i]/(1-self.evalues[i])
                      * np.outer(self.rbasis[:,i],self.lbasis[:,i])
                      for i in range(1,self.evalues.size)])
        log.mt(0,'computed M matrix')
        if False:
            pl.matshow(self.Ktilde)
            pl.title('Ktilde')
            pl.colorbar()
            pl.matshow(self.M)
            pl.title('M')
            pl.colorbar()
            pl.show()

    def compDdiff(self):
        """ 
        Returns the distance matrix in the diffusion metric.

        Is based on the M matrix. self.Ddiff[self.iroot,:] stores diffusion
        pseudotime as a vector.
        """
        self.Dbool = True
        self.Ddiff = sp.spatial.distance.pdist(self.M)
        self.Ddiff = sp.spatial.distance.squareform(self.Ddiff)
        log.mt(0,'computed Ddiff distance matrix')
        if False:
            pl.matshow(self.Ddiff)
            pl.title('Ddiff')
            pl.colorbar()
            pl.show()
        return self.Ddiff

    def find_root(self,xroot):
        """ 
        Determine the index of the root cell.

        Given an expression vector, find the observation index that is closest
        to this vector.
        
        Parameters
        ----------
        xroot : np.ndarray
            Vector that marks the root cell, the vector storing the initial
            condition, only relevant for computing pseudotime.
        """
        # this is the squared distance
        dsqroot = 1e10
        self.iroot = 0
        for i in range(self.N):
            diff = self.X[i,:]-xroot
            dsq = diff.dot(diff)
            if  dsq < dsqroot:
                dsqroot = dsq
                self.iroot = i
                if np.sqrt(dsqroot) < 1e-10:
                    log.m(2,'root found at machine prec')
                    break
        log.m(1,'sample',self.iroot,'has distance',np.sqrt(dsqroot),'from root')
        return self.iroot

    def set_pseudotimes(self):
        """
        Return pseudotimes with respect to root point.
        """
        self.pseudotimes = self.Ddiff[self.iroot]/np.max(self.Ddiff[self.iroot])

    def branchings_segments(self):
        """ 
        Detect branchings and partition the data into corresponding segments.

        Detect all branchings up to params['num_branchings'].

        Writes
        ------
        segs : np.ndarray
            Array of dimension (number of segments) x (number of data
            points). Each row stores a mask array that defines a segment.
        segstips : np.ndarray
            Array of dimension (number of segments) x 2. Each row stores the
            indices of the two tip points of each segment.
        segslabels : np.ndarray
            Array of dimension (number of data points). Stores an integer label
            for each segment.
        """
        self.detect_branchings()
        self.check_segments()
        self.postprocess_segments()
        self.order_segments()
        self.set_segslabels()
        self.order_pseudotime()

    def select_segment(self,segs,segstips):
        """
        Out of a list of line segments, choose segment that has the most
        distant second data point.

        Assume the distance matrix Ddiff is sorted according to seg_idcs.
        Compute all the distances.

        Returns
        -------
        iseg : int
            Index identifying the position within the list of line segments.
        tips3 : int
            Positions of tips within chosen segment.
        """
        scores_tips = np.zeros((len(segs),4))
        allindices = np.arange(self.N,dtype=int)
        for iseg, seg in enumerate(segs):
            # do not consider 'unproper segments'
            if segstips[iseg][0] == -1:
                continue
            # restrict distance matrix to points in segment
            Dseg = self.Ddiff[np.ix_(seg,seg)]
            # obtain the two indices that maximize distance in the segment
            # call them tips
            if False:
                # obtain the position within the segment by searching for
                # the maximum
                tips = list(np.unravel_index(np.argmax(Dseg),Dseg.shape))
            if True:
                # map the global position to the position within the segment
                tips = [np.where(allindices[seg] == tip)[0][0]
                        for tip in segstips[iseg]]
            # find the third point on the segment that has maximal
            # added distance from the two tip points
            dseg = Dseg[tips[0]] + Dseg[tips[1]]
            # add this point to tips, it's a third tip, we store it at the first
            # position in an array called tips3
            tips3 = np.insert(tips,0,np.argmax(dseg))
            # compute the score as ratio of the added distance to the third tip,
            # to what it would be if it were on the straight line between the
            # two first tips, given by Dseg[tips[:2]]
            # if we did not normalize with, there would be a danger of simply
            # assigning the highest score to the longest segment
            score = dseg[tips3[0]]/Dseg[tips3[1],tips3[2]]
            # write result
            scores_tips[iseg,0] = score
            scores_tips[iseg,1:] = tips3
        iseg = np.argmax(scores_tips[:,0])
        tips3 = scores_tips[iseg,1:].astype(int)
        return iseg, tips3

    def detect_branchings(self):
        """ 
        Detect all branchings up to params['num_branchings'].
        """
        log.m(0,'detect',self.params['num_branchings'],'branchings')
        # a segment is a subset of points of the data set
        # it's completely defined by the indices of the points in the segment
        # initialize the search for branchings with a single segment,
        # that is, get the indices of the whole data set
        indices_all = np.arange(self.Ddiff.shape[0],dtype=int)
        # let's keep a list of segments, the first segment to add is the 
        # whole data set
        segs = [indices_all]
        # a segment can as well be defined by the two points that have maximal
        # distance in the segment, the "tips" of the segment 
        #
        # the rest of the points in the segment is then defined by demanding
        # them to "be close to the line segment that connects the tips", that
        # is, for such a point, the normalized added distance to both tips is
        # smaller than one:
        #     (D[tips[0],i] + D[tips[1],i])/D[tips[0],tips[1] < 1
        # of course, this condition is fulfilled by the full cylindrical
        # subspace surrounding that line segment, where the radius of the
        # cylinder can be infinite
        #
        # if D denotes a euclidian distance matrix, a line segment is a linear
        # object, and the name "line" is justified. if we take the
        # diffusion-based distance matrix Ddiff, which approximates geodesic
        # distance, with "line", we mean the shortest path between two points,
        # which can be highly non-linear in the original space
        #
        # let us define the tips of the whole data set
        tips_all = list(np.unravel_index(np.argmax(self.Ddiff),self.Ddiff.shape))
        # we keep a list of the tips of each segment
        segstips = [tips_all]
        for ibranch in range(self.params['num_branchings']):
            # out of the list of segments, determine the segment
            # that most strongly deviates from a straight line
            # and provide the three tip points that span the triangle
            # of maximally distant points
            iseg, tips3 = self.select_segment(segs,segstips)
            log.m(0,'detected tip points',tips3)
            # detect branching and update segs and segstips
            segs, segstips = self.detect_branching(segs,segstips,iseg,tips3)
        # store as class members
        self.segs = segs
        self.segstips = segstips
        log.mt(0,'finished branching detection')

    def postprocess_segments(self):
        """
        Convert the format of the segment class members.
        """
        # make segs a list of mask arrays, it's easier to store 
        # as there is a hdf5 equivalent
        for iseg,seg in enumerate(self.segs):
            mask = np.zeros(self.Ddiff.shape[0],dtype=bool)
            mask[seg] = True
            self.segs[iseg] = mask
        # convert to arrays
        self.segs = np.array(self.segs)
        self.segstips = np.array(self.segstips)
        
    def check_segments(self):
        """
        Perform checks on segments and sort them according to pseudotime.
        """
        # find the segment that contains the root cell
        for iseg,seg in enumerate(self.segs):
            if self.iroot in seg:
                isegroot = iseg
                break
        # check whether the root cell is one of the tip cells of the 
        # segment, if not we need to introduce a new branching, directly
        # at the root cell
        if self.iroot not in self.segstips[iseg]:
            # if it's not exactly a tip, but very close to it, 
            # just keep it as it is
            dist_to_root = self.Ddiff[self.iroot,self.segstips[iseg]]
            # otherwise, allow branching at root
            if (np.min(dist_to_root) > 0.01*self.Ddiff[tuple(self.segstips[iseg])]
                and self.params['allow_branching_at_root']):
                allindices = np.arange(self.N,dtype=int)
                tips3_global = np.insert(self.segstips[iseg],0,self.iroot)
                # map the global position to the position within the segment
                tips3 = np.array([np.where(allindices[self.segs[iseg]] == tip)[0][0]
                                  for tip in tips3_global])
                # detect branching and update self.segs and self.segstips
                self.segs, self.segstips = self.detect_branching(self.segs,
                                                                 self.segstips,
                                                                 iseg,tips3)            
        
    def order_segments(self):
        """
        Order segments according to average pseudotime.
        """ 
        # there are different options for computing the score
        if False:
            # minimum of pseudotimes in the segment
            score = np.min
        if True:
            # average pseudotime
            score = np.average
        # score segments by minimal pseudotime
        seg_scores = []
        for seg in self.segs:
            seg_scores.append(score(self.pseudotimes[seg]))
        indices = np.argsort(seg_scores)
        # order segments by minimal pseudotime
        self.segs = self.segs[indices]
        self.segstips = self.segstips[indices]
        # within segstips, order tips according to pseudotime
        for itips, tips in enumerate(self.segstips):
            if tips[0] != -1:
                indices = np.argsort(self.pseudotimes[tips])
                self.segstips[itips] = self.segstips[itips][indices]

    def set_segslabels(self):
        """
        Return a single array that stores integer segment labels.
        """
        segslabels = np.zeros(self.Ddiff.shape[0],dtype=int)
        for iseg,seg in enumerate(self.segs):
            segslabels[seg] = iseg
        self.segslabels = segslabels

    def order_pseudotime(self):
        """
        Define ordering of points that reflects segment and pseudotime order.

        Writes
        ------
        indices : np.ndarray
            Index array of shape n, which stores an ordering of the data points
            with respect to increasing segment index and increasing pseudotime.
        changepoints : np.ndarray
            Index array of shape len(ssegs)-1, which stores the indices of
            points where the segment index changes, with respect to the ordering
            of indices.
        """ 
        # sort indices according to segments
        indices = np.argsort(self.segslabels)
        segslabels = self.segslabels[indices]
        # find changepoints of segments
        changepoints = np.arange(indices.size-1)[np.diff(segslabels)==1]+1
        pseudotimes = self.pseudotimes[indices]
        for iseg,seg in enumerate(self.segs):
            # only consider one segment, it's already ordered by segment
            seg_sorted = seg[indices]
            # consider the pseudotimes on this segment and sort them
            seg_indices = np.argsort(pseudotimes[seg_sorted])
            # within the segment, order indices according to increasing pseudotime
            indices[seg_sorted] = indices[seg_sorted][seg_indices]
        # define class members
        self.indices = indices
        self.changepoints = changepoints

    def detect_branching(self,segs,segstips,iseg,tips3):
        """ 
        Detect branching on given segment.

        Call function _detect_branching and perform bookkeeping on segs and
        segstips.

        Parameters
        ----------
        segs : list of np.ndarray
            Ddiff distance matrix restricted to segment.
        segstips : list of np.ndarray
            Stores all tip points for the segments in segs.
        iseg : int
            Position of segment under study in segs.
        tips3 : np.ndarray
            The three tip points. They form a 'triangle' that contains the data.
        
        Returns
        -------
        segs : list of np.ndarray
            Updated list of segments.
        segstips : list of np.ndarray
            Updated list of segstips.
        """
        seg = segs[iseg]
        # restrict distance matrix to points in chosen segment seg
        Dseg = self.Ddiff[np.ix_(seg,seg)]
        # given the three tip points and the distance matrix detect the
        # branching on the segment, return the list ssegs of segments that
        # are defined by splitting this segment
        ssegs, ssegs_tips = self._detect_branching(Dseg,tips3)            
        # map back to global indices
        for iseg_new,seg_new in enumerate(ssegs):
            ssegs[iseg_new] = seg[seg_new]
            if ssegs_tips[iseg_new][0] != -1:
                ssegs_tips[iseg_new] = seg[ssegs_tips[iseg_new]] 
        # remove previous segment
        segs.pop(iseg)
        segstips.pop(iseg)
        # append new segments
        segs += ssegs
        segstips += ssegs_tips
        return segs, segstips

    def _detect_branching(self,Dseg,tips):
        """ 
        Detect branching on given segment.

        Call function __detect_branching three times for all three orderings of
        tips. Points that do not belong to the same segment in all three
        orderings are assigned to a fourth segment. The latter is, by Haghverdi
        et al. (2016) referred to as 'undecided cells'.

        Parameters
        ----------
        Dseg : np.ndarray
            Ddiff distance matrix restricted to segment.
        tips : np.ndarray
            The three tip points. They form a 'triangle' that contains the data.
        
        Returns
        -------
        ssegs : list of np.ndarray
            List of segments obtained from splitting the single segment defined
            via the first two tip cells.
        ssegstips : list of np.ndarray
            List of tips of segments in ssegs.
        """
        if False:
            ssegs = self._detect_branching_versions(Dseg,tips)
        if True:
            ssegs = self._detect_branching_single(Dseg,tips)
        # make sure that each data point has a unique association with a segment
        masks = np.zeros((3,Dseg.shape[0]),dtype=bool)
        for iseg,seg in enumerate(ssegs):
            masks[iseg][seg] = True
        nonunique = np.sum(masks,axis=0) > 1
        # obtain the corresponding index arrays from masks
        ssegs = []
        for iseg,mask in enumerate(masks):
            mask[nonunique] = False
            ssegs.append(np.arange(Dseg.shape[0],dtype=int)[mask])
        # compute new tips within new segments
        ssegstips = []
        for inewseg, newseg in enumerate(ssegs):
            # get tip point position within segment
            tip = np.where(np.arange(Dseg.shape[0])[newseg]
                           == tips[inewseg])[0][0]
            # new tip within restricted distance matrix
            secondtip = np.argmax(Dseg[np.ix_(newseg,newseg)][tip])
            # map back to position within segment
            secondtip = np.arange(Dseg.shape[0])[newseg][secondtip]
            # add to list
            ssegstips.append([tips[inewseg],secondtip])
        # for the points that cannot be assigned to the three segments of the
        # branching, hence have no tip cells, but form a subset of their own,
        # add dummy tips [-1,-1]
        # this is not a good solution, but it ensures that we can easily write
        # to hdf5 as ssegstips can be transformed to np.ndarray with dtype = int
        ssegstips.append(np.array([-1,-1]))
        # the following would be preferrable, but then ssegstips results in 
        # a np.ndarray with dtype = object, for which there is no straight
        # forward hdf5 format, a solution via masks seems too much work
        #     ssegstips.append(np.array([],dtype=int))
        # also add the points not associated with a clear seg to ssegs
        mask = np.zeros(Dseg.shape[0],dtype=bool)
        # all points assigned to segments (flatten ssegs)
        mask[[i for l in ssegs for i in l]] = True
        # append all the points that have not been assigned. in Haghverdi et
        # al. (2016), we call them 'undecided cells'
        ssegs.append(np.arange(Dseg.shape[0],dtype=int)[mask==False])

        return ssegs, ssegstips

    def _detect_branching_single(self,Dseg,tips):
        """ 
        Detect branching on given segment.
        """
        # compute branchings using different starting points the first index of
        # tips is the starting point for the other two, the order does not
        # matter
        ssegs = []
        # permutations of tip cells
        ps = [[0,1,2], # start by computing distances from the first tip 
              [1,2,0], #             -"-                       second tip
              [2,0,1], #             -"-                       third tip
              ]
        for i,p in enumerate(ps):
            ssegs.append(self.__detect_branching(Dseg,
                                                 tips[p])[0])
        return ssegs

    def _detect_branching_versions(self,Dseg,tips):
        """ 
        Detect branching on given segment using three different versions.
        """
        # compute branchings using different starting points the first index of
        # tips is the starting point for the other two, the order does not
        # matter
        ssegs_versions = []
        # permutations of tip cells
        ps = [[0,1,2], # start by computing distances from the first tip 
              [1,2,0], #             -"-                       second tip
              [2,0,1], #             -"-                       third tip
              ]
        # invert permutations
        inv_ps = [[0,1,2],
                  [2,0,1],
                  [1,2,0],
                  ]
        for i,p in enumerate(ps):
            ssegs = self.__detect_branching(Dseg,
                                            tips[p])
            ssegs_versions.append(np.array(ssegs)[inv_ps[i]])
        ssegs = []
        # run through all three assignments of segments, and keep
        # only those assignments that were found in all three runs        
        for inewseg, newseg_versions in enumerate(np.array(ssegs_versions).T):
            if len(newseg_versions) == 3:
                newseg = np.intersect1d(np.intersect1d(newseg_versions[0],
                                                       newseg_versions[1]),
                                        newseg_versions[2])
            else:
                newseg = newseg_versions[0]
            ssegs.append(newseg)
        
        return ssegs
        

    def __detect_branching(self,Dseg,tips):
        """ 
        Detect branching on given segment.

        Compute point that maximizes kendall tau correlation of the sequences of
        distances to the second and the third tip, respectively, when 'moving
        away' from the first tip: tips[0]. 'Moving away' means moving in the
        direction of increasing distance from the first tip.

        Parameters
        ----------
        Dseg : np.ndarray
            Ddiff distance matrix restricted to segment.
        tips : np.ndarray
            The three tip points. They form a 'triangle' that contains the data.

        Returns
        -------
        ssegs : list of np.ndarray
            List of segments obtained from splitting the single segment defined
            via the first two tip cells.
        """
        # sort distance from first tip point
        idcs = np.argsort(Dseg[tips[0]])
        # then the sequence of distances Dseg[tips[0]][idcs] increases
        # consider now the sequence of distances from the other 
        # two tip points, which only increase when being close to tips[0]
        # where they become correlated
        # at the point where this happens, we define a branching point
        if True:
            imax = self.kendall_tau_split(Dseg[tips[1]][idcs],
                                             Dseg[tips[2]][idcs])
        if False:
            # if we were in euclidian space, the following should work
            # as well, but here, it doesn't because the scales in Dseg are
            # highly different, one would need to write the following equation
            # in terms of an ordering, such as exploited by the kendall 
            # correlation method above
            imax = np.argmin(Dseg[tips[0]][idcs] 
                                + Dseg[tips[1]][idcs]
                                + Dseg[tips[2]][idcs])
        # init list to store new segments
        ssegs = []
        # first new segment: all points until, but excluding the branching point
        ibranch = imax + 1
        # ibranch = int(0.95*imax) # more conservative here
        ssegs.append(idcs[:ibranch])
        # define nomalized distances to tip points for the rest of the data
        dist1 = Dseg[tips[1],idcs[ibranch:]]/Dseg[tips[1],idcs[ibranch-1]]
        dist2 = Dseg[tips[2],idcs[ibranch:]]/Dseg[tips[2],idcs[ibranch-1]]
        # assign points according to whether being closer to tip cell 1 or 2
        ssegs.append(idcs[ibranch:][dist1 <= dist2])
        ssegs.append(idcs[ibranch:][dist1 > dist2])

        return ssegs

    def kendall_tau_split(self,a,b):
        """
        Return splitting index that maximizes correlation in the sequences.

        Compute difference in Kendall tau for all splitted sequences. 

        For each splitting index i, compute the difference of the two
        correlation measures kendalltau(a[:i],b[:i]) and
        kendalltau(a[i:],b[i:]).

        Returns the splitting index that maximizes 
            kendalltau(a[:i],b[:i]) - kendalltau(a[i:],b[i:])        

        Parameters
        ----------
        a, b : np.ndarray
            One dimensional sequences. 

        Returns
        -------
        i : int
            Splitting index according to above description.
        """
        if a.size != b.size:
            raise ValueError('a and b need to have the same size')
        if a.ndim != b.ndim != 1:
            raise ValueError('a and b need to be one-dimensional arrays')            

        min_length = 5
        n = a.size
        idx_range = np.arange(min_length,a.size-min_length-1,dtype=int)
        corr_coeff = np.zeros(idx_range.size)
        pos_old = sp.stats.kendalltau(a[:min_length],b[:min_length])[0]
        neg_old = sp.stats.kendalltau(a[min_length:],b[min_length:])[0]
        for ii,i in enumerate(idx_range):
            if True:
                # compute differences in concordance when adding a[i] and b[i]
                # to the first subsequence, and removing these elements from
                # the second subsequence
                diff_pos, diff_neg = self._kendall_tau_diff(a,b,i)
                pos = pos_old + self._kendall_tau_add(i,diff_pos,pos_old)
                neg = neg_old + self._kendall_tau_subtract(n-i,diff_neg,neg_old)
                pos_old = pos
                neg_old = neg
            if False:
                # computation using sp.stats.kendalltau, takes much longer!
                # just for debugging purposes
                pos = sp.stats.kendalltau(a[:i+1],b[:i+1])[0]
                neg = sp.stats.kendalltau(a[i+1:],b[i+1:])[0]
            if False:
                # the following is much slower than using sp.stats.kendalltau,
                # it is only good for debugging because it allows to compute the
                # tau-a version, which does not account for ties, whereas
                # sp.stats.kendalltau computes tau-b version, which accounts for
                # ties
                pos = sp.stats.mstats.kendalltau(a[:i],b[:i],use_ties=False)[0]
                neg = sp.stats.mstats.kendalltau(a[i:],b[i:],use_ties=False)[0]
            corr_coeff[ii] = pos - neg
        iimax = np.argmax(corr_coeff)
        imax = min_length + iimax
        corr_coeff_max = corr_coeff[iimax]
        if corr_coeff_max < 0.3:
            log.m(1,'  -> is root itself, never obtain significant correlation')
        return imax
        
    def _kendall_tau_add(self,len_old,diff_pos,tau_old):
        """
        Compute Kendall tau delta.

        The new sequence has length len_old + 1.

        Parameters
        ----------
        len_old : int
            The length of the old sequence, used to compute tau_old.
        diff_pos : int
            Difference between concordant and non-concordant pairs.
        tau_old : float
            Kendall rank correlation of the old sequence.
        """
        return 2./(len_old+1)*(float(diff_pos)/len_old-tau_old)

    def _kendall_tau_subtract(self,len_old,diff_neg,tau_old):
        """
        Compute Kendall tau delta.

        The new sequence has length len_old - 1.

        Parameters
        ----------
        len_old : int
            The length of the old sequence, used to compute tau_old.
        diff_neg : int
            Difference between concordant and non-concordant pairs.
        tau_old : float
            Kendall rank correlation of the old sequence.
        """
        return 2./(len_old-2)*(-float(diff_neg)/(len_old-1)+tau_old)

    def _kendall_tau_diff(self,a,b,i):
        """
        Compute difference in concordance of pairs in split sequences.

        Consider splitting a and b at index i.

        Parameters
        ----------
        a, b : np.ndarray

        Returns
        -------
        diff_pos, diff_neg : int, int
            Difference between concordant and non-concordant pairs for both
            subsequences.
        """
        # compute ordering relation of the single points a[i] and b[i] 
        # with all previous points of the sequences a and b, respectively
        a_pos = np.zeros(a[:i].size,dtype=int)
        a_pos[a[:i]>a[i]] = 1
        a_pos[a[:i]<a[i]] = -1
        b_pos = np.zeros(b[:i].size,dtype=int)
        b_pos[b[:i]>b[i]] = 1
        b_pos[b[:i]<b[i]] = -1
        diff_pos = np.dot(a_pos,b_pos).astype(float)
        
        # compute ordering relation of the single points a[i] and b[i] 
        # with all later points of the sequences
        a_neg = np.zeros(a[i:].size,dtype=int)
        a_neg[a[i:]>a[i]] = 1
        a_neg[a[i:]<a[i]] = -1
        b_neg = np.zeros(b[i:].size,dtype=int)
        b_neg[b[i:]>b[i]] = 1
        b_neg[b[i:]<b[i]] = -1
        diff_neg = np.dot(a_neg,b_neg)

        return diff_pos, diff_neg

    def subsequences(self):
        """ 
        Extract subsequences. Experimental beta feature.
        """
        pt_order, pseudotimes, X_ordered = self.sort_pseudotime()
        Ddiff = np.array(self.Ddiff)
        rbasis = np.array(self.rbasis)
        lbasis = np.array(self.lbasis)
                
        i = 0
        indices = [[i],[]]
        for count in range(Ddiff.shape[0]):
            jmin = np.argmin(Ddiff[i,i+1:])
            i = i + 1 + jmin
            indices[0].append(i)
            if i >= Ddiff.shape[0]-1:
                break

        mask = np.ones(Ddiff.shape[0],dtype=bool)
        mask[indices[0]] = False
        tosort = np.arange(Ddiff.shape[0],dtype=int)[mask]

        isub = 0
        indices[1].append(isub) 
        for count in range(tosort.size):
            jmin = np.argmin(Ddiff[tosort[isub],tosort[isub+1:]])
            isub = isub + 1 + jmin
            indices[1].append(isub)
            if isub >= tosort.size-1:
                break            
        # make global indices from subindices
        indices[1] = tosort[indices[1]]

        iall = np.arange(Ddiff.shape[0])
        all = Ddiff[iall[:-1],iall[1:]]
        ibranch0 = indices[0]
        branch0 = Ddiff[indices[0][:-1],indices[0][1:]]
        pl.figure()
        pl.scatter(iall[:-1],all,marker='.')
        pl.scatter(ibranch0[:-1],branch0,marker='x')

        pl.figure()
        pl.plot(np.arange(len(indices[0])),X_ordered[indices[0]])
        pl.plot(len(indices[0])+np.arange(len(indices[1])),X_ordered[indices[1]])
        pl.show()

def run(function,args):
    """ 
    Perform DPT analysis or compute Diffusion Map.

    Writes an output file 'out/exkey..._function.ext', which stores the results
    of a dpt analysis. File type is determined by variable log.extd allowed are
    'h5' (hdf5), 'xlsx' (Excel) or 'csv' (comma separated value file).

    If called twice with the same settings, plotting is based on the existing
    output file.

    Parameters
    ----------
    function : str
        Either 'diffmap' or 'analysis'.
    args : dict
        Arguments supplied at the command line. Just two keys are needed.
        exkey : str
            String naming the output file, written to 'out/exkey...'
        subsample : int
            Integer giving the inverse fraction of sampled data points.
    """

    sig = log.suffix + '_ss{:02}'.format(args['subsample'])
    basename = args['exkey'] + sig + '_' + function 
    outfile = 'out/'+basename+'.'+log.extd

    # get the module for preprocessing the data
    if 'preprocess' in args:
        preprocess = args['preprocess']
    else:
        preprocess = scanpy.preprocess
    # always read the data from the datafile
    example, ddata = preprocess.get_data(args['exkey'])
    X = ddata['X']
    # if not specified otherwise, subsampling takes the same seed and yields a
    # reproducible result
    X, rows = utils.subsample(X,args['subsample'])

    if os.path.exists(outfile) and not log.recompute:
        log.m(0,'reading file',outfile)
        dfunc = utils.read_file_to_dict(outfile,log.extd)
    else:
        log.m(0,'X has shape',X.shape[0],'x',X.shape[1])
        if 'dpt' in example:
            params = example['dpt']
        else:
            params = {}
        if function == 'dpt':
            dfunc = dpt(X,ddata['xiroot'],params)
        elif function == 'diffmap':
            dfunc = diffmap(X,params)
        else:
            raise ValueError('either call "dpt" or "diffmap"')
        dfunc['rownames'] = ddata['rownames'][rows]
        dfunc['colnames'] = ddata['colnames']
        dfunc['basename'] = basename
        utils.write_dict_to_file(outfile,dfunc,log.extd)

    return X, dfunc

def plot_diffmap(ddmap):
    """
    Plots the results of a diffusion map calculation.
    """

    colors = 'blue'
    if 'moignard15' in ddmap['basename']:
        colors = scanpy.preprocess._moignard15_colors_labels(ddmap['rownames'])

    highlights = []

    utils.plot_diffmap(ddmap['rbasis'][:,[1,2,3]],
                       c=colors,
                       highlights=highlights,
                       layout='unfolded 3d')

    if log.savefigs:
        pl.savefig('figs/'+ddmap['basename']+'.'+log.extf)
    else:
        pl.show()

def plot_dpt(X,ddpt):
    """
    Plots the results of a DPT analysis.
    """

    layout = '2d'

    # color by pseudotime and by segments
    colors = [
        ddpt['pseudotimes'],
        ddpt['seglabels'],
    ]

    # example-specific settings
    if 'moignard15' in ddpt['basename']:
        # coloring according to experimental labels
        colors.append(scanpy.preprocess._moignard15_colors_labels(ddpt['rownames']))
    if 'krumsiek11' in ddpt['basename']:
        layout = '3d'
    if 'nestorowa16' in ddpt['basename']:
        layout = '2d'

    # highlight root
    highlights = list(ddpt['iroot'])
    # highlight tip points of each segment
    if False:
        highlights = [i for l in ddpt['segtips'] for i in l if l[0] != -1]

    # a single figure for all colors using 2 diffusion components
    if layout in ['2d','3d']:
        utils.plot_diffmap([ddpt['rbasis'][:,[1,2,3]]],
                            c=colors,
                            highlights=highlights,
                            layout=layout)
    # one figure per coloring, 3d and 2d projections
    else:
        for color in colors:
            utils.plot_diffmap([ddpt['rbasis'][:,[1,2,3]]],
                                c=color,
                                highlights=highlights,
                                layout='unfolded 3d')
    if log.savefigs:
        pl.savefig('figs/'+ddpt['basename']+'_diffmap.'+log.extf)

    if 11 < X.shape[1] < 50:
        # plot time series as heatmap, as in Haghverdi et al. (2016), Fig. 1d
        utils.plot_timeseries_as_heatmap(X[ddpt['indices'],:40],ddpt['colnames'])
        if False:
            pl.savefig('figs/'+ddpt['basename']+'_heatmap.pdf')

    if X.shape[1] <= 11:
        # plot time series as genes vs expression
        utils.plot_timeseries(X[ddpt['indices'],:12],ddpt['colnames'],
                              highlightsX=ddpt['changepoints'],
                              xlim=[0,1.3*X.shape[0]])
        if log.savefigs:
            pl.savefig('figs/'+ddpt['basename']+'_vsorder.'+log.extf)

    # plot segments and pseudotimes
    pl.figure()
    pl.subplot(211)
    utils.plot_timeseries_subplot(ddpt['seglabels'][ddpt['indices'],np.newaxis],
                          c=ddpt['seglabels'][ddpt['indices']],
                          highlightsX=ddpt['changepoints'],
                          ylabel='segments',
                          yticks=(np.arange(ddpt['segs'].shape[0],dtype=int) if 
                                  ddpt['segs'].shape[0] < 5 else None))

    pl.subplot(212)
    utils.plot_timeseries_subplot(ddpt['pseudotimes'][ddpt['indices'],np.newaxis],
                          c=ddpt['pseudotimes'][ddpt['indices']],
                          highlightsX=ddpt['changepoints'],
                          ylabel='pseudotime',
                          yticks=[0,1])
    if log.savefigs:
        pl.savefig('figs/'+ddpt['basename']+'_segpt.'+log.extf)
    else:
        pl.show()

def main(args):
    """
    Run and plot DPT analysis based on command line arguments.
    """

    X, ddpt = run('dpt',args)

    plot_dpt(X,ddpt)

def main_diffmap(args):
    """
    Run and plot Diffusion Maps based on command line arguments.
    """

    X, ddmap = run('diffmap',args)

    plot_diffmap(ddmap)


if __name__ == '__main__':

    args = utils.read_args(__doc__,preprocess.examples)
    main(args)



