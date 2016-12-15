"""

Plotting

From package scanpy (single-cell analysis in python).

Written in Python 3 (compatible with 2) using Numpy, Scipy, Matplotlib
(use Anaconda package for all modules).

Copyright (c) 2016 F. Alexander Wolf (http://falexwolf.de).
   
"""   

# standard modules
import os
import sys
import glob
import argparse
# scientific modules
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.lines as mlines
from matplotlib import rcParams
from matplotlib import ticker
    
# color dict
cDict = {
    '0':'black',
    # red
    '1':'#8B0000',  # darkred
    '11':'#EE0000', # red
    # orange
     '2':'#DD2F00',   # darkorange1                                        
     '21':'#DD8C00',  # darkorange                                         
     '22':'#FF7F50',  # darkorange                            
     '23':'#FFA500',  # orange                                             
     '24':'#FFBB00',  # orange                                             
    # green
    '3':'#006400',   # darkgreen
    '31':'#556B2F',  # DarkOliveGreen
    '32':'#228B22',  # forestgreen
    '33':'#66CD00',  # chartreuse3
    '34':'#42CD42',  # limegreen
    '35':'#7CFC00',  # LawnGreen
    # blue
    '4':'#00008B',   # darkblue
    '41':'#104E8B',  # DodgerBlue4
    '42':'#1874CD',  # dodgerblue3
    '43':'#1E90FF',  # DodgerBlue
    '44':'#1E90FF',  # DodgerBlue
    '45':'#00BFFF',  # DeepSkyBlue
    # violett
    '5':'#68228B',  # DarkOrchid4
    '51':'#9932CC', # darkorchid
    '52':'#8B008B', # darkmagenta
    '53':'#FF34B3', # maroon1
    # yellow
    '6':'#FFA500',   # orange
    '61':'#FFB90F',  # DarkGoldenrod1
    '62':'#FFD700',  # gold
    '63':'#FFFF00',  # yellow
    # other
    '7':'turquoise', # turquoise
    '8':'#212121',   # darkgrey
    '81':'#424242',
    '82':'grey',
    '99':'white',
}
# list of colors
cl = [ cDict[k] for k in ['4','3','1','21','42','33','11',
                          '53','23','35','7','81','0']]
# list of colors (rainbow)
clrb = [ cDict[k] for k in ['1','11',                 # red (dark to light)
                            '2','21','23',            # orange (dark to light)
                            '61','62',                # yellow (dark to light)
                            '35','34','33','3',       # green (light to dark)
                            '45','43','42','41','4',  # blue (light to dark)
                            '5',                      # violet
                            '8'                       # dark grey 
                            ]]
# standard linewidth
lw0 = rcParams['lines.linewidth']
# list of markers
ml = ['o','s','^','d']

# init default values for rcParams
def init_figpar():
    # args figure
    rcParams['figure.figsize'] = (5,4)
    rcParams['figure.subplot.left'] = 0.18
    rcParams['figure.subplot.right'] = 0.96
    rcParams['figure.subplot.bottom'] = 0.15
    rcParams['figure.subplot.top'] = 0.91

    # args font
    rcParams['font.family'] = ['Arial','Helvetica']
    fontsize = 14
    rcParams['font.size'] = fontsize
    rcParams['legend.fontsize'] = 0.92*fontsize
    rcParams['axes.titlesize'] = fontsize

    # legend
    rcParams['legend.numpoints'] = 1

    rcParams['savefig.dpi'] = 200


init_figpar()



# for use as plotting tool of data files
if __name__ == '__main__':

    p = argparse.ArgumentParser(
        description='Plot data files.')
    aa = p.add_argument
    aa('filenames',type=str,nargs='+',
                       help='filenames of datafiles (wildcards are allowed)')
    aa('--type',type=str,default='standard',
                       help='choose type of plot [standard, matrix]')
    aa('--full',dest='full',action='store_const',const=True,default=False,
                       help='plot all columns of the data file versus the first')
    aa('--crainbow',dest='crainbow',action='store_const',const=True,default=False,
                       help='order colors as in rainbow')
    aa('--transpose',action='store_const',const=True,default=False,
                       help='transpose data array before manipulating it')
    aa('--scatter',dest='scatter',action='store_const',const=True,default=False,
                       help='use scatter plot to plot array that represents grid')
    aa('--mataspect', dest='mataspect', default=[1,1,1,1,1,1], type = float, nargs='*',
                       help='aspect ratio for mataspect')
    aa('--cshrink',dest='cshrink',type=float,default=1,
                       help='shrink the colorbar')
    aa('--cticks',dest='cticks',type=float,default=[],nargs='*',
                       help='list of cticks')
    aa('--yticks',dest='yticks',type=float,default=[],nargs='*',
                       help='list of yticks')
    aa('--xticks',dest='xticks',type=float,default=[],nargs='*',
                       help='list of xticks')
    aa('--xticklabels',dest='xticklabels',type=str,default=[],nargs='*',
                       help='list of xticklabels that matches xticks')
    aa('--noshow',dest='noshow',action='store_const',const=True,default=False,
                       help='do not actually plot anything' + 
                            '[only meaningful if computation is done instead]')
    aa('--coly', dest='coly', type=int, default=[], nargs='*',
                       help='specify columns for the y values [1]')
    aa('--colslice', dest='colslice', type=int, default=[], nargs='*',
                       help='if the value changes in the column colslice, then start new line')
    aa('--contourlevels', dest='contourlevels', type=float, default=[], nargs='*',
                       help='specify contourlevels')
    aa('--colx', dest='colx', type=int, default=[], nargs='*',
                       help='specify a column for the x values [0]')
    aa('--pts', dest='pts',  action='store_const', const=True, default=False,
                       help='use points instead of lines for display')
    aa('--ls', dest='ls', default=[], type = str, nargs='*',
                       help='specify line style' +  
                            '[allowed values {-}, . or list like [{-} {.} {:} {--}] ]')
    aa('--lw', dest='lw', default=[1], type = float, nargs='*',
                       help='specify line width in list like fashion [1 1.5 1 0.7 ]')
    aa('--ms', dest='ms', default=[], type = float, nargs='*',
                       help='marker size')
    aa('--c', dest='c', default=[0], type = int, nargs='*',
                       help='specify color in list like fashion [1 2 1 4 ]')
    aa('--zorder', dest='zorder', default=[], type = int, nargs='*',
                       help='specifyzorder in list like fashion [ 0 2 1 1 4 ]')
    aa('--log', dest='log',action='store_const', const=True, default=False,
                       help='log log plot')
    aa('--logy', dest='logy',  action='store_const', const=True, default=False,
                       help='set log scale for y axis')
    aa('--logx', dest='logx',  action='store_const', const=True, default=False,
                       help='set log scale for x axis')
    aa('--abs', dest='abs',  action='store_const', const=True, default=False,
                       help='plot absolute values')
    aa('--save', dest='save', default='', type=str,
                       help='save plot under filename specfied after save ')
    aa('--legend', dest='legend', type=str, default=[], nargs='*',
                       help='specify the legend entries, enforce no legend by "none"')
    aa('--legendloc', dest='legendloc', type=str, default='upper right',
                       help='specify legend location, specify "none" if no legend should be displayed')
    aa('--legendframe',dest='legendframe',action='store_const',const=True,default=False,
                       help='put white background for legend')
    aa('--handlelen', dest='handlelen', type=float, default=None,
                       help='specify handlelength of legend')
    aa('--handletextpad', dest='handletextpad', type=float, default=None,
                       help='specify handletextpad')
    aa('--xlabel', dest='xlabel', type=str, default='',
                       help='specify xlabel')
    aa('--ylabel', dest='ylabel', type=str, default='',
                       help='specify ylabel')
    aa('--clabel', dest='clabel', type=str, default='',
                       help='specify clabel (colorbar)')
    aa('--figsize', dest='figsize', type=float, default=[], nargs=2,
                       help='x and y length of figure')
    aa('--xlim', dest='xlim', type=float, default=[], nargs='*',
                       help='limits for x [expects 2 values]')
    aa('--ylim', dest='ylim', type=float, default=[], nargs='*',
                       help='limits for y [expects 2 values]')
    aa('--lb', dest='lb', type=float, default=0.,
                       help='left border to change rcParams')
    aa('--rb', dest='rb', type=float, default=0.,
                       help='right border to change rcParams')
    aa('--bb', dest='bb', type=float, default=0.,
                       help='bottom border')
    aa('--tb', dest='tb', type=float, default=0.,
                       help='top border')
    aa('--xnum', dest='xnum',  action='store_const', 
                        const=True, default=False, 
                       help='plot versus data index')
    aa('--add', dest='add',  action='store_const', const=True, default=False, 
                       help='add the first two ydata values')
    aa('--label', dest='label', type=str, default=[], nargs='*',
                       help='specify one or more labels with position [LABEL1 posx1 posy1 LABEL2 posx2 posy2 ...]')
    aa('--fs', dest='fs', type=float, default=18,
                       help='global font size')
    aa('--labelfs', dest='labelfs', type=float, default=1.1,
                       help='specfiy relative label fontsize as float [default 1.1]')
    aa('--legendfs', dest='legendfs', type=float, default=0.92,
                       help='specfiy relative label fontsize as float [default 0.92]')

    args = p.parse_args()

    filenames = args.filenames
    if '*' in filenames[0]:
        filenames = glob.glob(filenames[0])
        filenames = sorted(filenames)
        for f in filenames:
            print(f),
        print

    def loadtxtH(filename):
        """ return header as string """
        header=''
        for line in open(filename):
            if line.startswith("#"):
                header += line
            else:
                break
        return np.loadtxt(filename), header

    colx = args.colx
    coly = args.coly

    if args.figsize != []:
            rcParams['figure.figsize'] = args.figsize[0],args.figsize[1]
    rcParams['font.size'] = args.fs
    rcParams['legend.fontsize'] = args.legendfs*args.fs

    # define line styles
    if len(args.ls) == 1:
            lis = [ args.ls[0].strip('{}') for i in range(30)] 
    elif len(args.ls) > 1: 
            lis = [ l.strip('{}') for l in args.ls]
    elif args.full:
            lis = ['.' for i in range(30)]
    elif args.pts: 
            lis = ['.','+','x','d','s','>']
    else: 
            lis = ['-' for i in range(30)]

    # define boundarys of image
    if args.lb > 0.: rcParams['figure.subplot.left'] = args.lb
    if args.rb > 0.: rcParams['figure.subplot.right'] = args.rb
    if args.bb > 0.: rcParams['figure.subplot.bottom'] = args.bb
    if args.tb > 0.: rcParams['figure.subplot.top'] = args.tb

    # define line width
    lw0 = rcParams['lines.linewidth']
    lws = [args.lw[0]*lw0 for i in range(50)]
    if len(args.lw) > 1: lws = [a*lw0 for a in args.lw]

    # define marker size
    ms0 = 7
    msl = [ms0 for i in range(50)]
    if args.ms != []:
        msl = args.ms

    #new
    if args.crainbow:
        cl = clrb
    if len(args.c)>1: cl = [mpl_plot.cl[i] for i in args.c]

    shiftx = 0
    shifty = 0
    scalex = 1
    scaley = 1

    checkAllRequiredDataFilesExistent = True

    #####################################################################
    # postprocessing of plot
    #####################################################################

    def pimpaxis():
        ylabelunit = ''
        if not args.log and not args.logy:
            if abs(pl.ylim()[1])>10000 or (abs(pl.ylim()[1]) < 3e-2 and pl.ylim()[1] > 0):
                #print(pl.ylim()[1])
                scale_pow = -int(np.floor(np.log10(pl.ylim()[1])))
                pl.gca().get_yaxis().set_major_formatter(
                    ticker.FuncFormatter(lambda x,p : 
                                                             ("%.2f"%(x*(10**scale_pow))).rstrip('0').rstrip('.')))
                ylabelunit = r'$\times\,10^{'+str(-scale_pow)+'}$'
                pl.ylabel(pl.gca().get_ylabel()+ylabelunit)
            else:
                # replace trailing zeros for ticks
                pl.gca().get_yaxis().set_major_formatter(
                    ticker.FuncFormatter(lambda x,p : ("%.3f"%(x)).rstrip('0').rstrip('.')))
        if not args.logx:
            pl.gca().get_xaxis().set_major_formatter(
                ticker.FuncFormatter(lambda x,p : ("%.3f"%(x)).rstrip('0').rstrip('.')))

    def postprocess(ax):

        if args.label != []:
            for i in range(len(args.label)/3):
                fig.text(float(args.label[3*i+1]),float(args.label[3*i+2]),
                            args.label[3*i],fontsize=args.labelfs*args.fs)

        if args.yticks != []:
            pl.yticks(args.yticks)

        if args.xticks != []:
            pl.xticks(args.xticks)

        if args.xticklabels != []:
            ax.set_xticklabels(args.xticklabels)

    def showorsavefig(block=True):

        if args.save != '':
            # save the python script
            if args.savescript: os.system('cp ' + sys.argv[0] + ' ' + args.save[:-4]+'.py')        
            # make directory
            if not os.path.exists(os.path.dirname(args.save)):
                os.system('mkdir ' + os.path.dirname(args.save))
            # save the figure
            pl.savefig(args.save,dpi=args.dpi)
            # save the data
            tarfile = args.save[:-4]+'.tar'
            tarfileCompress = tarfile.replace('.tar','.tar.gz')
            if os.path.exists(tarfileCompress) and checkAllRequiredDataFilesExistent:
                print('remove tarfile and build it again')
                os.system('rm ' + tarfileCompress)
            # only if the tarfile does not exist (either as it is the first time
            # this script is called, or as it has been removed in the preceding line)
            if not os.path.exists(tarfileCompress): 
                #    os.system('tar -rf ' + tarfile 
                #                         + ' -C '+ os.path.dirname(args.save) 
                #                         + ' ' + os.path.basename(args.save))    
                for filename in filenames:
                    os.system('tar -rf ' + tarfile 
                                        + ' -C '+ os.path.dirname(filename) 
                                        + ' ' + os.path.basename(filename))
                # zip tarfile
                os.system('gzip ' + tarfile) 
            else:
                print('original data is no longer there, but stored in tarfile')
            plotfile = args.save[:-4]+'.plot'
            f = file(plotfile,'w') 
            string = ''
            for arg in sys.argv: 
                string += (arg + ' ' 
                                     if (sum([c in arg for c in ['(',')','$',' ']])==0 
                                             and arg!='') else '\"' + arg + '\" ')
            # f.write('echo ' + string + '\n')
            f.write(string+'\n')
            f.close()
            os.system('chmod +x ' + plotfile)
        #    os.system('tar -rf ' + tarfile + ' -C '+ os.path.dirname(plotfile) + ' ' + os.path.basename(plotfile))    
            if sys.platform == "linux" or sys.platform == "linux2":
                if '.pdf' == args.save[-4:]:
                    os.system('xpdf ' + args.save + ' & ')
                else:
                    os.system('eog ' + args.save + ' & ')
            elif sys.platform == "darwin":
                os.system('open -a Preview ' + args.save)

        if not args.noshow and args.save == '': pl.show(block=block)
        quit()

    def standardplot():
        """
        Generate standard plot.
        """
        fig = pl.figure()
        ax = pl.axes()
        axMain = ax
        #
        data, header = loadtxtH(filenames[0])
        if len(header) > 0:
            labels = header.split('\n')[-2].strip('#').split()
        else:
            labels = []
        if args.full:
            filename = filenames[0]
            if not os.path.exists(filename) and args.save!= '':
                checkAllRequiredDataFilesExistent = False
                tarfile = args.save[:-4]+'.tgz'
                if not os.path.exists(tarfile): 
                    tarfile = args.save[:-4]+'.tar.gz'
                    if not os.path.exists(tarfile): 
                        print('missing file '+filename+' and the tarfile',tarfile)
                        quit()
                if not os.path.exists(tarfile[:-4]): 
                    os.system('../10_toolkit/gunzip ' + tarfile)
                filename = args.save[:-4] + '/' + os.path.basename(filename)
            numberOfCurves = len(data[0,:])-1
        else:
            numberOfCurves = max(len(filenames),len(args.coly))
            if len(args.coly) > 1 and len(filenames) > 1 and len(args.coly) != len(filenames):            
                raise ValueError(('provide either one coly for all files' 
                                                    +'or the same number of coly arguments as filenames'))
        #
        for i in range(numberOfCurves):
            # run through list of filenames
            if not args.full and numberOfCurves==len(filenames):
                filename = filenames[i]
            # for each plot, take the first filename
            elif not args.full:
                filename = filenames[0]
            # check that all data files are there if we want to save
            if not os.path.exists(filename) and args.save != '':
                checkAllRequiredDataFilesExistent = False
                tarfile = args.save[:-4]+'.tgz'
                if not os.path.exists(tarfile): 
                    tarfile = args.save[:-4]+'.tar.gz'
                    if not os.path.exists(tarfile): 
                        print('missing file '+filename+' and the tarfile',tarfile)
                        quit()
                if not os.path.exists(tarfile[:-4]): 
                    os.system('../10_toolkit/gunzip ' + tarfile)
                filename = args.save[:-4] + '/' + os.path.basename(filename)
            # define parameters to manage data
            zorder = (args.zorder[i] if i<len(args.zorder) else 0)
            colx = (args.colx[i] if len(args.colx) > 1
                                                    else args.colx[0] if len(args.colx) == 1 else 0)
            coly, color_index = ( (args.coly[i], args.coly[i]-1) if len(args.coly) > 1 
                                   else (args.coly[0], i) if len(args.coly) == 1 else (1, i))
            if args.full:
                coly = i+1
            # load data
            if not args.full:
                data = np.loadtxt(filename)
                # treat case of only one row in the file
                if len(data.shape)==1:
                    data = np.array([data])
                if len(data)==0:
                    break
            # slice data according to values provided in args.colslice
            if args.colslice != []:
                dold = data[0,args.colslice[0]]
                itold = 0
                j = 0
                for it,d in enumerate(data[:,args.colslice[0]]):
                    if d != dold or it == len(data)-1:                     
                        if j == i:
                            print('consider slice',itold,'to',it,'value',dold)
                            data = np.array(data[itold : (it if it != len(data)-1 else len(data))])
                            break
                        itold = it
                        j += 1
                        dold = d
            # transpose data if rows should be plotted
            if args.transpose:
                data = data.T
            # plot versus integers
            if not args.xnum:
                datax = scalex*data[::,colx] + shiftx
            if scaley != 1. or shifty != 0.:
                    data[::,coly] = scaley*data[::,coly] + shifty
            if args.xnum:
                datax = np.arange(0,len(data[:,coly]),) + shiftx
            if args.abs: 
                data[::,coly] = abs(data[::,coly])
            # label
            lb = str(i)
            if i < len(args.coly): lb = coly - 1
            if i < len(args.legend): lb = args.legend[i].strip('{}')
            elif len(args.legend) > 0: lb = ''        
            if args.legend == [] and labels != []:
                lb = labels[coly]
            # enforce no legend
            if len(args.legend)>0 and args.legend[0] == 'none':
                lb = ''
            # initial definition of datay
            datay = data[::,coly]
            # actual plotting commands
            c = cl[color_index]
            mfc = 'white'
            if lis[i] == '.':
                mfc = c
            if args.logy:
                mask = abs(datay) > 1e-12
                ax.semilogy(datax[mask],datay[mask],
                                        lis[i],ms=msl[i],c=c,lw=lws[i],label=lb,mec=c,mfc=mfc,zorder=zorder)
            elif args.logx: 
                ax.semilogx(datax,datay,
                                        lis[i],ms=msl[i],c=c,label=lb,mec=c,mfc=mfc,zorder=zorder)
            elif args.log: 
                ax.loglog(datax,datay,
                                    lis[i],c=c,lw=lws[i],label=lb,mec=c,mfc=mfc,zorder=zorder)
            else: 
                ax.plot(datax,datay,
                        lis[i],ms=msl[i],c=c,lw=lws[i],label=lb,mec=c,mfc=mfc,zorder=zorder)
            # post process plots (set labels ...)
            if i == numberOfCurves-1:
                ax.set_xlabel(args.xlabel)
                ax.set_ylabel(args.ylabel)
                if len(args.xlim)>0:
                    pl.xlim(args.xlim[0],args.xlim[1])
                if len(args.ylim)>0:
                    pl.ylim(args.ylim[0],args.ylim[1])
                pimpaxis()
                postprocess(ax)
        if args.legendloc != 'none' or args.full:
             leg = axMain.legend(loc=args.legendloc,handlelength=args.handlelen,
                                 handletextpad=args.handletextpad,frameon=args.legendframe)
             frame = leg.get_frame()
             frame.set_color('white')

    def matrixplot():
            for filename in args.filenames:
                data = np.loadtxt(filename)
                if args.transpose:
                        data = data.T
                for i,row in enumerate(data):
                        data[i,i] = 0

                pl.matshow(data)

    if args.type == 'standard':
        standardplot()
    else:
        matrixplot()

    showorsavefig()
