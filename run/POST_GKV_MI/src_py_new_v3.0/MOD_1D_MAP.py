from time import time as get_time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import LogFormatter
from matplotlib.ticker import LogFormatterExponent
from matplotlib.ticker import ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.animation import FuncAnimation
from netCDF4 import Dataset
from matplotlib.ticker import ScalarFormatter
from IPython.display import HTML
import json
import os
matplotlib.rcParams['animation.embed_limit'] = 2**128
import xarray as xr
import glob

from MOD_CONST import *
from MOD_MISC  import cmap_sym

writer = 'ffmpeg'
fps    = 15
dpi    = 200

fsize_ticks = 11
fsize_label = 13
fsize_legnd = 14

#################
def PLOT_1D( tit, dls, time, dfuncs, dlabel, dcolor, xlabel, ylabel, xmin, xmax, ymin, ymax, pdf_pages ):
#################

    plt.xlabel( xlabel, fontsize=fsize_label )
    plt.ylabel( ylabel, fontsize=fsize_label )

    plt.xticks(fontsize = fsize_ticks )
    plt.yticks(fontsize = fsize_ticks )

    if xmax>=0:
        plt.xlim( xmin, xmax )
    if ymax >=0:
        plt.ylim( ymin, ymax )

    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    plt.gca().yaxis.offsetText.set_fontsize(fsize_ticks)

    plt.title( tit )

    iend = len(dfuncs)

    for i in range(0,iend):
        plt.plot( time, dfuncs[i], label=dlabel[i], lw = 2, color=dcolor[i], ls=dls[i] )

    plt.legend( fontsize=fsize_legnd )

    pdf_pages.savefig()

    plt.close()



#=================
def PLOT_1D_MOVIE( DFIG, PROJ, f1, f2, time, tag1, tag2, tlabel, flabel, cadj ):
#=================

    print( f" # Plotting {flabel}: started ", "\n" )

    frames = len(time)

    fig, ax = plt.subplots()   #, gridspec_kw={'width_ratios': [1, 1]}, figsize=(20, 8))
    title_text = r'time: {:.2f} $[l/v_A]$'  # Initial title text with placeholder for time
    ##suptitle = fig.suptitle(title_text.format(0.0), fontsize=16)  # Initialize suptitle object
    tcnt = ax.set_title( tlabel + ', time: {:.2f} ($l/v_A$)'.format(time[0]))
    nzg = f1.shape[1]

    line1, = ax.plot([],[],'-', label=tag1, color="b")
    line2, = ax.plot([],[],'-', label=tag2, color="r")

    if cadj>0:
       fmin1, fmax1 = cmap_sym( f1[0] )
       fmin2, fmax2 = cmap_sym( f2[0] )
    else:
       fmin1, fmax1 = cmap_sym( f1 )
       fmin2, fmax2 = cmap_sym( f2 )
    fmax   = max(fmax1,fmax2)
    fmin   = min(fmin1,fmin2) 
    ax.set_ylim([fmin,fmax]) 
   

    xticks       = [0, nzg/4, nzg/2, 3*nzg/4, nzg-1]
    xticks_label = ["0","0.25","0.5", "0.75","1"]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks_label, fontsize=14)
    ax.set_xlabel( r'$z/l$' ,fontsize=14)

    ax.legend(fontsize=14,loc="upper right" )

    x = np.arange(0,nzg)
    def update( frame ):
        y1 = f1[frame]
        y2 = f2[frame]

        if cadj>0:
           fmin1, fmax1 = cmap_sym( f1[frame] )
           fmin2, fmax2 = cmap_sym( f2[frame] )
           fmax   = max(fmax1,fmax2)
           fmin   = min(fmin1,fmin2)
           ax.set_ylim([fmin,fmax])

        line1.set_data(x,y1)
        line2.set_data(x,y2)
        ##suptitle.set_text(title_text.format(time[frame]))
        tcnt.set_text( tlabel + ', time: {:.2f} ($l/v_A$)'.format(time[frame]))

    ani =  FuncAnimation(fig, update, frames=frames, interval=100)

    ofnm = DFIG + "fig_" + PROJ + f"_1D-zprof_{flabel}.mp4"
    ani.save(ofnm, writer=writer, fps=fps, dpi=dpi)

    plt.close()
