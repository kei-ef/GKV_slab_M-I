#====== Libraries ======
import sys
#ys.path.append("/home/fujitak/M-I_MODEL/POST_GKV_MI/py_src_dipole_v1")
sys.path.append("./src_py_new_v3.0")
from time import time as get_time
import numpy as np
import json
import os

from MOD_DIRS import *
from MOD_CONST import *
from MOD_PARAM import *
from MOD_HST  import MAIN_HST
from MOD_IONO import MAIN_IONO
from MOD_PHI  import MAIN_PHI
from MOD_MOM  import MAIN_MOM
from MOD_TRN  import MAIN_TRN
from MOD_CNSV import MAIN_CNSV
from MOD_FZV  import MAIN_FZV
from MOD_EFLX import MAIN_EFLX

#====== Flags ======

# 1: read NetCDF; 2: read saved file

ftag      = ""
flag_hst  = 0
flag_iono = 0
flag_phi  = 0
flag_mom  = 0
flag_trn  = 0
flag_cnsv = 0
flag_fzv  = 2
flag_eflx = 0

ist, iend  = 0, -1 # time step 
tst, tend  = 0, -1 # time in MHD unit
iskip = 1
tave_int, tave_end = 25.0, -1 # in the MHD unit, <0 for auto
jst, jend  = -1, -1
xmin, xmax = 0, -1
ymin, ymax = -1, -1
fmin, fmax = 0, 0.5 # down and up energy flux 
kappa   = 0.0
fct_min = 1e-3
ispc = 0
ipos = 2
pos0 = 0.5
ns = 2
mx_plt, my_plt = -2, 1
t4x4   = [ 6.0, 15.6, 16.5, 18 ]

plt_2D_x = 1
plt_2D_k = 0
plt_1D_x = 0
plt_1D_k = 0
plt_snap = 1
plt_spct = 0
cadj     = 1 # color bar adjustment; 0:fix; 1:dynamic

dict_plt = { "ist"      : ist,
             "iend"     : iend,
             "iskip"    : iskip,
             "jst"      : jst,
             "jend"     : jend,
             "tst"      : tst,
             "tend"     : tend,
             "t4x4"     : t4x4,
             "xmin"     : xmin,
             "xmax"     : xmax,
             "ymin"     : ymin,
             "ymax"     : ymax,
             "fmin"     : fmin,
             "fmax"     : fmax,
             "kappa"    : kappa,
             "fct_min"  : fct_min,
             "tave_int" : tave_int,
             "tave_end" : tave_end,
             "mx_plt"   : mx_plt,
             "my_plt"   : my_plt,
             "ispc"     : ispc,
             "ipos"     : ipos,
             "pos0"     : pos0,
             "ftag"     : ftag,
             "plt_2D_x" : plt_2D_x,
             "plt_2D_k" : plt_2D_k,
             "plt_1D_x" : plt_1D_x,
             "plt_1D_k" : plt_1D_k,
             "plt_snap" : plt_snap,
             "plt_spct" : plt_spct,
             "cadj"     : cadj }




#====== Start ======

print('*** Post process (dipole ver.1) ***', '\n' )

#------ Set directories and file names ------

CDIR = os.getcwd()
HEAD = "gkvp-slabMI_"

dict_dir = SET_DIRS(HEAD,CDIR)

DCDF = dict_dir["DCDF"]
DFIG = dict_dir["DFIG"]
PROJ = dict_dir["PROJ"]
DHST = dict_dir["DHST"]
DOUT = dict_dir["DOUT"]

#------ Construct normalization constants ------

dict_resc = Construct_rescaling_factors()

####################

if flag_hst>0:
   MAIN_HST( flag_hst, dict_plt, dict_dir, dict_resc )

if flag_iono>0:
   MAIN_IONO( flag_iono, dict_plt, dict_dir, dict_resc, ist, iend, iskip )

if flag_phi>0:
   MAIN_PHI( flag_phi, dict_plt, dict_dir, dict_resc, ist, iend, iskip )

if flag_mom>0:
   MAIN_MOM( flag_mom, dict_plt, dict_dir, dict_resc, ist, iend, iskip, ns )

if flag_trn>0:
   MAIN_TRN( flag_trn, dict_plt, dict_dir, dict_resc, ist, iend, iskip, ns )

if flag_cnsv>0:
   MAIN_CNSV( flag_cnsv, dict_plt, dict_dir, dict_resc, ist, iend, iskip, ns )

if flag_fzv>0:
   MAIN_FZV( flag_fzv, dict_plt, dict_dir, dict_resc, ist, iend, iskip )

if flag_eflx>0:
   MAIN_EFLX( flag_eflx, dict_plt, dict_dir, dict_resc, ist, iend, iskip )
