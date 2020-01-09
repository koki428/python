from . import read
from . import google
from . import resolution
class R2D2_data():
    '''
    R2D2_data is a class for reading R2D2 data

    Attributes:
        p (dic): basic parameters
        qs (dic): 2D data at selected height
        qq (dic): 3D full data 
        qt (dic): 2D data at constant optical depths
        t (float): time
        vc (dic): data of on the fly analysis
    '''

R2D2_data.__init__       = read.init        
R2D2_data.read_qq_select = read.read_qq_select
R2D2_data.read_qq        = read.read_qq
R2D2_data.read_qq_tau    = read.read_qq_tau
R2D2_data.read_time      = read.read_time
R2D2_data.read_vc        = read.read_vc
R2D2_data.read_qq_check  = read.read_qq_check

R2D2_data.out_gspread    = google.out_gspread
R2D2_data.upgrade_resolution = resolution.upgrade_resolution

