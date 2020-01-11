######################################################
######################################################
def gen_coord(xmax,xmin,ix,margin):
    '''
    This fucntion defines uniform geometry

    Parameters:
        xmax (float): location of upper boundary
        xmin (float): location of lower boundary
        ix (int): number of grid without margin
        margin (int): number of margin

    Return:
        x (float) [ix + 2*margin]: generated geometry

    '''
    import numpy as np
    dx = (xmax - xmin)/ix
    x = np.arange(xmin - (margin - 0.5)*dx,xmax + (margin + 0.5)*dx,dx)

    return x

######################################################
######################################################
def gen_coord_ununiform(xmax,xmin,ix,margin,dx00,ix_ununi):
    '''
    This function defines ununiform geometry

    Parameters: 
        xmax (float): location of upper boundary
        xmin (float): location of lower boundary
        ix (int): number of grid without margin
        margin (int): number of margin
        dx00 (float): grid spacing in uniform grid region
        ix_uniuni (int): number of uniform grid

    Return:
        x (float) [ix + 2*margin]: generated geometry
   
    '''
    import numpy as np
    
    xrange = xmax - xmin
    xrange0 = dx00*ix_ununi
    xrange1 = xrange - xrange0
    nxx = ix - ix_ununi

    ixg = ix + 2*margin

    fdx = 2*(xrange1 - dx00*nxx)/(nxx-4)/nxx
    x = np.zeros(ixg)
    x[ixg - margin - 1] = xmax - 0.5*dx00
    for i in range(ixg - margin,ixg):
        x[i] = x[i-1] + dx00

    for i in range(ixg - margin - 2, ixg - margin - ix_ununi - 3, -1):
        x[i] = x[i+1] - dx00

    dx11 = dx00
    for i in range(ixg - margin - ix_ununi-3,3,-1):
        x[i] = x[i+1] - dx11
        dx11 = dx11 + fdx

    for i in range(3,-1,-1):
        x[i] = x[i+1] - dx11
    
    return x

######################################################
######################################################
def upgrade_resolution(
        self,caseid,n
        ,xmin,xmax,ymin,ymax,zmin,zmax 
        ,ixf=2,jxf=2,kxf=2
        ,ix_ununi=32,dx00=4.8e6,x_ununif=False
        ,endian='<',end_step=False
        ):
    '''
    This function chabges the resolution and output for the next
    calculation.

    Parameters:
        self (R2D2_data): instance of R2D2_data
        caseid (str): caseid of destination directory e.g. 'd002'
        xmax (float): max location in x direction
        xmin (float): min location in x direction
        ymax (float): max location in y direction
        ymin (float): min location in y direction
        zmax (float): max location in z direction
        zmin (float): min location in z direction

        endian (str): endian
    
        ixf (int): increase factor of grid point in x direction
        jxf (int): increase factor of grid point in y direction
        kxf (int): increase factor of grid point in z direction
        
        below prameters are effective only when x_ununif=True
        ix_ununi (int): number of grid in uniform grid region
        dx00 (float): grid spacing in uniform grid region
        x_ununif (bool): whethere ununiform grid is used
    
    '''
    import os
    import os.path
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator
    from . import common

    self.up = {}

    ## read Model S based stratification
    ## data is stored in self.models dictionary
    self.models_init() 

    ## check if the destination directory
    ## if not this method finishes
    if not os.path.exists('../run/'+caseid):
        print('### Please download software to '+caseid+' first ###')
        return

    ## number of grid after upgrade
    self.up['ix'] = self.p['ix']*ixf
    self.up['jx'] = self.p['jx']*jxf
    self.up['kx'] = self.p['kx']*kxf

    ## number of grid with margin after upgrade
    self.up['ixg'] = self.up['ix'] + 2*self.p['margin']
    self.up['jxg'] = self.up['jx'] + 2*self.p['margin']
    self.up['kxg'] = self.up['kx'] + 2*self.p['margin']

    if x_ununif:
        ## generate upgraded coordinate in ununiform geometry
        self.up['x'] = gen_coord_ununiform(xmax,xmin,self.up['ix'],self.p['margin'],dx00,ix_ununi)
    else:
        ## generate upgraded coordinate in uniform geometry
        self.up['x'] = gen_coord(xmax,xmin,self.up['ix'],self.p['margin'])

    ## generate upgraded coordinate
    self.up['y'] = gen_coord(ymax,ymin,self.up['jx'],self.p['margin'])
    self.up['z'] = gen_coord(zmax,zmin,self.up['kx'],self.p['margin'])

    ## background density and entropy in original setting
    RO0, tmp, tmp = np.meshgrid(self.p['ro0g'],self.p['yg'],self.p['zg'],indexing='ij')
    SE0, tmp, tmp = np.meshgrid(self.p['se0g'],self.p['yg'],self.p['zg'],indexing='ij')

    ## background density and entropy in upgraded setting
    self.up['ro0'] = np.interp(self.up['x'],self.models['x']*self.p['rsun'],self.models['ro0'])
    self.up['se0'] = np.interp(self.up['x'],self.models['x']*self.p['rsun'],self.models['se0'])

    ## 1D arrays are converted to 3D
    XU, YU, ZU = np.meshgrid(self.up['x'],self.up['y'],self.up['z'],indexing='ij')
    RO0U, tmp, tmp = np.meshgrid(self.up['ro0'],self.up['y'],self.up['z'],indexing='ij')
    SE0U, tmp, tmp = np.meshgrid(self.up['se0'],self.up['y'],self.up['z'],indexing='ij')

    ## read original checkpoint data
    self.read_qq_check(n,silent=True,end_step=end_step)

    ## generate upgraded checkpoint data
    self.qu = np.zeros((self.p['mtype'],self.up['ixg'],self.up['jxg'],self.up['kxg']))
    rob = np.zeros((self.up['ixg'],self.up['jxg'],self.up['kxg']))
    seb = np.zeros((self.up['ixg'],self.up['jxg'],self.up['kxg']))
    
    for m in [0,7]:
        regrid_function = \
                RegularGridInterpolator((self.p['xg'],self.p['yg'],self.p['zg']) \
                                        ,self.qc[m,:,:,:],fill_value=0,bounds_error=False)
        tmp = regrid_function((XU,YU,ZU))
        if m == 0:
            rob = tmp
        if m == 7:
            seb = tmp
    
    self.qc[0,:,:,:] = np.log(self.qc[0,:,:,:] + RO0)
    self.qc[7,:,:,:] = self.qc[7,:,:,:] + SE0
    
    for m in range(0,self.p['mtype']):
        regrid_function = \
                RegularGridInterpolator((self.p['xg'],self.p['yg'],self.p['zg']) \
                                        ,self.qc[m,:,:,:],fill_value=0,bounds_error=False)
        tmp = regrid_function((XU,YU,ZU))
        self.qu[m,:,:,:] = tmp
        print(str(m)+' finished...')

    ros = np.exp(self.qu[0,:,:,:]) - RO0U
    ses = self.qu[7,:,:,:] - SE0U

    ros[np.abs(rob)/RO0U < 1.e-3] = rob[np.abs(rob)/RO0U < 1.e-3]
    ses[np.abs(seb)/SE0U < 1.e-3] = seb[np.abs(seb)/SE0U < 1.e-3]

    self.qu[0,:,:,:] = ros
    self.qu[7,:,:,:] = ses    

    os.makedirs('../run/'+caseid+'/data/param/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/qq/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/remap/qq/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/remap/vl/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/time/mhd/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/time/tau/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/tau/',exist_ok=True)

    self.qu.reshape([self.p['mtype']*self.up['ixg']*self.up['jxg']*self.up['kxg']] \
            ,order='F').astype(endian+'d').tofile('../run/'+caseid+'/data/qq/qq.dac.e')

    t = np.zeros(1)
    t.reshape([1] \
            ,order='F').astype(endian+'d').tofile('../run/'+caseid+'/data/time/mhd/t.dac.e')

    f = open('../run/'+caseid+'/data/param/nd.dac',mode='w')
    f.write('{0:08d}'.format(0)+'{0:08d}'.format(0))
    f.close()

    print(' ')
    print('### Data upgrade funished ###')
    print('Please use following parameters')
    print(' ')

    def sign_judge(value):
        if np.sign(value) == 1.0:
            sign = "+"
        else:
            sign = "-"
        return sign

    def change_judge(value,self,key):
        if value == self.p[key]:
            change = ' (unchange)'
        else:
            change = ' \033[31m(changed)\033[0m'
        return change
    
    value = xmin - self.p['rsun']
    print('xmin = rsun '+sign_judge(value),'{:.4e}'.format(abs(value)),change_judge(xmin,self,'xmin'))
    value = xmax - self.p['rsun']
    print('xmax = rsun '+sign_judge(value),'{:.4e}'.format(abs(value)),change_judge(xmax,self,'xmax'))
    print('ymin = '+'{:.4e}'.format(ymin),change_judge(ymin,self,'ymin'))
    print('ymax = '+'{:.4e}'.format(ymax),change_judge(ymax,self,'ymax'))
    print('zmin = '+'{:.4e}'.format(zmin),change_judge(zmin,self,'ymin'))
    print('zmax = '+'{:.4e}'.format(zmax),change_judge(zmax,self,'ymax'))
    print(' ')
    
    print('nx0*ix0 = '+str(self.up['ix']),change_judge(self.up['ix'],self,'ix'))
    print('ny0*jx0 = '+str(self.up['jx']),change_judge(self.up['jx'],self,'jx'))
    print('nz0*kx0 = '+str(self.up['kx']),change_judge(self.up['kx'],self,'kx'))

    print(' ')

    if x_ununif:
        uniform_flag = '.true.'
    else:
        uniform_flag = '.false.'

    print('uniform_flag = '+uniform_flag \
          ,change_judge(x_ununif,self,'ununiform_flag'))

    print('ix_ununi = '+str(ix_ununi),change_judge(ix_ununi,self,'ix_ununi'))
    print('dx00 = ''{:.4e}'.format(dx00),change_judge(dx00,self,'dx00'))
        
