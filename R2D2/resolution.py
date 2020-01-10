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
        ,endian='<'
        ,ixf=2,jxf=2,kxf=2
        ,ix_ununi=32,dx00=4.8e6,x_ununif=False):
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
    from scipy.interpolate import Rbf
    from scipy.interpolate import RegularGridInterpolator

    if not os.path.exists('../run/'+caseid):
        print('### Please download software to '+caseid+' first ###')
        return
    
    ixu = self.p['ix']*ixf
    jxu = self.p['jx']*jxf
    kxu = self.p['kx']*kxf

    ixug = ixu + 2*self.p['margin']
    jxug = jxu + 2*self.p['margin']
    kxug = kxu + 2*self.p['margin']

    import numpy as np
    if x_ununif:
        self.xu = gen_coord_ununiform(xmax,xmin,ixu,self.p['margin'],dx00,ix_ununi)
    else:
        self.xu = gen_coord(xmax,xmin,ixu,self.p['margin'])

    self.yu = gen_coord(ymax,ymin,jxu,self.p['margin'])
    self.zu = gen_coord(zmax,zmin,kxu,self.p['margin'])

    XU, YU, ZU = np.meshgrid(self.xu,self.yu,self.zu,indexing='ij')
        
    self.read_qq_check(n,silent=True)
    self.qu = np.zeros((self.p['mtype'],ixug,jxug,kxug))
    for m in range(0,self.p['mtype']):
        regrid_function = \
                RegularGridInterpolator((self.p['xg'],self.p['yg'],self.p['zg']) \
                                        ,self.qc[m,:,:,:],fill_value=0,bounds_error=False)
        tmp = regrid_function((XU,YU,ZU))
        #tmp[XU < self.p['xmin']] = 0
        #tmp[XU > self.p['xmax']] = 0
        #tmp[YU < self.p['ymin']] = 0
        #tmp[YU > self.p['ymax']] = 0
        #tmp[ZU < self.p['zmin']] = 0
        #tmp[ZU > self.p['zmax']] = 0
        self.qu[m,:,:,:] = tmp
        print(str(m)+' finished...')


    os.makedirs('../run/'+caseid+'/data/param/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/qq/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/remap/qq/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/remap/vl/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/time/mhd/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/time/tau/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/tau/',exist_ok=True)

    #f = open('../run/'+caseid+'/data/qq/qq.dac.e',mode='wb')
    #f.write(self.qu.reshape([self.p['mtype']*ixug*jxug*kxug],order='F'))    
    #f.close()
    #f = open('../run/'+caseid+'/data/time/mhd/t.dac.e',mode='wb')
    #f.write(
    #f.close()
    self.qu.reshape([self.p['mtype']*ixug*jxug*kxug] \
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
    
    print('nx0*ix0 = '+str(ixu),change_judge(ixu,self,'ix'))
    print('ny0*jx0 = '+str(jxu),change_judge(jxu,self,'jx'))
    print('nz0*kx0 = '+str(kxu),change_judge(kxu,self,'kx'))

    print(' ')

    if x_ununif:
        uniform_flag = '.true.'
    else:
        uniform_flag = '.false.'

    print('uniform_flag = '+uniform_flag \
          ,change_judge(x_ununif,self,'ununiform_flag'))

    print('ix_ununi = '+str(ix_ununi),change_judge(ix_ununi,self,'ix_ununi'))
    print('dx00 = ''{:.4e}'.format(dx00),change_judge(dx00,self,'dx00'))
        
