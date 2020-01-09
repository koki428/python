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
        ,x_ununif=False):
    '''
    '''

    ixu = self.p['ix']*ixf
    jxu = self.p['jx']*jxf
    kxu = self.p['kx']*kxf

    ixug = ixu + 2*margin
    jxug = jxu + 2*margin
    kxug = kxu + 2*margin
        
    import numpy as np
    if x_ununif:
        xu = gen_coord(xmax,xmin,ixu,self.p['margin'])        
    else:
        xu = gen_coord(xmax,xmin,ixu,self.p['margin'])

    yu = gen_coord(ymax,ymin,jxu,self.p['margin'])
    zu = gen_coord(zmax,zmin,kxu,self.p['kx']*kxf,p['margin'])

    XU, YU, ZU = np.meshgrid(xu,yu,zu,indexing='ij',sparse=True)
    
    qu = np.zeros((mtype,))
