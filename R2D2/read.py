def init(self, datadir):
    '''
    This method reads basic data for calculation setting
    The data is stored in self.p dictionary

    Parameters:
        datadir (str): location of data
    '''
    import numpy as np
    import sys
        
    self.p = {}
    self.qs = {}
    self.qq = {}
    self.qt = {}
    self.t = 0
    self.vc = {}
        
    self.p['datadir'] = datadir 

    # read parameters
    f = open(self.p['datadir']+"param/nd.dac","r")
    nn = f.read().split()
    nd = int(nn[0])
    ni = int(nn[1])
    f.close()
    
    self.p['nd'] = nd
    self.p['ni'] = ni
    
    R2D2_py_ver = 1.2
    f = open(self.p['datadir']+"param/params.dac","r")
    line = f.readline().split()
    if R2D2_py_ver != float(line[2]):
        print("#######################################################")
        print("#######################################################")
        print("### Current R2D2 Python version is ",R2D2_py_ver,".")
        print("### You use the data from R2D2 version ",float(line[2]),".")
        print("### Please use the same version of fortran R2D2.")
        print("#######################################################")
        print("#######################################################")
        sys.exit()
        
    line = f.readline()
    while line:
        if line.split()[2] == 'i':
            self.p[line.split()[1]] = int(line.split()[0])
        if line.split()[2] == 'd':
            self.p[line.split()[1]] = float(line.split()[0])
        if line.split()[2] == 'c':
            self.p[line.split()[1]] = line.split()[0]
        if line.split()[2] == 'l':
            self.p[line.split()[1]] = line.split()[0]            
        line = f.readline()
    f.close()

    # transform endiant for python
    if self.p["swap"] == 0: # little
        self.p["endian"] = "<"
    else:
        self.p["endian"] = ">"

    self.p["ix"] = self.p["ix0"]*self.p["nx"]
    self.p["jx"] = self.p["jx0"]*self.p["ny"]
    self.p["kx"] = self.p["kx0"]*self.p["nz"]
   
    ixg = self.p["ix"] + 2*self.p["margin"]
    jxg = self.p["jx"] + 2*self.p["margin"]
    kxg = self.p["kx"] + 2*self.p["margin"]
        
    endian = self.p["endian"]
    dtyp=np.dtype([ \
                    ("head",endian+"i"),\
                    ("x",endian+str(ixg)+"d"),\
                    ("y",endian+str(jxg)+"d"),\
                    ("z",endian+str(kxg)+"d"),\
                    ("pr0",endian+str(ixg)+"d"),\
                    ("te0",endian+str(ixg)+"d"),\
                    ("ro0",endian+str(ixg)+"d"),\
                    ("se0",endian+str(ixg)+"d"),\
                    ("en0",endian+str(ixg)+"d"),\
                    ("op0",endian+str(ixg)+"d"),\
                    ("tu0",endian+str(ixg)+"d"),\
                    ("dsedr0",endian+str(ixg)+"d"),\
                    ("dtedr0",endian+str(ixg)+"d"),\
                    ("dprdro",endian+str(ixg)+"d"),\
                    ("dprdse",endian+str(ixg)+"d"),\
                    ("dtedro",endian+str(ixg)+"d"),\
                    ("dtedse",endian+str(ixg)+"d"),\
                    ("dendro",endian+str(ixg)+"d"),\
                    ("dendse",endian+str(ixg)+"d"),\
                    ("gx",endian+str(ixg)+"d"),\
                    ("kp",endian+str(ixg)+"d"),\
                    ("cp",endian+str(ixg)+"d"),\
                    ("fa",endian+str(ixg)+"d"),\
                    ("sa",endian+str(ixg)+"d"),\
                    ("xi",endian+str(ixg)+"d"),\
                    ("tail",endian+"i")\
    ])
    f = open(self.p['datadir']+"param/back.dac",'rb')
    back = np.fromfile(f,dtype=dtyp,count=1)
    f.close()

    self.p['xg'] = back['x'].reshape((ixg),order='F')
    self.p['yg'] = back['y'].reshape((jxg),order='F')
    self.p['zg'] = back['z'].reshape((kxg),order='F')
    
    for key in back.dtype.names:
        if back[key].size == ixg:
            self.p[key] = back[key].reshape((ixg),order="F")[self.p["margin"]:ixg-self.p["margin"]]
        elif back[key].size == jxg:
            self.p[key] = back[key].reshape((jxg),order="F")[self.p["margin"]:jxg-self.p["margin"]]
        elif back[key].size == kxg:
            self.p[key] = back[key].reshape((kxg),order="F")[self.p["margin"]:kxg-self.p["margin"]]
            
    self.p["rsun"] = 6.9598947e+10
    self.p["xr"] = self.p["x"]/self.p["rsun"]
    self.p["xn"] = (self.p["x"]-self.p["rsun"])*1.e-8
    
    if self.p['zdcheck'] == 2:
        dimension = '3d'
    else:
        dimension = '2d'
        
    ##############################
    # read value information
    if dimension == "3d":
        f = open(self.p['datadir']+"remap/vl/c.dac","r")
        value = f.read().split('\n')
        self.p["m2da"] = int(value[0])
        del value[0]
        self.p["cl"] = list(map(str.strip,value)) ## strip space from character
        f.close()

        ##############################
        # read mpi information
        npe = self.p["npe"]
        dtyp=np.dtype([ \
                    ("iss",endian+str(npe)+"i4"),\
                    ("iee",endian+str(npe)+"i4"),\
                    ("jss",endian+str(npe)+"i4"),\
                    ("jee",endian+str(npe)+"i4"),\
                    ("iixl",endian+str(npe)+"i4"),\
                    ("jjxl",endian+str(npe)+"i4"),\
                    ("np_ijr",endian+str(self.p["ixr"]*self.p["jxr"])+"i4"),\
                    ("ir",endian+str(npe)+"i4"),\
                    ("jr",endian+str(npe)+"i4"),\
                    ("i2ir",endian+str(ixg)+"i4"),\
                    ("j2jr",endian+str(jxg)+"i4"),\
        ])
        
        f = open(datadir+"remap/remap_info.dac",'rb')
        mpi = np.fromfile(f,dtype=dtyp,count=1)    
        f.close()
        
        for key in mpi.dtype.names:
            if key == "np_ijr":
                self.p[key] = mpi[key].reshape((self.p["ixr"],self.p["jxr"]),order="F")
            else:
                self.p[key] = mpi[key].reshape((mpi[key].size),order="F")
    
        self.p["i2ir"] = self.p["i2ir"][self.p["margin"]:ixg-self.p["margin"]]
        self.p["j2jr"] = self.p["j2jr"][self.p["margin"]:jxg-self.p["margin"]]
        
        self.p["iss"] = self.p["iss"] - 1
        self.p["iee"] = self.p["iee"] - 1
        self.p["jss"] = self.p["jss"] - 1
        self.p["jee"] = self.p["jee"] - 1

##############################
def read_qq_select(self,xs,n,silent=False):
    '''
    This method reads 2D slice data at a selected height.
    The data is stored in self.qs dictionary

    Parameters:
        xs (float): a selected height for data
        n (int): a setected time step for data
        silent (logic): True suppresses a message of store
        
    '''

    import numpy as np
    i0 = np.argmin(np.abs(self.p["x"]-xs))
    ir0 = self.p["i2ir"][i0]
    mtype = self.p["mtype"]
    iixl = self.p["iixl"]
    jjxl = self.p["jjxl"]
    jx = self.p["kx"]
    kx = self.p["kx"]
    iss = self.p["iss"]
    jss = self.p["jss"]
    jee = self.p["jee"]
    
    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in self.qs:
        memflag = not self.qs['ro'].shape == (jx,kx)
    if 'ro' not in self.qs or memflag:
        print('memory is newly allocated')
        self.qs["ro"] = np.zeros((jx,kx))
        self.qs["vx"] = np.zeros((jx,kx))
        self.qs["vy"] = np.zeros((jx,kx))
        self.qs["vz"] = np.zeros((jx,kx))
        self.qs["bx"] = np.zeros((jx,kx))
        self.qs["by"] = np.zeros((jx,kx))
        self.qs["bz"] = np.zeros((jx,kx))
        self.qs["se"] = np.zeros((jx,kx))
        self.qs["pr"] = np.zeros((jx,kx))
        self.qs["te"] = np.zeros((jx,kx))
        self.qs["op"] = np.zeros((jx,kx))
    for jr0 in range(1,self.p["jxr"]+1):
        np0 = self.p["np_ijr"][ir0-1,jr0-1]
        dtyp=np.dtype([ \
                ("qq",self.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("pr",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("te",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("op",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
        ])
        f = open(self.p['datadir']+"remap/qq/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
        qqq = np.fromfile(f,dtype=dtyp,count=1)
        self.qs["ro"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[0,i0-iss[np0],:,:]
        self.qs["vx"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[1,i0-iss[np0],:,:]
        self.qs["vy"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[2,i0-iss[np0],:,:]
        self.qs["vz"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[3,i0-iss[np0],:,:]
        self.qs["bx"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[4,i0-iss[np0],:,:]
        self.qs["by"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[5,i0-iss[np0],:,:]
        self.qs["bz"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[6,i0-iss[np0],:,:]
        self.qs["se"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[7,i0-iss[np0],:,:]
        self.qs["pr"][jss[np0]:jee[np0]+1,:] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        self.qs["te"][jss[np0]:jee[np0]+1,:] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        self.qs["op"][jss[np0]:jee[np0]+1,:] = qqq["op"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        f.close()
        
    if not silent :
        print('### variales are stored in self.qs ###')
            
##############################
def read_qq(self,n,silent=False):
    '''
    This method reads 3D full data
    The data is stored in self.qq dictionary

    Parameters:
        n (int): a selected time step for data
        silent (logic): True suppresses a message of store
    '''
    import numpy as np
    
    mtype = self.p["mtype"]
    iixl = self.p["iixl"]
    jjxl = self.p["jjxl"]
    kx = self.p["kx"]
    iss = self.p["iss"]
    iee = self.p["iee"]
    jss = self.p["jss"]
    jee = self.p["jee"]
    ix = self.p["ix"]
    jx = self.p["jx"]
    kx = self.p["kx"]

    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in self.qq:
        memflag = not self.qq['ro'].shape == (ix,jx,kx)
    if 'ro' not in self.qq or memflag:
        print('memory is newly allocated')
        self.qq["ro"] = np.zeros((ix,jx,kx))
        self.qq["vx"] = np.zeros((ix,jx,kx))
        self.qq["vy"] = np.zeros((ix,jx,kx))
        self.qq["vz"] = np.zeros((ix,jx,kx))
        self.qq["bx"] = np.zeros((ix,jx,kx))
        self.qq["by"] = np.zeros((ix,jx,kx))
        self.qq["bz"] = np.zeros((ix,jx,kx))
        self.qq["se"] = np.zeros((ix,jx,kx))
        self.qq["pr"] = np.zeros((ix,jx,kx))
        self.qq["te"] = np.zeros((ix,jx,kx))
        self.qq["op"] = np.zeros((ix,jx,kx))

    for ir0 in range(1,self.p["ixr"]+1):
        for jr0 in range(1,self.p["jxr"]+1):
            np0 = self.p["np_ijr"][ir0-1,jr0-1]
            dtyp=np.dtype([ \
                ("qq",self.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("pr",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("te",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("op",self.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
            ])
            f = open(self.p['datadir']+"remap/qq/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
            qqq = np.fromfile(f,dtype=dtyp,count=1)
            self.qq["ro"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[0,:,:,:]
            self.qq["vx"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[1,:,:,:]
            self.qq["vy"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[2,:,:,:]
            self.qq["vz"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[3,:,:,:]
            self.qq["bx"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[4,:,:,:]
            self.qq["by"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[5,:,:,:]
            self.qq["bz"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[6,:,:,:]
            self.qq["se"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[7,:,:,:]
            self.qq["pr"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            self.qq["te"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            self.qq["op"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["op"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            f.close()

    if not silent :
        print('### variales are stored in self.qq ###')


##############################
def read_qq_tau(self,n,silent=False):
    '''
    This method reads 2D data at certain optical depths.
    The data is stored in self.qt dictionary.
    In this version the selected optical depth is 1, 0.1, and 0.01

    Parameters:
        n (int): a setected time step for data
        silent (logic): True suppresses a message of store
    '''
    import numpy as np

    f = open(self.p['datadir']+"tau/qq.dac."+'{0:08d}'.format(n),"rb")
    qq_in0 = np.fromfile(f,self.p["endian"]+'f',self.p["m_tu"]*self.p["m_in"]*self.p["jx"]*self.p["kx"])
    f.close()

    m_tu = self.p['m_tu']
    m_in = self.p['m_in']
    jx = self.p['jx']
    kx = self.p['kx']
        
    self.qt["in"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,0,:,:]
    self.qt["ro"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,1,:,:]
    self.qt["se"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,2,:,:]
    self.qt["pr"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,3,:,:]
    self.qt["te"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,4,:,:]
    self.qt["vx"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,5,:,:]
    self.qt["vy"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,6,:,:]
    self.qt["vz"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,7,:,:]
    self.qt["bx"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,8,:,:]
    self.qt["by"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,9,:,:]
    self.qt["bz"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,10,:,:]
    self.qt["he"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,11,:,:]
    self.qt["fr"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[0,12,:,:]

    self.qt["in01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,0,:,:]
    self.qt["ro01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,1,:,:]
    self.qt["se01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,2,:,:]
    self.qt["pr01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,3,:,:]
    self.qt["te01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,4,:,:]
    self.qt["vx01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,5,:,:]
    self.qt["vy01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,6,:,:]
    self.qt["vz01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,7,:,:]
    self.qt["bx01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,8,:,:]
    self.qt["by01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,9,:,:]
    self.qt["bz01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,10,:,:]
    self.qt["he01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,11,:,:]
    self.qt["fr01"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[1,12,:,:]

    self.qt["in001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,0,:,:]
    self.qt["ro001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,1,:,:]
    self.qt["se001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,2,:,:]
    self.qt["pr001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,3,:,:]
    self.qt["te001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,4,:,:]
    self.qt["vx001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,5,:,:]
    self.qt["vy001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,6,:,:]
    self.qt["vz001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,7,:,:]
    self.qt["bx001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,8,:,:]
    self.qt["by001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,9,:,:]
    self.qt["bz001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,10,:,:]
    self.qt["he001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,11,:,:]
    self.qt["fr001"] = np.reshape(qq_in0,(m_tu,m_in,jx,kx),order="F")[2,12,:,:]

    if not silent :
        print('### variales are stored in self.qt ###')
                
##############################
def read_time(self,n,tau=False,silent=True):
    '''
    This method reads time at a selected time step
    The data is stored in self.t

    Parameters:
        n (int): a setected time step for data
        tau (logic): if True time for optical depth

    '''

    import numpy as np

    if tau:
        f = open(self.p['datadir']+"time/tau/t.dac."+'{0:08d}'.format(n),"rb")
        self.t = np.fromfile(f,self.p['endian']+'d',1).reshape((1),order='F')[0]            
        f.close()    
    else:
        f = open(self.p['datadir']+"time/mhd/t.dac."+'{0:08d}'.format(n),"rb")
        self.t = np.fromfile(f,self.p['endian']+'d',1).reshape((1),order='F')[0]            
        f.close()    
    
    if not silent :
        print('### variales are stored in self.t ###')

    return self.t
##############################
def read_vc(self,n,silent=False):
    '''
    This method reads on the fly analysis data from fortran.
    The data is stored in self.vc dictionary

    Parameters:
        silent (logic): True suppresses a message of store
    '''

    import numpy as np

    f = open(self.p['datadir']+"remap/vl/vla.dac."+'{0:08d}'.format(n),"rb")
    vl = np.fromfile(f,self.p["endian"]+'f',self.p['m2da']*self.p['ix']*self.p['jx']) \
           .reshape((self.p['ix'],self.p['jx'],self.p['m2da']),order="F")
    f.close()

    for m in range(self.p["m2da"]):
        self.vc[self.p["cl"][m]] = vl[:,:,m]

    if not silent :
        print('### variales are stored in self.vc ###')

##############################
def read_qq_check(self,n,silent=False):
    '''
    This method reads 3D full data for checkpoint
    The data is stored in self.qc dictionary

    Parameters:
        n (int): a setected time step for data
        silent (logic): True suppresses a message of store

    '''

    import numpy as np

    mtype = self.p['mtype']
    ix = self.p['ix']
    jx = self.p['jx']
    kx = self.p['kx']
    margin = self.p['margin']

    ixg = ix + 2*margin
    jxg = jx + 2*margin
    kxg = kx + 2*margin
    
    f = open(self.p['datadir']+"qq/qq.dac."+'{0:08d}'.format(n),'rb')
    self.qc = np.fromfile(f,self.p['endian']+'d',mtype*ixg*jxg*kxg).reshape((mtype,ixg,jxg,kxg),order="F")    
    f.close()
    
    if not silent :
        print('### variales are stored in self.qc ###')
