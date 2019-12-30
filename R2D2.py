'''
R2D2 python module provides functions for reading R2D2 simulation results.
'''

# Global variables
p = {}
qc = {}
q2 = {}
q3 = {}
qi = {}
vc = {}

######################################################
def init(datadir):
    '''
    Thie function reads basic data for the calculation setting.
    The data is stored in R2D2.p dictionary

    Parameters:
        datadir (str): data location

    Returnes:
        None

    Examples:
        >>> import R2D2
        >>> datadir = 'data'
        >>> R2D2.read_init(dir)
        >>> print(R2D2.p['ix'])

    '''

    import numpy as np
    import sys

    
    p['datadir'] = datadir

    f = open(p['datadir']+"param/nd.dac","r")
    nn = f.read().split()
    nd = int(nn[0])
    ni = int(nn[1])
    f.close()
    
    p['nd'] = nd
    p['ni'] = ni
    
    R2D2_py_ver = 1.2
    f = open(p['datadir']+"param/params.dac","r")
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
            p[line.split()[1]] = int(line.split()[0])
        if line.split()[2] == 'd':
            p[line.split()[1]] = float(line.split()[0])
        line = f.readline()

    f.close()
            
    if p["swap"] == 0: # little
        p["endian"] = "<"
    else:
        p["endian"] = ">"

    p["ix"] = p["ix0"]*p["nx"]
    p["jx"] = p["jx0"]*p["ny"]
    p["kx"] = p["kx0"]*p["nz"]

    marginx = p["margin"]*(p["xdcheck"]-1)
    marginy = p["margin"]*(p["ydcheck"]-1)
    marginz = p["margin"]*(p["zdcheck"]-1)
   
    ixg = p["ix"] + 2*marginx
    jxg = p["jx"] + 2*marginy
    kxg = p["kx"] + 2*marginz
        
    endian = p["endian"]
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
    f = open(p['datadir']+"param/back.dac",'rb')
    back = np.fromfile(f,dtype=dtyp,count=1)
    f.close()

    for key in back.dtype.names:
        if back[key].size == ixg:
            p[key] = back[key].reshape((ixg),order="F")[marginx:ixg-marginx]
        elif back[key].size == jxg:
            p[key] = back[key].reshape((jxg),order="F")[marginy:jxg-marginy]
        elif back[key].size == kxg:
            p[key] = back[key].reshape((kxg),order="F")[marginz:kxg-marginz]

    p["rsun"] = 6.9598947e+10
    p["xr"] = p["x"]/p["rsun"]
    p["xn"] = (p["x"]-p["rsun"])*1.e-8
    
    if p['zdcheck'] == 2:
        dimension = '3d'
    else:
        dimension = '2d'

    ##############################
    # read value information
    if dimension == "3d":
        f = open(p['datadir']+"remap/vl/c.dac","r")
        value = f.read().split('\n')
        p["m2da"] = int(value[0])
        del value[0]
        p["cl"] = list(map(str.strip,value)) ## strip space from character
        f.close()

        ##############################
        # read mpi information
        npe = p["npe"]
        dtyp=np.dtype([ \
                ("iss",endian+str(npe)+"i4"),\
                    ("iee",endian+str(npe)+"i4"),\
                    ("jss",endian+str(npe)+"i4"),\
                    ("jee",endian+str(npe)+"i4"),\
                    ("iixl",endian+str(npe)+"i4"),\
                    ("jjxl",endian+str(npe)+"i4"),\
                    ("np_ijr",endian+str(p["ixr"]*p["jxr"])+"i4"),\
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
                p[key] = mpi[key].reshape((p["ixr"],p["jxr"]),order="F")
            else:
                p[key] = mpi[key].reshape((mpi[key].size),order="F")
    
        p["i2ir"] = p["i2ir"][marginx:ixg-marginx]
        p["j2jr"] = p["j2jr"][marginy:jxg-marginy]

        p["iss"] = p["iss"] - 1
        p["iee"] = p["iee"] - 1
        p["jss"] = p["jss"] - 1
        p["jee"] = p["jee"] - 1

######################################################
def read_qq_select(xs,n,silent=False,out=False):
    '''
    This funcition read 2D slice data at a selected height.
    The data is stored in R2D2.q2 dictionary

    Parameters:
        xs (float): a selected height for data
        n (int): a setected time step for data
        silent (logic): True suppresses a message of store
        out (logic): True returns stored data, otherwise stored only in R2D2.q2

    Returnes:
        None, but data dictionary is returned when out=False is specified

    '''

    import numpy as np
    i0 = np.argmin(np.abs(p["x"]-xs))
    ir0 = p["i2ir"][i0]
    mtype = p["mtype"]
    iixl = p["iixl"]
    jjxl = p["jjxl"]
    jx = p["kx"]
    kx = p["kx"]
    iss = p["iss"]
    jss = p["jss"]
    jee = p["jee"]

    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in q2:
        memflag = not q2['ro'].shape == (jx,kx)
    if 'ro' not in q2 or memflag:
        print('memory is newly allocated')
        q2["ro"] = np.zeros((jx,kx))
        q2["vx"] = np.zeros((jx,kx))
        q2["vy"] = np.zeros((jx,kx))
        q2["vz"] = np.zeros((jx,kx))
        q2["bx"] = np.zeros((jx,kx))
        q2["by"] = np.zeros((jx,kx))
        q2["bz"] = np.zeros((jx,kx))
        q2["se"] = np.zeros((jx,kx))
        q2["pr"] = np.zeros((jx,kx))
        q2["te"] = np.zeros((jx,kx))
        q2["op"] = np.zeros((jx,kx))
    for jr0 in range(1,p["jxr"]+1):
        np0 = p["np_ijr"][ir0-1,jr0-1]
        dtyp=np.dtype([ \
                ("qq",p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("pr",p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("te",p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("op",p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ])
        f = open(p['datadir']+"remap/qq/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
        qqq = np.fromfile(f,dtype=dtyp,count=1)
        q2["ro"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[0,i0-iss[np0],:,:]
        q2["vx"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[1,i0-iss[np0],:,:]
        q2["vy"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[2,i0-iss[np0],:,:]
        q2["vz"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[3,i0-iss[np0],:,:]
        q2["bx"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[4,i0-iss[np0],:,:]
        q2["by"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[5,i0-iss[np0],:,:]
        q2["bz"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[6,i0-iss[np0],:,:]
        q2["se"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[7,i0-iss[np0],:,:]
        q2["pr"][jss[np0]:jee[np0]+1,:] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        q2["te"][jss[np0]:jee[np0]+1,:] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        q2["op"][jss[np0]:jee[np0]+1,:] = qqq["op"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        f.close()

    if (not silent) and (not out) :
        print('### variales are stored in R2D2.q2 ###')
    
    if out:
        return q2

######################################################
######################################################
######################################################
### read horizontal variable
def read_qq(n,silent=False,out=False):
    '''
    This funcition read 3D full data.
    The data is stored in R2D2.q3 dictionary

    Parameters:
        n (int): a setected time step for data
        silent (logic): True suppresses a message of store
        out (logic): True returns stored data, otherwise stored only in R2D2.q3

    Returnes:
        None, but data dictionary is returned when out=False is specified

    '''

    import numpy as np
    
    mtype = p["mtype"]
    iixl = p["iixl"]
    jjxl = p["jjxl"]
    kx = p["kx"]
    iss = p["iss"]
    iee = p["iee"]
    jss = p["jss"]
    jee = p["jee"]
    ix = p["ix"]
    jx = p["jx"]
    kx = p["kx"]

    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in q3:
        memflag = not q3['ro'].shape == (ix,jx,kx)
    if 'ro' not in q3 or memflag:
        print('memory is newly allocated')
        q3["ro"] = np.zeros((ix,jx,kx))
        q3["vx"] = np.zeros((ix,jx,kx))
        q3["vy"] = np.zeros((ix,jx,kx))
        q3["vz"] = np.zeros((ix,jx,kx))
        q3["bx"] = np.zeros((ix,jx,kx))
        q3["by"] = np.zeros((ix,jx,kx))
        q3["bz"] = np.zeros((ix,jx,kx))
        q3["se"] = np.zeros((ix,jx,kx))
        q3["pr"] = np.zeros((ix,jx,kx))
        q3["te"] = np.zeros((ix,jx,kx))
        q3["op"] = np.zeros((ix,jx,kx))

    for ir0 in range(1,p["ixr"]+1):
        for jr0 in range(1,p["jxr"]+1):
            np0 = p["np_ijr"][ir0-1,jr0-1]
            dtyp=np.dtype([ \
                    ("qq",p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ("pr",p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ("te",p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ("op",p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ])
            f = open(p['datadir']+"remap/qq/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
            qqq = np.fromfile(f,dtype=dtyp,count=1)
            q3["ro"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[0,:,:,:]
            q3["vx"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[1,:,:,:]
            q3["vy"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[2,:,:,:]
            q3["vz"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[3,:,:,:]
            q3["bx"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[4,:,:,:]
            q3["by"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[5,:,:,:]
            q3["bz"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[6,:,:,:]
            q3["se"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[7,:,:,:]
            q3["pr"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            q3["te"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            q3["op"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["op"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            f.close()

    if (not silent) and (not out) :
        print('### variales are stored in R2D2.q3 ###')

    if out:
        return q3
    else:
        return
    
    return

######################################################
######################################################
### read all 2D variable
def read_qq_2d(n):
    '''
    This funcition read result of 2D simulation.

    Parameters:
        n (int): a setected time step for data

    Returnes:
        (dictionary): full 2D data

    '''
    import numpy as np

    ix = p["ix"]
    jx = p["jx"]

    qq1 = np.zeros((p["mtype"]+5,ix,jx))

    f = open(p['datadir']+"remap/qq.dac."+'{0:08d}'.format(n),"rb")
    qq0 = np.fromfile(f,p["endian"]+'f',(p["mtype"]+5)*p["ix"]*p["jx"])
    f.close()
    qq1 = np.reshape(qq0,(p["mtype"]+5,p["ix"],p["jx"]),order="F")

    qq = {"a":0}
    qq["ro"] = qq1[0,:,:]
    qq["vx"] = qq1[1,:,:]
    qq["vy"] = qq1[2,:,:]
    qq["vz"] = qq1[3,:,:]
    qq["bx"] = qq1[4,:,:]
    qq["by"] = qq1[5,:,:]
    qq["bz"] = qq1[6,:,:]
    qq["se"] = qq1[7,:,:]
    qq["en"] = qq1[8,:,:]
    qq["te"] = qq1[9,:,:]
    qq["pr"] = qq1[10,:,:]
    qq["tu"] = qq1[11,:,:]
    qq["fr"] = qq1[12,:,:]

    return qq
    
##############################
# read intensity related value
def read_tau(n,silent=False,out=False):
    '''
    This funcition read 2D data at certain optical depths.
    The data is stored in R2D2.qi dictionary.
    In this version the selected optical depth is 1, 0.1, and 0.01

    Parameters:
        n (int): a setected time step for data
        silent (logic): True suppresses a message of store
        out (logic): True returns stored data, otherwise stored only in R2D2.qi

    Returnes:
        None, but data dictionary is returned when out=False is specified

    '''
    import numpy as np

    f = open(p['datadir']+"tau/qq.dac."+'{0:08d}'.format(n),"rb")
    qq_in0 = np.fromfile(f,p["endian"]+'f',p["m_tu"]*p["m_in"]*p["jx"]*p["kx"])
    f.close()

    qi["in"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,0,:,:]
    qi["ro"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,1,:,:]
    qi["se"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,2,:,:]
    qi["pr"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,3,:,:]
    qi["te"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,4,:,:]
    qi["vx"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,5,:,:]
    qi["vy"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,6,:,:]
    qi["vz"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,7,:,:]
    qi["bx"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,8,:,:]
    qi["by"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,9,:,:]
    qi["bz"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,10,:,:]
    qi["he"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,11,:,:]
    qi["fr"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[0,12,:,:]

    qi["in01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,0,:,:]
    qi["ro01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,1,:,:]
    qi["se01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,2,:,:]
    qi["pr01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,3,:,:]
    qi["te01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,4,:,:]
    qi["vx01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,5,:,:]
    qi["vy01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,6,:,:]
    qi["vz01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,7,:,:]
    qi["bx01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,8,:,:]
    qi["by01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,9,:,:]
    qi["bz01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,10,:,:]
    qi["he01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,11,:,:]
    qi["fr01"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[1,12,:,:]

    qi["in001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,0,:,:]
    qi["ro001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,1,:,:]
    qi["se001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,2,:,:]
    qi["pr001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,3,:,:]
    qi["te001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,4,:,:]
    qi["vx001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,5,:,:]
    qi["vy001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,6,:,:]
    qi["vz001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,7,:,:]
    qi["bx001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,8,:,:]
    qi["by001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,9,:,:]
    qi["bz001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,10,:,:]
    qi["he001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,11,:,:]
    qi["fr001"] = np.reshape(qq_in0,(p["m_tu"],p["m_in"],p["jx"],p["kx"]),order="F")[2,12,:,:]

    if (not silent) and (not out):
        print('### variales are stored in R2D2.qi ###')

    if out:
        return qi

##############################
def read_time(n):
    '''
    This funcition read time at a selected time step

    Parameters:
        n (int): a setected time step for data

    Returnes:
        (float): time at a selected time step

    '''

    import numpy as np

    f = open(p['datadir']+"time/mhd/t.dac."+'{0:08d}'.format(n),"rb")
    t = np.fromfile(f,p['endian']+'d',1)
    f.close()    
    t = np.reshape(t,(1),order="F")[0]

    return t
    
##############################
# read remap_calc variable
def read_vc(n,silent=False,out=True):
    '''
    Thie function reads on the fly analysis data from fortran.
    The data is stored in R2D2.vc dictionary

    Parameters:
        dimension (str): 2D or 3D
        silent (logic): True suppresses a message of store
        out (logic): True returns stored data, otherwise stored only in R2D2.vc

    Returnes:
        None, but data dictionary is returned when out=False is specified

    '''
    import numpy as np

    f = open(p['datadir']+"remap/vl/vla.dac."+'{0:08d}'.format(n),"rb")
    vl0 = np.fromfile(f,p["endian"]+'f',p['m2da']*p['ix']*p['jx'])
    f.close()

    vl = np.reshape(vl0,(p['ix'],p['jx'],p['m2da']),order="F")

    for m in range(p["m2da"]):
        vc[p["cl"][m]] = vl[:,:,m]

    if (not silent) and (not out):
        print('### variales are stored in R2D2.vc ###')
    
    if out:
        return vc

    return

######################################################
######################################################
######################################################
### read full 3D variable for checkpoint

def read_qq_check(n,silent=False,out=False):
    '''
    Thie function reads 3D full data for checkpoint
    The data is stored in R2D2.qc dictionary

    Parameters:
        n (int): a setected time step for data
        silent (logic): True suppresses a message of store
        out (logic): True returns stored data, otherwise stored only in R2D2.qc

    Returnes:
        None, but data dictionary is returned when out=False is specified

    '''
    import numpy as np

    mtype = p['mtype']
    ix = p['ix']
    jx = p['jx']
    kx = p['kx']
    margin = p['margin']

    ixg = ix + 2*margin
    jxg = jx + 2*margin
    kxg = kx + 2*margin

    f = open(p['datadir']+"qq/qq.dac."+'{0:08d}'.format(n),'rb')
    qc = np.fromfile(f,p['endian']+'d',mtype*ixg*jxg*kxg).reshape((mtype,ixg,jxg,kxg),order="F")
    
    f.close()

    if (not silent) and (not out) :
        print('### variales are stored in R2D2.qc ###')
    
    if out:
        return qc

