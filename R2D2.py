# R2D2.py
######################################################
######################################################
######################################################
def read_init(dir,dimension):
    import numpy as np
    import config as c
    import sys
    c.p = {"a":0}
    f = open(dir+"param/nd.dac","r")
    nn = f.read().split()
    nd = int(nn[0])
    ni = int(nn[1])
    f.close()
    
    c.p['nd'] = nd
    c.p['ni'] = ni
    
    R2D2_py_ver = 1.1
    f = open(dir+"param/params.dac","r")
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
            c.p[line.split()[1]] = int(line.split()[0])
        if line.split()[2] == 'd':
            c.p[line.split()[1]] = float(line.split()[0])
        line = f.readline()

    f.close()
            
    if c.p["swap"] == 0: # little
        c.p["endian"] = "<"
    else:
        c.p["endian"] = ">"

    c.p["ix"] = c.p["ix0"]*c.p["nx"]
    c.p["jx"] = c.p["jx0"]*c.p["ny"]
    c.p["kx"] = c.p["kx0"]*c.p["nz"]

    marginx = c.p["margin"]*(c.p["xdcheck"]-1)
    marginy = c.p["margin"]*(c.p["ydcheck"]-1)
    marginz = c.p["margin"]*(c.p["zdcheck"]-1)
   
    ixg = c.p["ix"] + 2*marginx
    jxg = c.p["jx"] + 2*marginy
    kxg = c.p["kx"] + 2*marginz
        
    endian = c.p["endian"]
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
    f = open(dir+"back.dac",'rb')
    back = np.fromfile(f,dtype=dtyp,count=1)
    f.close()

    for key in back.dtype.names:
        if back[key].size == ixg:
            c.p[key] = back[key].reshape((ixg),order="F")[marginx:ixg-marginx]
        elif back[key].size == jxg:
            c.p[key] = back[key].reshape((jxg),order="F")[marginy:jxg-marginy]
        elif back[key].size == kxg:
            c.p[key] = back[key].reshape((kxg),order="F")[marginz:kxg-marginz]

    c.p["rsun"] = 6.9598947e+10
    c.p["xr"] = c.p["x"]/c.p["rsun"]
    c.p["xn"] = (c.p["x"]-c.p["rsun"])*1.e-8
    
    ##############################
    # read value information
    if dimension == "3d":
        f = open(dir+"remap/c.dac","r")
        value = f.read().split('\n')
        c.p["m2da"] = int(value[0])
        del value[0]
        c.p["cl"] = list(map(str.strip,value)) ## strip space from character
        f.close()

        ##############################
        # read mpi information
        npe = c.p["npe"]
        dtyp=np.dtype([ \
                ("iss",endian+str(npe)+"i4"),\
                    ("iee",endian+str(npe)+"i4"),\
                    ("jss",endian+str(npe)+"i4"),\
                    ("jee",endian+str(npe)+"i4"),\
                    ("iixl",endian+str(npe)+"i4"),\
                    ("jjxl",endian+str(npe)+"i4"),\
                    ("np_ijr",endian+str(c.p["ixr"]*c.p["jxr"])+"i4"),\
                    ("ir",endian+str(npe)+"i4"),\
                    ("jr",endian+str(npe)+"i4"),\
                    ("i2ir",endian+str(ixg)+"i4"),\
                    ("j2jr",endian+str(jxg)+"i4"),\
                    ])
        
        f = open(dir+"remap/remap_info.dac",'rb')
        mpi = np.fromfile(f,dtype=dtyp,count=1)    
        f.close()

        for key in mpi.dtype.names:
            if key == "np_ijr":
                c.p[key] = mpi[key].reshape((c.p["ixr"],c.p["jxr"]),order="F")
            else:
                c.p[key] = mpi[key].reshape((mpi[key].size),order="F")
    
        c.p["i2ir"] = c.p["i2ir"][marginx:ixg-marginx]
        c.p["j2jr"] = c.p["j2jr"][marginy:jxg-marginy]

        c.p["iss"] = c.p["iss"] - 1
        c.p["iee"] = c.p["iee"] - 1
        c.p["jss"] = c.p["jss"] - 1
        c.p["jee"] = c.p["jee"] - 1

######################################################
######################################################
######################################################
### read horizontal variable
### prepare array qq = np.zeros((mtype+3,jx,kx))
def read_qq_select(dir,xs,n,silent=False):
    import numpy as np
    import config as c
    i0 = np.argmin(np.abs(c.p["x"]-xs))
    ir0 = c.p["i2ir"][i0]
    mtype = c.p["mtype"]
    iixl = c.p["iixl"]
    jjxl = c.p["jjxl"]
    jx = c.p["kx"]
    kx = c.p["kx"]
    iss = c.p["iss"]
    jss = c.p["jss"]
    jee = c.p["jee"]

    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in c.q2:
        memflag = not c.q2['ro'].shape == (jx,kx)
    if 'ro' not in c.q2 or memflag:
        print('memory is newly allocated')
        c.q2["ro"] = np.zeros((jx,kx))
        c.q2["vx"] = np.zeros((jx,kx))
        c.q2["vy"] = np.zeros((jx,kx))
        c.q2["vz"] = np.zeros((jx,kx))
        c.q2["bx"] = np.zeros((jx,kx))
        c.q2["by"] = np.zeros((jx,kx))
        c.q2["bz"] = np.zeros((jx,kx))
        c.q2["se"] = np.zeros((jx,kx))
        c.q2["pr"] = np.zeros((jx,kx))
        c.q2["te"] = np.zeros((jx,kx))
        c.q2["op"] = np.zeros((jx,kx))
    for jr0 in range(1,c.p["jxr"]+1):
        np0 = c.p["np_ijr"][ir0-1,jr0-1]
        dtyp=np.dtype([ \
                ("qq",c.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("pr",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("te",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("op",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ])
        f = open(dir+"remap/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
        qqq = np.fromfile(f,dtype=dtyp,count=1)
        c.q2["ro"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[0,i0-iss[np0],:,:]
        c.q2["vx"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[1,i0-iss[np0],:,:]
        c.q2["vy"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[2,i0-iss[np0],:,:]
        c.q2["vz"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[3,i0-iss[np0],:,:]
        c.q2["bx"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[4,i0-iss[np0],:,:]
        c.q2["by"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[5,i0-iss[np0],:,:]
        c.q2["bz"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[6,i0-iss[np0],:,:]
        c.q2["se"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[7,i0-iss[np0],:,:]
        c.q2["pr"][jss[np0]:jee[np0]+1,:] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        c.q2["te"][jss[np0]:jee[np0]+1,:] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        c.q2["op"][jss[np0]:jee[np0]+1,:] = qqq["op"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        f.close()

    if not silent:
        print('### variales are stored in c.q2 ###')
    
    return

######################################################
######################################################
######################################################
### read horizontal variable
### prepare array qq = np.zeros((mtype+3,ix,jx,kx))
def read_qq(dir,n,silent=False):
    import numpy as np
    import config as c
    
    mtype = c.p["mtype"]
    iixl = c.p["iixl"]
    jjxl = c.p["jjxl"]
    kx = c.p["kx"]
    iss = c.p["iss"]
    iee = c.p["iee"]
    jss = c.p["jss"]
    jee = c.p["jee"]
    ix = c.p["ix"]
    jx = c.p["jx"]
    kx = c.p["kx"]

    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in c.q3:
        memflag = not c.q3['ro'].shape == (ix,jx,kx)
    if 'ro' not in c.q3 or memflag:
        print('memory is newly allocated')
        c.q3["ro"] = np.zeros((ix,jx,kx))
        c.q3["vx"] = np.zeros((ix,jx,kx))
        c.q3["vy"] = np.zeros((ix,jx,kx))
        c.q3["vz"] = np.zeros((ix,jx,kx))
        c.q3["bx"] = np.zeros((ix,jx,kx))
        c.q3["by"] = np.zeros((ix,jx,kx))
        c.q3["bz"] = np.zeros((ix,jx,kx))
        c.q3["se"] = np.zeros((ix,jx,kx))
        c.q3["pr"] = np.zeros((ix,jx,kx))
        c.q3["te"] = np.zeros((ix,jx,kx))
        c.q3["op"] = np.zeros((ix,jx,kx))

    for ir0 in range(1,c.p["ixr"]+1):
        for jr0 in range(1,c.p["jxr"]+1):
            np0 = c.p["np_ijr"][ir0-1,jr0-1]
            dtyp=np.dtype([ \
                    ("qq",c.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ("pr",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ("te",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ("op",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ])
            f = open(dir+"remap/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
            qqq = np.fromfile(f,dtype=dtyp,count=1)
            c.q3["ro"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[0,:,:,:]
            c.q3["vx"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[1,:,:,:]
            c.q3["vy"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[2,:,:,:]
            c.q3["vz"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[3,:,:,:]
            c.q3["bx"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[4,:,:,:]
            c.q3["by"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[5,:,:,:]
            c.q3["bz"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[6,:,:,:]
            c.q3["se"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[7,:,:,:]
            c.q3["pr"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            c.q3["te"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            c.q3["op"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["op"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            f.close()

    if not silent:
        print('### variales are stored in c.q3 ###')
    
    return

######################################################
######################################################
### read all 2D variable
### prepare array qq = np.zeros((mtype+5,ix,jx))
def read_qq_2d(dir,n):
    import numpy as np
    import config as c

    ix = c.p["ix"]
    jx = c.p["jx"]

    qq1 = np.zeros((c.p["mtype"]+5,ix,jx))

    f = open(dir+"remap/qq.dac."+'{0:08d}'.format(n),"rb")
    qq0 = np.fromfile(f,c.p["endian"]+'f',(c.p["mtype"]+5)*c.p["ix"]*c.p["jx"])
    f.close()
    qq1 = np.reshape(qq0,(c.p["mtype"]+5,c.p["ix"],c.p["jx"]),order="F")

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
def read_tau_one(dir,n,silent=False):
    import numpy as np
    import config as c
    f = open(dir+"remap/qq_tu.dac."+'{0:08d}'.format(n),"rb")
    qq_in0 = np.fromfile(f,c.p["endian"]+'f',c.p["m_tu"]*c.p["m_in"]*c.p["jx"]*c.p["kx"])
    f.close()

    c.qi["in"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,0,:,:]
    c.qi["ro"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,1,:,:]
    c.qi["se"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,2,:,:]
    c.qi["pr"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,3,:,:]
    c.qi["te"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,4,:,:]
    c.qi["vx"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,5,:,:]
    c.qi["vy"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,6,:,:]
    c.qi["vz"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,7,:,:]
    c.qi["bx"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,8,:,:]
    c.qi["by"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,9,:,:]
    c.qi["bz"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,10,:,:]
    c.qi["he"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,11,:,:]
    c.qi["fr"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[0,12,:,:]

    c.qi["in01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,0,:,:]
    c.qi["ro01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,1,:,:]
    c.qi["se01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,2,:,:]
    c.qi["pr01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,3,:,:]
    c.qi["te01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,4,:,:]
    c.qi["vx01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,5,:,:]
    c.qi["vy01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,6,:,:]
    c.qi["vz01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,7,:,:]
    c.qi["bx01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,8,:,:]
    c.qi["by01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,9,:,:]
    c.qi["bz01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,10,:,:]
    c.qi["he01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,11,:,:]
    c.qi["fr01"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[1,12,:,:]

    c.qi["in001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,0,:,:]
    c.qi["ro001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,1,:,:]
    c.qi["se001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,2,:,:]
    c.qi["pr001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,3,:,:]
    c.qi["te001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,4,:,:]
    c.qi["vx001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,5,:,:]
    c.qi["vy001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,6,:,:]
    c.qi["vz001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,7,:,:]
    c.qi["bx001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,8,:,:]
    c.qi["by001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,9,:,:]
    c.qi["bz001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,10,:,:]
    c.qi["he001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,11,:,:]
    c.qi["fr001"] = np.reshape(qq_in0,(c.p["m_tu"],c.p["m_in"],c.p["jx"],c.p["kx"]),order="F")[2,12,:,:]

    if not silent:
        print('### variales are stored in c.qi ###')

    return 

def read_time(dir,n):
    import numpy as np
    import config as c
    f = open(dir+"time/t.dac."+'{0:08d}'.format(n),"rb")
    t = np.fromfile(f,c.p['endian']+'d',1)
    f.close()    
    t = np.reshape(t,(1),order="F")[0]

    return t
    
##############################
# read remap_calc variable
def read_vc(dir,n,silent=False):
    import numpy as np
    import config as c

    f = open(dir+"remap/vla.dac."+'{0:08d}'.format(n),"rb")
    vl0 = np.fromfile(f,c.p["endian"]+'f',c.p['m2da']*c.p['ix']*c.p['jx'])
    f.close()

    vl = np.reshape(vl0,(c.p['ix'],c.p['jx'],c.p['m2da']),order="F")

    for m in range(c.p["m2da"]):
        c.vc[c.p["cl"][m]] = vl[:,:,m]

    if not silent:
        print('### variales are stored in c.vc ###')
    
    return
