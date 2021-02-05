import R2D2
try:
    caseid
except NameError:
    print("input caseid id (3 digit)")
    caseid = 0
    caseid = input()
    caseid = "d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"
casedir="../figs/"+caseid
os.makedirs(casedir,exist_ok=True)

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd"]:
    n0 = d.p["nd"]

n = 200
ixc = 3072
jxc = 3072
kxc = 3072
var = 'se'

READ    = False
CONVERT = True
WRITE   = True

file = 'test.vtk'

print('### Time step =',n)
if READ:
    print('- Read data')
    d.read_qq_variable(200,var)
    d.read_vc(200)

    qqm0  ,tmp = np.meshgrid(d.vc[var+'m']  ,z,indexing='ij')
    qqrms0,tmp = np.meshgrid(d.vc[var+'rms'],z,indexing='ij')

    qqm   = qqm0.reshape((ix,jx,kx))
    qqrms = qqrms0.reshape((ix,jx,kx))

    qqs = (d.qv - qqm)/qqrms

if CONVERT:
    print('- Convert Spherical to Cartesian ###')
    qqc,xc,yc,zc = R2D2.geometry_convert.spherical2cartesian(qqs,x/rsun,y,z,ixc,jxc,kxc)

if WRITE:
    print('- Write data in VTK format')
    R2D2.vtk.write_3D(qqc,xc,yc,zc,file,var)