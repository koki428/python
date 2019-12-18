import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import R2D2
import config as c
import sys

try:
    caseid
except NameError:
    print("input caseid id (3 digit)")
    caseid = 0
    caseid = input()
    caseid = "d"+caseid.zfill(3)

dir="../run/"+caseid+"/data/"

R2D2.read_init(dir,"3d")
for key in c.p:
    exec('%s = %s%s%s' % (key, 'c.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > c.p["nd"]:
    n0 = c.p["nd"]

print("Maximum time step= ",nd," time ="\
      ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

vxrmst = np.zeros((ix,nd-n0+1))
vyrmst = np.zeros((ix,nd-n0+1))
vzrmst = np.zeros((ix,nd-n0+1))

bxrmst = np.zeros((ix,nd-n0+1))
byrmst = np.zeros((ix,nd-n0+1))
bzrmst = np.zeros((ix,nd-n0+1))

rormst = np.zeros((ix,nd-n0+1))
sermst = np.zeros((ix,nd-n0+1))
prrmst = np.zeros((ix,nd-n0+1))
termst = np.zeros((ix,nd-n0+1))

romt = np.zeros((ix,nd-n0+1))
semt = np.zeros((ix,nd-n0+1))
prmt = np.zeros((ix,nd-n0+1))
temt = np.zeros((ix,nd-n0+1))

fet = np.zeros((ix,nd-n0+1))
fdt = np.zeros((ix,nd-n0+1))
fkt = np.zeros((ix,nd-n0+1))
frt = np.zeros((ix,nd-n0+1))
ftt = np.zeros((ix,nd-n0+1))

for n in range(n0,nd+1):
#for n in range(0,1):
    print(n)
    ##############################
    # read time
    f = open(dir+"time/t.dac."+'{0:08d}'.format(n),"rb")
    t = np.fromfile(f,endian+'d',1)
    f.close()    
    t = np.reshape(t,(1),order="F")
    
    ##############################
    # read value
    R2D2.read_vc(dir,n,silent=True)
        
    ##############################    
    fsun = 6.318e10
    fe = np.average(c.vc["fe"],axis=1)
    fd = np.average(c.vc["fd"],axis=1)
    fk = np.average(c.vc["fk"],axis=1)
    fr = np.average(c.vc["fr"],axis=1)

    xs = c.p["rsun"] - 2.e8
    ds = 2.e7
    sr = 0.5e0*(1.e0 + np.tanh((c.p["x"]-xs)/ds))
    SR, sry = np.meshgrid(sr,y,indexing="ij")
    
    ff = fd*sr + fe*(1.e0-sr)
    ft = ff + fk + fr

    vxrmst[:,n-n0] = np.sqrt(np.average(c.vc["vxrms"]**2,axis=1))
    vyrmst[:,n-n0] = np.sqrt(np.average(c.vc["vyrms"]**2,axis=1))
    vzrmst[:,n-n0] = np.sqrt(np.average(c.vc["vzrms"]**2,axis=1))

    bxrmst[:,n-n0] = np.sqrt(np.average(c.vc["bxrms"]**2,axis=1))
    byrmst[:,n-n0] = np.sqrt(np.average(c.vc["byrms"]**2,axis=1))
    bzrmst[:,n-n0] = np.sqrt(np.average(c.vc["bzrms"]**2,axis=1))

    rormst[:,n-n0] = np.sqrt(np.average(c.vc["rorms"]**2,axis=1))
    sermst[:,n-n0] = np.sqrt(np.average(c.vc["serms"]**2,axis=1))
    prrmst[:,n-n0] = np.sqrt(np.average(c.vc["prrms"]**2,axis=1))
    termst[:,n-n0] = np.sqrt(np.average(c.vc["terms"]**2,axis=1))

    romt[:,n-n0] = np.sqrt(np.average(c.vc["rom"]**2,axis=1))
    semt[:,n-n0] = np.sqrt(np.average(c.vc["sem"]**2,axis=1))
    prmt[:,n-n0] = np.sqrt(np.average(c.vc["prm"]**2,axis=1))
    temt[:,n-n0] = np.sqrt(np.average(c.vc["tem"]**2,axis=1))

    fet[:,n-n0] = fe
    fdt[:,n-n0] = fd
    fkt[:,n-n0] = fk
    frt[:,n-n0] = fr
    ftt[:,n-n0] = ft

    fontsize = 12
    fmax = 3.0
    fmin = -2.0

    plt.rcParams["font.size"] = 15
    fig1 = plt.figure(num=1,figsize=(12,5))
    ax1 = fig1.add_subplot(121)
    ax2 = fig1.add_subplot(122)
    ax1.plot(c.p["xr"],ff/fsun,color="red",label="$F_\mathrm{e}$")
    ax1.plot(c.p["xr"],fk/fsun,color="green",label="$F_\mathrm{k}$")
    ax1.plot(c.p["xr"],fr/fsun,color="blue",label="$F_\mathrm{r}$")
    ax1.plot(c.p["xr"],ft/fsun,color="black",label="$F_\mathrm{t}$")
    ax1.set_xlim(c.p["xmin"]/c.p["rsun"],c.p["xmax"]/c.p["rsun"])
    ax1.set_ylim(fmin,fmax)
    ax1.set_xlabel("$x/R_{\odot}$")
    ax1.set_ylabel("$F/F_{\odot}$")
    ax1.set_title("Full convection zone")
    ax1.legend(loc='upper left',prop={'size': 10})
    ax1.annotate(s="t="+"{:.2f}".format(t[0]/3600./24.)+" [day]"\
                     ,xy=[0.01,0.01],xycoords="figure fraction",fontsize=18)
    
    xtmp = 0.5*(x[ix//4*3] + x[ix//4*3-1])
    if deep_top_flag == 1:
        ax1.vlines(xtmp/rsun,-2,3)

    #####################
    ax2.plot(c.p["xn"],ff/fsun,color="red")
    ax2.plot(c.p["xn"],fk/fsun,color="green")
    ax2.plot(c.p["xn"],fr/fsun,color="blue")
    ax2.plot(c.p["xn"],ft/fsun,color="black")
    ax2.set_xlim(-20,1)
    ax2.set_ylim(fmin,fmax)
    ax2.set_xlabel("$x - R_{\odot} \ [\mathrm{Mm}]$")
    ax2.set_ylabel("$F/F_{\odot}$")
    ax2.set_title("Around photosphere")

    if n == n0:
        plt.tight_layout()
    plt.pause(0.001)
    #plt.draw()
    
    if n != nd:
        plt.clf() # clear figure

    # loop end
    ###############################################################################
    ###############################################################################
    ###############################################################################

vxrms = np.sqrt(np.average(vxrmst**2,axis=1))
vyrms = np.sqrt(np.average(vyrmst**2,axis=1))
vzrms = np.sqrt(np.average(vzrmst**2,axis=1))

bxrms = np.sqrt(np.average(bxrmst**2,axis=1))
byrms = np.sqrt(np.average(byrmst**2,axis=1))
bzrms = np.sqrt(np.average(bzrmst**2,axis=1))

rorms = np.sqrt(np.average(rormst**2,axis=1))
serms = np.sqrt(np.average(sermst**2,axis=1))
prrms = np.sqrt(np.average(prrmst**2,axis=1))
terms = np.sqrt(np.average(termst**2,axis=1))

rom = np.average(romt,axis=1)
sem = np.average(semt,axis=1)
prm = np.average(prmt,axis=1)
tem = np.average(temt,axis=1)

fe = np.average(fet,axis=1)
fd = np.average(fdt,axis=1)

ff = fd*sr + fe*(1.e0-sr)

fk = np.average(fkt,axis=1)
fr = np.average(frt,axis=1)
ft = ff + fk + fr

np.savez(dir+"est.npz"\
             ,x=x,y=y,z=z,rsun=rsun\
             ,ro0=ro0,pr0=pr0,te0=te0,se0=se0\
             ,vxrms=vxrms,vyrms=vyrms,vzrms=vzrms\
             ,bxrms=bxrms,byrms=byrms,bzrms=bzrms\
             ,rorms=rorms,prrms=prrms,serms=serms,terms=terms\
             ,rom=rom,prm=prm,sem=sem,tem=tem\
             ,ff=ff,fk=fk,fr=fr,ft=ft\
             )
         
plt.rcParams["font.size"] = 15
fig2 = plt.figure(num=2,figsize=(12,5))
ax3 = fig2.add_subplot(121)
ax4 = fig2.add_subplot(122)
ax3.plot(xr,ff/fsun,color="red",label="$F_\mathrm{e}$")
ax3.plot(xr,fk/fsun,color="green",label="$F_\mathrm{k}$")
ax3.plot(xr,fr/fsun,color="blue",label="$F_\mathrm{r}$")
ax3.plot(xr,ft/fsun,color="black",label="$F_\mathrm{t}$")
ax3.set_xlim(xmin/rsun,xmax/rsun)
ax3.set_ylim(fmin,fmax)
ax3.set_xlabel("$x/R_{\odot}$")
ax3.set_ylabel("$F/F_{\odot}$")
ax3.set_title("Full convection zone")
ax3.legend(loc='upper left',prop={'size': 15})
ax3.annotate(s="t="+"{:.2f}".format(t[0]/3600./24.)+" [day]"\
                 ,xy=[0.01,0.01],xycoords="figure fraction",fontsize=18)

ax4.plot(xn,ff/fsun,color="red")
ax4.plot(xn,fk/fsun,color="green")
ax4.plot(xn,fr/fsun,color="blue")
ax4.plot(xn,ft/fsun,color="black")
ax4.set_xlim(-20,1)
ax4.set_ylim(fmin,fmax)
ax4.set_xlabel("$x - R_{\odot} \ [\mathrm{Mm}]$")
ax4.set_ylabel("$F/F_{\odot}$")
ax4.set_title("Around photosphere")
fig2.tight_layout()
plt.pause(0.001)
plt.ion()
    
