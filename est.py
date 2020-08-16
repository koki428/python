import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import R2D2
import sys

try:
    caseid
except NameError:
    print("input caseid id (3 digit)")
    caseid = 0
    caseid = input()
    caseid = "d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd"]:
    n0 = d.p["nd"]

print("Maximum time step= ",nd," time ="\
      ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

#n0 = 27
#nd = 27

#n0 = 1
#nd = 1

if geometry == 'Spherical':
    xx,yy = np.meshgrid(x,y,indexing='ij')
    sinyy = sin(yy)
    sinyym = np.average(sinyy,axis=1)

    xx,yy_flux = np.meshgrid(x_flux,y,indexing='ij')
    sinyy_flux = sin(yy_flux)
    sinyym_flux = np.average(sinyy_flux,axis=1)
else:
    sinyy = 1
    sinyym = 1.
    
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

fet = np.zeros((ix+1,nd-n0+1))
fdt = np.zeros((ix+1,nd-n0+1))
fkt = np.zeros((ix+1,nd-n0+1))
frt = np.zeros((ix+1,nd-n0+1))
ftt = np.zeros((ix+1,nd-n0+1))

plt.close('all')
plt.clf()

for n in range(n0,nd+1):
    print(n)
    ##############################
    # read time
    t = d.read_time(n)
    
    ##############################
    # read value
    d.read_vc(n,silent=True)

    print(d.vc['sem'].mean())
    #print(d.vc['rzm'].mean())
    print(d.vc['rom'].mean())

    ##############################    
    fsun = 6.306e10
    if geometry == 'Spherical':
        fsun = 3.86e33/pi/4
        fe = np.average(d.vc["fe"]*sinyy_flux,axis=1)/sinyym_flux*x_flux**2
        fd = np.average(d.vc["fd"]*sinyy_flux,axis=1)/sinyym_flux*x_flux**2
        fk = np.average(d.vc["fk"]*sinyy_flux,axis=1)/sinyym_flux*x_flux**2
        fr = np.average(d.vc["fr"]*sinyy_flux,axis=1)/sinyym_flux
    else:
        fe = np.average(d.vc["fe"],axis=1)/sinyym
        fd = np.average(d.vc["fd"],axis=1)/sinyym
        fk = np.average(d.vc["fk"],axis=1)/sinyym
        fr = np.average(d.vc["fr"],axis=1)/sinyym
    #fc = np.average(d.vc["fa"],axis=1)
    
    xs = rsun - 2.e8
    ds = 2.e7
    sr = 0.5e0*(1.e0 + np.tanh((x_flux-xs)/ds))
    SR, sry = np.meshgrid(sr,y,indexing="ij")
    
    ff = fd*sr + fe*(1.e0-sr)
    ft = ff + fk + fr 

    vxrmst[:,n-n0] = np.sqrt(np.average(d.vc["vxrms"]**2*sinyy,axis=1))/sinyym
    vyrmst[:,n-n0] = np.sqrt(np.average(d.vc["vyrms"]**2*sinyy,axis=1))/sinyym
    vzrmst[:,n-n0] = np.sqrt(np.average(d.vc["vzrms"]**2*sinyy,axis=1))/sinyym

    bxrmst[:,n-n0] = np.sqrt(np.average(d.vc["bxrms"]**2*sinyy,axis=1))/sinyym
    byrmst[:,n-n0] = np.sqrt(np.average(d.vc["byrms"]**2*sinyy,axis=1))/sinyym
    bzrmst[:,n-n0] = np.sqrt(np.average(d.vc["bzrms"]**2*sinyy,axis=1))/sinyym

    rormst[:,n-n0] = np.sqrt(np.average(d.vc["rorms"]**2*sinyy,axis=1))/sinyym
    sermst[:,n-n0] = np.sqrt(np.average(d.vc["serms"]**2*sinyy,axis=1))/sinyym
    prrmst[:,n-n0] = np.sqrt(np.average(d.vc["prrms"]**2*sinyy,axis=1))/sinyym
    termst[:,n-n0] = np.sqrt(np.average(d.vc["terms"]**2*sinyy,axis=1))/sinyym

    romt[:,n-n0] = np.average(d.vc["rom"]*sinyy,axis=1)/sinyym
    semt[:,n-n0] = np.average(d.vc["sem"]*sinyy,axis=1)/sinyym
    prmt[:,n-n0] = np.average(d.vc["prm"]*sinyy,axis=1)/sinyym
    temt[:,n-n0] = np.average(d.vc["tem"]*sinyy,axis=1)/sinyym

    fet[:,n-n0] = fe
    fdt[:,n-n0] = fd
    fkt[:,n-n0] = fk
    frt[:,n-n0] = fr
    ftt[:,n-n0] = ft

    fontsize = 12
    fmax = 2.0
    fmin = -1.0

    plt.rcParams["font.size"] = 15
    fig1 = plt.figure(num=200,figsize=(12,8))
    ax1 = fig1.add_subplot(221)
    ax2 = fig1.add_subplot(222)
    ax3 = fig1.add_subplot(223)
    ax4 = fig1.add_subplot(224)


    #####################
    if geometry == 'Spherical':
        xp = x_flux/rsun
        xlabel = r'$r/R_\odot$'
        xpp = x/rsun
    else:
        xp = (x_flux - rsun)*1.e-8
        xpp = (x - rsun)*1.e-8
        xlabel = r'$x-R_\odot\ \mathrm{[Mm]}$'

        
    ax1.plot(xp,ff/fsun,label=r'$F_\mathrm{e}$',color="red")
    ax1.plot(xp,fk/fsun,label=r'$F_\mathrm{k}$',color="green")
    ax1.plot(xp,fr/fsun,label=r'$F_\mathrm{r}$',color="blue")
    ax1.plot(xp,ft/fsun,label=r'$F_\mathrm{t}$',color="black")

    ax1.hlines(y=1,xmin=xp.min(),xmax=xp.max(),linestyle='--',color='black')
    ax1.set_ylim(fmin,fmax)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("$F/F_{\odot}$")
    ax1.set_title("Energy fluxes")
    ax1.legend()

    #####################
    vxrms = np.sqrt((d.vc['vxrms']**2).mean(axis=1))
    vhrms = np.sqrt((d.vc['vyrms']**2 + d.vc['vzrms']**2).mean(axis=1))
    ax2.plot(xpp,vxrms*1.e-5,label=r'$v_{x\mathrm{(rms)}}$',color='blue')
    ax2.plot(xpp,vhrms*1.e-5,label=r'$v_\mathrm{h(rms)}$',color='red')
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(r"velocities [km/s]")
    ax2.set_label('RMS velocities')
    if deep_flag == 0:
        ax2.set_yscale('log')
    ax2.legend()

    #####################
    bxrms = np.sqrt((d.vc['bxrms']**2).mean(axis=1))
    bhrms = np.sqrt((d.vc['byrms']**2 + d.vc['bzrms']**2).mean(axis=1))
    ax3.plot(xpp,bxrms,label=r'$B_{x\mathrm{(rms)}}$',color='blue')
    ax3.plot(xpp,bhrms,label=r'$B_\mathrm{h(rms)}$',color='red')
    ax3.set_xlabel(xlabel)
    ax3.set_ylabel(r"Magnetic field [G]")
    ax3.set_label('RMS magnetic field')
    #ax3.set_yscale('log')
    ax3.legend()

    hp = pr0/ro0/gx
    semx = np.zeros(ix)
    for i in range(1,ix):
        semx[i] = (semt[i,n-n0] - semt[i-1,n-n0])/(x[i] - x[i-1])

    ax4.plot(semx*hp/cp)
    
    
    if n == n0:
        plt.tight_layout()

    ax3.annotate(s="t="+"{:.2f}".format(t/3600./24.)+" [day]"\
                     ,xy=[0.01,0.01],xycoords="figure fraction",fontsize=18)

    plt.pause(0.001)

    
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

np.savez(d.p['datadir']+"est.npz"\
             ,x=x,y=y,z=z,rsun=rsun\
             ,ro0=ro0,pr0=pr0,te0=te0,se0=se0\
             ,vxrms=vxrms,vyrms=vyrms,vzrms=vzrms\
             ,bxrms=bxrms,byrms=byrms,bzrms=bzrms\
             ,rorms=rorms,prrms=prrms,serms=serms,terms=terms\
             ,rom=rom,prm=prm,sem=sem,tem=tem\
             ,ff=ff,fk=fk,fr=fr,ft=ft\
             )
         
plt.rcParams["font.size"] = 15
fig2 = plt.figure(num=100,figsize=(12,5))
ax23 = fig2.add_subplot(121)
ax24 = fig2.add_subplot(122)
ax23.plot(xp,ff/fsun,color="red",label="$F_\mathrm{e}$")
ax23.plot(xp,fk/fsun,color="green",label="$F_\mathrm{k}$")
ax23.plot(xp,fr/fsun,color="blue",label="$F_\mathrm{r}$")
ax23.plot(xp,ft/fsun,color="black",label="$F_\mathrm{t}$")
#ax23.set_xlim(xmin/rsun,xmax/rsun)
ax23.set_ylim(fmin,fmax)
ax23.set_xlabel(xlabel)
ax23.set_ylabel("$F/F_{\odot}$")
ax23.set_title("Full convection zone")
ax23.legend(loc='upper left',prop={'size': 15})
ax23.annotate(s="t="+"{:.2f}".format(t/3600./24.)+" [day]"\
                 ,xy=[0.01,0.01],xycoords="figure fraction",fontsize=18)

ax23.hlines(y=1,xmin=xmin/rsun,xmax=xmax/rsun,linestyle='--',color='black')

ax24.plot(x_flux,ff/fsun,color="red")
ax24.plot(x_flux,fk/fsun,color="green")
ax24.plot(x_flux,fr/fsun,color="blue")
ax24.plot(x_flux,ft/fsun,color="black")
ax24.set_ylim(fmin,fmax)
ax24.set_xlabel("$x - R_{\odot} \ [\mathrm{Mm}]$")
ax24.set_ylabel("$F/F_{\odot}$")
ax24.set_title("Around photosphere")
fig2.tight_layout()
plt.pause(0.001)
plt.ion()
    
