import R2D2
import glob
import os

try:
    caseid
except NameError:
    print("input caseid id (3 digit)")
    caseid = 0
    caseid = input()
    caseid = "d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"

R2D2.init(datadir)
for key in R2D2.p:
    exec('%s = %s%s%s' % (key, 'R2D2.p["',key,'"]'))

json_key = glob.glob(os.environ['HOME']+'/json/*')[0]
project = os.getcwd().split('/')[-2]

R2D2.init_gspread(json_key,project)
R2D2.out_gspread(caseid,json_key,project)
