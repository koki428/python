import R2D2
import glob
import os

print("input caseid id (3 digit) or just type enter")
caseid_tmp = 0
caseid_tmp = input()
print(caseid_tmp)
if caseid_tmp != '':
    caseid = "d"+caseid_tmp.zfill(3)

datadir="../run/"+caseid+"/data/"

ReadFlag = False
try:
    d
except NameError:
    ReadFlag = True
else:
    if d.p['datadir'] != datadir:
        ReadFlag = True

if ReadFlag:
    d = R2D2.R2D2_data(datadir)
    for key in d.p:
        exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

json_key = glob.glob(os.environ['HOME']+'/json/*')[0]
project = os.getcwd().split('/')[-2]

R2D2.google.init_gspread(json_key,project)
R2D2.google.out_gspread(d,caseid,json_key,project)
