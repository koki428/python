######################################################
######################################################
######################################################
### for initialize google spread sheet

def init_gspread(json_key,project):
    '''
    This function initialize the utility of google spread 

    Parameters:
        json_key (str): file of json key to access Google API
        projet (str): project name, typically name of upper directory

    Returnes:
        None
    '''
    import gspread
    from oauth2client.service_account import ServiceAccountCredentials

    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']

    credentials = ServiceAccountCredentials.from_json_keyfile_name(json_key, scope)
    gc = gspread.authorize(credentials)
    wks = gc.open(project).sheet1

    wks.update_acell('A1', 'Case ID')
    wks.update_acell('B1', 'Server')
    wks.update_acell('C1', '(ix,jx,kx)')
    wks.update_acell('D1', 'xmin [Mm]')
    wks.update_acell('E1', 'xmax [Mm]')
    wks.update_acell('F1', 'ymin [Mm]')
    wks.update_acell('G1', 'ymax [Mm]')
    wks.update_acell('H1', 'zmin [Mm]')
    wks.update_acell('I1', 'zmax [Mm]')
    wks.update_acell('J1', 'uniform')
    wks.update_acell('K1', 'dx [km]')
    wks.update_acell('L1', 'm ray')
    wks.update_acell('M1', 'dtout [s]')
    wks.update_acell('N1', 'dtout_tau [s]')
    wks.update_acell('O1', 'alpha')
    wks.update_acell('P1', 'RSST')
    wks.update_acell('Q1', 'Omega[Sun]')
    wks.update_acell('R1', 'upodate time')
    wks.update_acell('S1', 'origin')

######################################################
######################################################
def out_gspread(self,caseid,json_key,project):
    '''
    This function output parameters to 
    Google spread sheet

    Parameters:
        caseid (str): caseid
        json_key (str): file of json key to access Google API
        projet (str): project name, typically name of upper directory
        d (R2D2_data): instance of R2D2_data class
    
    '''
    import datetime
    import gspread
    
    from oauth2client.service_account import ServiceAccountCredentials

    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']

    credentials = ServiceAccountCredentials.from_json_keyfile_name(json_key, scope)
    gc = gspread.authorize(credentials)
    wks = gc.open(project).sheet1

    str_id = str(int(caseid[1:])+1)
    
    wks.update_acell('A'+str_id, caseid)
    wks.update_acell('B'+str_id, self.p['server'])
    wks.update_acell('C'+str_id, str(self.p['ix'])+' '+str(self.p['jx'])+' '+str(self.p['kx']))
    wks.update_acell('D'+str_id, '{:6.2f}'.format((self.p['xmin']-self.p['rsun'])*1.e-8))
    wks.update_acell('E'+str_id, '{:6.2f}'.format((self.p['xmax']-self.p['rsun'])*1.e-8))
    wks.update_acell('F'+str_id, '{:6.2f}'.format(self.p['ymin']*1.e-8))
    wks.update_acell('G'+str_id, '{:6.2f}'.format(self.p['ymax']*1.e-8))
    wks.update_acell('H'+str_id, '{:6.2f}'.format(self.p['zmin']*1.e-8))
    wks.update_acell('I'+str_id, '{:6.2f}'.format(self.p['zmax']*1.e-8))
    if ((self.p['x'][1] - self.p['x'][0]) == (self.p['x'][self.p['ix']-1] - self.p['x'][self.p['ix']-2])):
        wks.update_acell('J'+str_id,'T')
    else:
        wks.update_acell('J'+str_id,'F')
    dx0 = (self.p['x'][1] - self.p['x'][0])*1.e-5
    dx1 = (self.p['x'][self.p['ix']-1] - self.p['x'][self.p['ix']-2])*1.e-5
    wks.update_acell('K'+str_id, '{:6.2f}'.format(dx0)+' '+'{:6.2f}'.format(dx1))
    wks.update_acell('L'+str_id, self.p['rte'])
    wks.update_acell('M'+str_id, '{:6.2f}'.format(self.p['dtout']))
    wks.update_acell('N'+str_id, '{:6.2f}'.format(self.p['dtout_tau']))
    wks.update_acell('O'+str_id, '{:5.2f}'.format(self.p['potential_alpha']))
    if self.p['xi'].max() == 1.0:
        wks.update_acell('P'+str_id,'F')
    else:
        wks.update_acell('P'+str_id,'T')
    
    wks.update_acell('Q'+str_id, '{:5.1f}'.format(self.p['omfac']))
    wks.update_acell('R'+str_id,str(datetime.datetime.now()).split('.')[0])
    wks.update_acell('S'+str_id,self.p['origin'])
