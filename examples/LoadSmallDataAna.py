import sys
import numpy as np
import argparse
from matplotlib import pyplot as plt
import socket
import os
import time
#temporary 'fix'?
#try:
#    import xarray as xr
#except:
#    import xarray as xr
from smalldata_tools.SmallDataUtils import getUserData
from smalldata_tools.utilities import hist2d,rebin,addToHdf5
from IPython.terminal.prompts import Prompts,Token
import itertools

class myPrompt(Prompts):
    def in_prompt_tokens(self, cli=None):
        return [(Token.Prompt, 'SDAna In ['), (Token.PromptNum, str(self.shell.execution_count)), (Token.Prompt, ']: ' ),]
    def out_prompt_tokens(self):
        return [(Token.Prompt, 'SDAna Out ['), (Token.PromptNum, str(self.shell.execution_count)), (Token.Prompt, ']: ' ),]

ip=get_ipython()
ip.prompts = myPrompt(ip)

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run",type=int)
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--dir", help="output directory")
parser.add_argument("--file", help="file name")
args = parser.parse_args()

if not args.run:
    run=int(raw_input('provide a run number of experiment %s:'%expname))
else:
    run=args.run

hutches=['amo','sxr','xpp','xcs','mfx','cxi','mec']
if not args.exp:
    hostname=socket.gethostname()
    foundHutch=False
    for ihutch in hutches:
        if hostname.find(ihutch)>=0:
            hutch=ihutch.upper()
            foundHutch=True
            break
    if not foundHutch and hostname.find('psusr')>=0:
        if hostname.find('psusr11')>=0:
            hutch='AMO'
        elif hostname.find('psusr12')>=0:
            hutch='SXR'
        elif hostname.find('psusr13')>=0:
            hutch='XPP'
        elif hostname.find('psusr21')>=0:
            hutch='XCS'
        elif hostname.find('psusr22')>=0:
            hutch='CXI'
        elif hostname.find('psusr23')>=0:
            hutch='MEC'
        elif hostname.find('psusr24')>=0:
            hutch='MFX'
        if hutch!='':
            foundHutch=True
    else:
        #then check current path
        path=os.getcwd()
        for ihutch in hutches:
            if path.find(ihutch)>=0:
                hutch=ihutch.upper()
                foundHutch=True
                break
        if not foundHutch:
            print('I cannot figure out which hutch we are in, so cannot determine the current experiment')
            sys.exit()

    try:
        import logging
        import requests
        ws_url = "https://pswww.slac.stanford.edu/ws/lgbk"
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        if hutch == 'cxi':
            print('Will assume the first CXI station, if this is wrong, please  -e <expname> on commandline')
        resp = requests.get(ws_url + "/lgbk/ws/activeexperiment_for_instrument_station", {"instrument_name": hutch, "station": 0})
        expname = resp.json().get("value", {}).get("name")
    except:
        print('could not determine experiment name, will quit')
        sys.exit()
else:
    expname=args.exp

dirname=''
if args.dir:
    dirname = args.dir
    if dirname.find('/')<=0:
        dirname+='/'
fname=''
if args.file:
    fname = args.file

#from smalldata_tools import SmallDataAna
sys.path.append('./smalldata_tools')
from smalldata_tools.SmallDataAna import SmallDataAna
ana = None
anaps = None
t0 = time.time()
try:
    from smalldata_tools.SmallDataAna_psana import SmallDataAna_psana
    #from smalldata_tools import SmallDataAna_psana
    anaps = SmallDataAna_psana(expname,run,dirname,fname)
except:
    print('failed to import & create SmallDataAna_psana object')
    pass

if anaps is not None and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
    print('use ana module from anaps')
    ana = anaps.sda
else:
    print('we will now try to open the littleData file directly')
    ana = SmallDataAna(expname,run, dirname, fname)
    if 'fh5' not in ana.__dict__.keys():
        ana = None

if ana is not None:
    ana.printRunInfo()

    #anaps.AvImage('cspad',numEvts=20)
    #from skimage.feature import canny
    #from utilities_FitCenter import applyCanny
    #ar = anaps.cspad.det.image(320, anaps.AvImg_pedSub_cspad)
    #mask = anaps.cspad.det.image(320, anaps.AvImg_pedSub_cspad,anaps.cspad.mask)
    #arThres, arSparse = applyCanny(ar, mask.astype(bool), sigma=1, thres=95)
    
    ana.addCut('lightStatus/xray',-0.5,0.5,'xoff')
    ana.addCut('lightStatus/xray',0.5,1.5,'xon')

    #ana.addCut('gas_detector/f_22_ENRC',-0.5,0.05,'gdetLow')
    #anaps.makePedestal('epix',filterName='gdetLow', dirname='calib',numEvts=1000)
    #ana.addCube('ipm3Cube','ipm3/sum',[0.,2.,4],useFilter='xon')
    #ana.addCube('ipm3Cube','ipm3/sum',[0.,2.,4],useFilter='xon',addBinVars=['gas_detector/f_22_ENRC',[2.6,4.,3]])
    #ana.addToCube('ipm3Cube',['ipm3/sum', 'diodeU/channels'])
    #epixDict = {'source':'epix_2','full':1,'image':1}
    #ana.addToCube('ipm3Cube',['ipm3/sum', 'diodeU/channels', epixDict])
    #ana.addToCube('ipm3Cube',['ipm2/sum'], isIdxVar=True)
    #ana.addToCube('ipm3Cube',['droplet:epix_2/droplet:image'])
    #ana.addToCube('ipm3Cube',['droplet:epix_2/droplet:array'])

    #cd = ana.makeCubeData('ipm3Cube')
    #cd,rd = ana.makeCubeData('ipm3Cube', toHdf5='h5')
    #anaps.makeCubeData('ipm3Cube', nEvtsPerBin=1)

    #ana.addCut('ipm2/sum',1.2,5.,'good')
    #ana.addCube('delayCube','delay',[0.,0.8,80],useFilter='good')
    #ana.addToCube('delayCube',['ipm2/sum', 'diodeU/channels', epixDict])
    #ana.addToCube('delayCube',['droplet:epix_vonHamos/droplet:image', epixDict])
    #ana.addToCube('delayCube',['droplet:epix_vonHamos/droplet:array', epixDict])
    #cd,rd = ana.makeCubeData('delayCube', toHdf5='h5', onoff=0)
    #cd,rd = ana.makeCubeData('delayCube', toHdf5='h5', onoff=1)

    #from CubeAna import CubeAna
    #cube = CubeAna(ana.expname, ana.run, cubeName='snelsonTest')
    #cube.addAzInt(detname='cspad', center=[70000.,11000.],dis_to_sam=110., eBeam=9.5)
    #cube.addDataFromRun(3)
    #cube.applyAzav()
    #cube.cubeSumToHdf5()

    ##cubeOff = CubeAna(ana.expname, ana.run, cubeName='ccmE_laserOff')
    ##offData = cubeOff.plotCube(sig=['diodeU__channels',3], plotWith='returnData')
    #cube = CubeAna(ana.expname, ana.run, cubeName='ipm3Cube')
    #cube.addDataFromRun(178)
    #cube.extendData()

    #ana.addCut('lightStatus/xray',0.5,1.5,'good')
    #ana.addCut('ipm5/sum',10, 1000.,'good')
    ##varList = ['ipm5/sum','snd_dio/t4d','snd_dio/dco','scan/_delay']
    #epixDict = {'source':'epix_ladm_1','full':1, 'image':1}
    #varList = ['ipm5/sum','snd_dio/t4d','snd_dio/dco','scan/_delay', epixDict]

    ##artificial cube to speed data reading up.
    #myrand=np.random.rand(ana.getVar('fiducials').shape[0])
    #ana.addVar('random/random',myrand)

    #ana.addCube('lowBinCube_1d4','snd_dio/t4d',np.linspace(1000.,5000.,4),'good')
    #ana.addCube('lowBinCube_1d9','snd_dio/t4d',np.linspace(1000.,5000.,9),'good') #1-d: need bound-1 bins, for 2-d NOT true! Fix that too!
    #ana.addCube('lowBinCube_2d','snd_dio/t4d',np.linspace(1000.,5000.,4),'good', addBinVars={'user/random':[0.,1.,2]})
    #ana.addCube('lowBinCube_2dB','snd_dio/t4d',np.linspace(1000.,5000.,4),'good', addBinVars={'random/random':[0.,1.,2]})
    #ana.addToCube('lowBinCube_1d4',varList) #list variables to be binned.
    #ana.addToCube('lowBinCube_1d9',varList)
    #ana.addToCube('lowBinCube_2d',varList) 
    #ana.addToCube('lowBinCube_2dB',varList) 

    #ana.addCube('DelayLineScan_new','scan/_delay',[],'good')
    ##ana.addCube('DelayLineScan_new','scan/_delay',np.linspace(0,0.1,11),'good')
    #ana.addToCube('DelayLineScan_new',varList)

    #cubeName='ipm5_laser'
    #binVarName='lightStatus/laser' #always there. Need to make two cubes.
    #binBounds=[-0.5,1.5,2]
    
    #ipm5 = ana.getVar('ipm5/sum')
    #binBounds_ipm5=np.percentile(ipm5,np.linspace(0.,1.,11)*100)
    #print 'ipm5: ',binBounds_ipm5
    ##we can't really do this. needs to be the same across runs.
    #binBounds_ipm5=[  -50 ,    80.  ,   130. ,   180. ,   230.,   300.,\
    #                  380. ,   490.  ,   650.  ,   1000. ,  1e6]

    #ana.addCube(cubeName,'ipm5/sum',binBounds_ipm5,'xon', addBinVars={'lightStatus/laser':[-0.5,1.5,2]}) #set event selection
    #ana.addToCube(cubeName,varList) #list variables to be binned.

    #anaps.AvImage('icarus_pink',numEvts=10)                        
    #anaps.AvImage('jungfrau512k',numEvts=10)

    #xoff = ana.getVar('lightStatus/xray')
    #pre_xoff = np.append([0], xoff[0:-1])
    #ana.addVar('pre_xoff',pre_xoff)
    #ana.addCut('pre_xoff',-0.5,0.5,'pre_xoff')

    #offsel = ana.getFilter('xoff')
    #print 'we have %d off events'%(offsel.sum())

    #fids = ana.getVar('fiducials')
    #epixDict = {'source':'epix10ka2m','common_mode':180,'full':1}
    #ana.addCut('ipm5/sum',200,1e6,'good')
    #ana.addCube('ringCube','fiducials',[fids.min(), fids.max(), fids.shape[0]/120],useFilter='good')
    #ana.addToCube('ringCube',['ipm5/sum', epixDict])
    #cd = ana.makeCubeData('ringCube')#, toHdf5='h5', onoff=0)
    
    #detector monitoring
    #ana.addCut('gas_detector/f_22_ENRC',-0.5,0.15,'xoff_gdet')
    #anaps.makePedestal('epix',i0Check='ipm',dirname='/reg/d/psdm/xpp/xpplp7515/results/detector_monitoring_new/', useFilter='xoff_gdet', numEvts=1000)

    ana.addCut('lightStatus/xray',0.5,1.5,'on')
    ana.addCut('lightStatus/laser',0.5,1.5,'on')

    
t1 = time.time()
print('this took %g seconds '%(t1-t0))
