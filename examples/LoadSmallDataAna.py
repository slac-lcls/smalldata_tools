import sys
import numpy as np
import argparse
from matplotlib import pyplot as plt
import socket
import os
from smalldata_tools import droplets,hist2d,getUserData,rebin,addToHdf5
from smalldata_tools import dropObject
from IPython.terminal.prompts import Prompts,Token
import itertools

class myPrompt(Prompts):
    def in_prompt_tokens(self, cli=None):
        return [(Token.Prompt, 'SDAna In ['), (Token.PromptNum, str(self.shell.execution_count)), (Token.Prompt, ']: ' ),]
    def out_prompt_tokens(self):
        return [(Token.Prompt, 'SDAna Out ['), (Token.PromptNum, str(self.shell.execution_count)), (Token.Prompt, ']: ' ),]
        #return [(Token.Prompt, 'LDAna Out:' ),]

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

if not args.exp:
    hostname=socket.gethostname()
    if hostname.find('xpp')>=0:
        hutch='XPP'
    elif hostname.find('xcs')>=0:
        hutch='XCS'
    elif hostname.find('mfx')>=0:
        hutch='MFX'
    else:
        #then check current path
        path=os.getcwd()
        if path.find('xcs')>=0:
            hutch='XCS'
        elif path.find('xpp')>=0:
            hutch='XPP'
        elif path.find('mfx')>=0:
            hutch='MFX'
        else:
            sys.exit()
    try:
        import RegDB.experiment_info
        expname=RegDB.experiment_info.active_experiment(hutch)[1]
    except:
        print 'could not determine experiment name, will quit'
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

from smalldata_tools import SmallDataAna
ana = None
anaps = None
try:
    from smalldata_tools import SmallDataAna_psana
    anaps = SmallDataAna_psana(expname,run,dirname,fname)
except:
    pass

if anaps is not None and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
    print 'create ana module from anaps'
    ana = anaps.sda
else:
    print 'we will now try to open the littleData file directly'
    ana = SmallDataAna(expname,run, dirname, fname)
    if 'fh5' not in ana.__dict__.keys():
        ana = None

if ana is not None:
    ana.printRunInfo()

    #somewhat non-obvious: if lower==upper; REJECT this point (e.g. FLTPOS==0)
    ana.addCut('lightStatus/xray',0.5,1.5,'xon')
    ana.addCut('lightStatus/xray',0.5,1.5,'on')
    ana.addCut('lightStatus/laser',0.5,1.5,'on')
    ana.addCut('lightStatus/xray',0.5,1.5,'off')
    ana.addCut('lightStatus/laser',-0.5,0.5,'off')
    
    ana.addCut('ipm3/sum',0.03,10.,'on')
    ana.addCut('ipm3/sum',0.03,10.,'off')
    ana.addCut('tt/AMPL',0.025,10.,'on')

    binRange = np.arange(13.,15.5,0.03)
    if run==65:
        binRange = np.arange(8.5,13,0.03)
    if run==153 or run == 304:
        binRange = np.arange(12,15,0.05)
    ana.addCube('cube','delay',binRange,'on')
    ana.addToCube('cube',['ipm2/sum','ipm3/sum','diodeU/channels','cs140_rob/ROI_0_sum'])
    ana.addCube('cubeOff','delay',binRange,'off')
    ana.addToCube('cubeOff',['ipm2/sum','ipm3/sum','diodeU/channels','cs140_rob/ROI_0_sum'])
    #ana.getVar('cs140_rob/ROI')
    #ana.addToCube('cube',['ipm2/sum','ipm3/sum','diodeU/channels','cs140_rob/ROI_sum','cs140_rob/ROI'])

    cubeData = ana.makeCubeData('cube')
    cubeDataOff = ana.makeCubeData('cubeOff')

    #works for run 65.
    if run == 65:
        plt.plot(cubeData['binVar_bins'],cubeData['cs140_rob__ROI_0_sum']/cubeData['ipm2__sum'],'o-')
        plt.plot(cubeDataOff['binVar_bins'],cubeDataOff['cs140_rob__ROI_0_sum']/cubeDataOff['ipm2__sum'],'o-')

    #works for run 153.
    if run == 153 or run == 304:
        plt.plot(cubeData['binVar_bins'],cubeData['diodeU__channels'][:,1]/cubeData['ipm2__sum'],'o-')
        plt.plot(cubeDataOff['binVar_bins'],cubeDataOff['diodeU__channels'][:,1]/cubeDataOff['ipm2__sum'],'o-')


    ##stuff for tim's data.
    ##this here is specifically for xpp00316, Tim's data.
    #ana.addCut('ipm3/sum',0.5,10.,'on')
    #ana.addCut('ipm3/sum',0.5,10.,'off')
    #ana.addCut('tt/AMPL',0.03,10.,'on')

    ##delay=ana.getDelay()
    ##encOffset=1306521.2
    #encOffset=1343965.
    #conv=20e-9
    #speedLight = 2.99792458e8
    #encScale=2 * conv / speedLight
    ##ana.xrData = ana.xrData.assign(newDelay = (ana.xrData.enc__ch0-encOffset)*encScale*1e12 + ana.xrData.tt_ttCorr))
    #newDelay = (ana.xrData.enc__ch0-encOffset)*encScale*1e12 + ana.xrData.tt__ttCorr
    #ana.addVar('newDelay', newDelay)

    #ana.addCube('cube','newDelay',np.arange(-2.,2.,0.025),'on')
    #ana.xrData = ana.xrData.assign(cspad__ROI1 = (ana.xrData.cspad__ROI_0_sum + ana.xrData.cspad__ROI_1_sum)) 
    #ana.xrData = ana.xrData.assign(cspad__ROI2 = (ana.xrData.cspad__ROI_2_sum + ana.xrData.cspad__ROI_3_sum  + ana.xrData.cspad__ROI_4_sum)) 
    ##ana.addToCube('cube',['ipm2/sum','ipm3/sum','diodeU/channels','/cspad/ROI_0_sum','/cspad/ROI_1_sum','/cspad/ROI_2_sum','/cspad/ROI_3_sum','/cspad/ROI_4_sum'])
    #azav = ana.getVar('cspad/azav')
    #ana.addToCube('cube',['ipm2/sum','ipm3/sum','diodeU/channels','/cspad/ROI_0_sum','/cspad/ROI_1_sum','/cspad/ROI_2_sum','/cspad/ROI_3_sum','/cspad/ROI_4_sum','cspad/azav'])

    #cubeData = ana.makeCubeData('cube')

    #azav = cubeData['cspad__azav'].mean(axis=1)

    #azavt0_2d = cubeData['cspad__azav'][:10].mean(axis=0)
    #azavt0_2d = azavt0_2d / azavt0_2d.max()

    #azavt10t20_2d = cubeData['cspad__azav'][10:20].mean(axis=0)
    #azavt10t20_2d = azavt10t20_2d / azavt10t20_2d.max()
    #azavt10t20_2d = azavt10t20_2d / azavt0_2d
    #azavt50t60_2d = cubeData['cspad__azav'][50:60].mean(axis=0)
    #azavt50t60_2d = azavt50t60_2d / azavt50t60_2d.max()
    #azavt50t60_2d = azavt50t60_2d / azavt0_2d
    #azavtm10_2d = cubeData['cspad__azav'][-10:].mean(axis=0)
    #azavtm10_2d = azavtm10_2d / azavtm10_2d.max()
    #azavtm10_2d = azavtm10_2d / azavt0_2d

    #f, ax = plt.subplots(nrows=2, ncols=2)
    #ax[0,0].imshow(azavt10t20_2d, interpolation='none', aspect='auto')
    #ax[1,0].imshow(azavt50t60_2d, interpolation='none', aspect='auto')
    #ax[1,1].imshow(azavtm10_2d, interpolation='none', aspect='auto')

    #azavnorm = np.array([ azavb/azavb.max() for azavb in azav])
    #azavt0 = azav[:10].mean(axis=0)
    #azavt0 = azavt0/azavt0.max()
    #azavnrel = np.array([azavn/azavt0 for azavn in azavnorm])
    #plt.imshow(azavnrel, interpolation='none', aspect='auto')
    #plt.figure()
    #plt.imshow(azavnrel,  interpolation='none', aspect='auto')

