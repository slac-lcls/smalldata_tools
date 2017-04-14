import sys
import numpy as np
import argparse
from matplotlib import pyplot as plt
import socket
import os
from smalldata_tools import droplets,hist2d
from IPython.terminal.prompts import Prompts,Token

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

if anaps and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
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

    #drop = droplets(ana.fh5,'epix_alcove','droplet')
    ana.addCube('cube','ipm2/sum',np.arange(0.,2.,0.5),'on')
    ana.addToCube('cube',['ipm2/sum','ipm3/sum','diodeU/channels'])
    ana.getVar('/cs140_rob/ROI')
    ana.addToCube('cube',['cs140_rob/ROI'])
    #delay=ana.getDelay()
    #vals = ana.fh5.get_node('/cs140_rob','ROI').read().squeeze()
    #ana.addVar(name='newVar',data=vals)
