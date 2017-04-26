# importing generic python modules
import psana
import abc
import numpy as np
#import more specific psana functions here
import xtcav.ShotToShotCharacterization as Xtcav  

#
# classes for default detector types
#
class defaultDetector(object):
    __metaclass__ = abc.ABCMeta
    def __init__(self, detname, name):
        self.name=name
        self.detname=detname
        self.det=psana.Detector(detname)
    def inRun(self):
        dNames=[]
        for dn in psana.DetNames():
            for dnn in dn:
                if dnn!='':
                    dNames.append(dnn)
        if self.detname in dNames:
            return True
        return False
    @abc.abstractmethod
    def data(self,evt):
        """method that should return a dict of values from event"""
        

class lightStatus(defaultDetector):
    def __init__(self, detname='evr0', codes=[[162],[]]):
        defaultDetector.__init__(self, detname, 'lightStatus')
        self.xrayCodes = codes[0]
        self.laserCodes = codes[1]

    def data(self,evt):
        xfel_status, laser_status = (1,1) # default if no EVR code matches
        dl={}
        evtCodes = self.det.eventCodes(evt)
        if evtCodes is not None:
            for xOff in self.xrayCodes:
                if xOff in evtCodes:
                    xfel_status = 0
            for lOff in self.laserCodes:
                if lOff in evtCodes:
                    laser_status = 0
        dl['xray']=xfel_status
        dl['laser']=laser_status
        return dl
        
class ipmDetector(defaultDetector):
    def __init__(self, detname, name=None, savePos=False):
        if name is None:
            self.name = detname
        else:
            self.name = name
        defaultDetector.__init__(self, detname, name)
        self.savePos = savePos
    def data(self, evt):
        dl={}
        if self.det.sum(evt) is not None:
            dl['sum']=self.det.sum(evt)
            dl['channels']=self.det.channel(evt)
            if self.savePos:
                dl['xpos']=self.det.xpos(evt)
                dl['ypos']=self.det.ypos(evt)
        return dl

class bmmonDetector(defaultDetector):
    def __init__(self, detname, name=None, savePos=False):
        if name is None:
            self.name = detname
        else:
            self.name = name
        defaultDetector.__init__(self, detname, name)
        self.savePos = savePos
    def data(self, evt):
        dl={}
        data = self.det.get(evt)
        dl['sum']=data.TotalIntensity()
        dl['channels']=data.peakA()
        if self.savePos:
            dl['xpos']=data.X_Position()
            dl['ypos']=data.Y_Position()
        return dl

class epicsDetector(defaultDetector):
    def __init__(self, name='epics', PVlist=[]):
        self.name = name
        self.detname='epics'
        self.PVlist = PVlist
        self.pvs=[]
        for pv in PVlist:
            try:
                self.pvs.append(psana.Detector(pv))
            except:
                print 'could not find EPICS PV %s in data'%pv
    def inRun(self):
        if len(self.pvs)>0:
            return True
        return False

    def data(self,evt):
        dl={}
        for pvname,pv in zip(self.PVlist,self.pvs):
            dl[pvname]=pv()
        return dl

class encoderDetector(defaultDetector):
    def __init__(self, detname, name=None):
        if name is None:
            self.name = detname
        else:
            self.name = name
        defaultDetector.__init__(self, detname, name)
    def data(self, evt):
        dl={}
        for desc,value in zip(self.det.descriptions(), self.det.values(evt)):
            if desc!='':
                dl[desc]=value
        return dl

class controlDetector(defaultDetector):
    def __init__(self, name='scan'):
        defaultDetector.__init__(self, 'ControlData', 'scan')
        try:
            self.stepPV = psana.Detector('scan_current_step')
        except:
            try:
                self.stepPV = psana.Detector('Scan_current_step')
            except:
                self.stepPV = None
    def data(self, evt):
        dl={}
        if self.stepPV is not None:
            dl['varStep']=self.stepPV()
        for icpv,cpv in enumerate(self.det().pvControls()):
            dl['var%d'%icpv]=cpv.value()
            dl[cpv.name()]=cpv.value()
        return dl

class aiDetector(defaultDetector):
    def __init__(self, detname, name=None):
        if name is None:
            self.name = detname
        else:
            self.name = name
        defaultDetector.__init__(self, detname, name)
        self.aioInfo = [[ i for i in range(0,16)], [ 'ch%02d'%i for i in range(0,16)], [ 1. for i in range(0,16)], [ 0. for i in range(0,16)]]

    def setPars(self, AIOPars):
        if len(AIOPars)<2:
            print 'need 2/3 lists: channel#, user-friendly names & conversion factors (optional)'
            return
        self.aioInfo[0] = AIOPars[0]
        self.aioInfo[1] = AIOPars[1]
        if len(AIOPars)==3:
            self.aioInfo[2] = AIOPars[2]
            if len(AIOPars)==4:
                self.aioInfo[3] = AIOPars[3]
            else:
                self.aioInfo[3] = [0. for entry in AIOPars[0]]
        else:
            self.aioInfo[2] = [1. for entry in AIOPars[0]]

    def data(self, evt):
        dl={}
        for ichn,chName,chnScale,chnOffset in zip(self.aioInfo[0], self.aioInfo[1], self.aioInfo[2], self.aioInfo[3]):
            dl[chName]=self.det.get(evt).channelVoltages()[ichn]*chnScale+chnOffset
        return dl

class ttDetector(defaultDetector):
    def __init__(self, name='tt', baseName='TTSPEC:'):
        self.name = name
        self.detname='epics'
        self.ttNames = ['FLTPOS','FLTPOS_PS','AMPL','FLTPOSFWHM','REFAMPL','AMPLNXT']
        self.PVlist = [ baseName+pvname for pvname in self.ttNames ]
        self.pvs=[]
        for pv in self.PVlist:
            try:
                self.pvs.append(psana.Detector(pv))
            except:
                print 'could not find timetool EPICS PV %s in data'%pv
        self.ttCalib=None
    def inRun(self):
        if len(self.pvs)>0:
            return True
        return False
    def setPars(self, calibPars):
      if calibPars != None:
        self.ttCalib = calibPars

    def data(self,evt):
        dl={}
        for ttname,pvname,pv in zip(self.ttNames,self.PVlist,self.pvs):
            dl[ttname]=pv()
        ttOrg = dl[self.ttNames[1]]
        if self.ttCalib is None:
            dl['ttCorr']=ttOrg
        else:
            dl['ttCorr']=self.ttCalib[0] + self.ttCalib[1]*ttOrg
            if len(self.ttCalib)>2:
                dl['ttCorr']+=ttOrg*ttOrg*self.ttCalib[2]
        return dl

class damageDetector(defaultDetector):
    def __init__(self, name='damage'):
        self.name = name
        self.detNames=[]
        for dn in psana.DetNames():
            if dn[1]!='':
                self.detNames.append(dn[1])
            else:
                self.detNames.append(dn[0])
        self.detAlias=[ det for det in self.detNames]
    def inRun(self):
        return True
    def setPars(self, detList):
        for det in detList:
            try:
                if det.detname in self.detNames:
                    self.detAlias[self.detNames.index(det.detname)]=det.name
            except:
                pass    
    def data(self,evt):
        #check if detectors are in event
        dl={}
        aliases = [ k.alias() for k in evt.keys() ]
        srcNames = [ k.src().__str__().replace(')','').replace('BldInfo(','').replace('DetInfo(','') for k in evt.keys() ]
        for det,alias in zip(self.detNames,self.detAlias):
            val=0
            if det in aliases:
                val=1
            elif det in srcNames:
                val=1
            dl[det.replace('-','_')]=val
            if alias!=det:
                dl[alias.replace('-','_')]=val
        return dl

#
# needs testing with data.
#
class ttRawDetector(defaultDetector):
    def __init__(self, name='ttRaw', env=None):
        self.name = name
        self.detname = ''
        self.kind='stepDown'
        self.weights=None
        self.ttROI_signal=None
        self.ttROI_sideband=None
        self.ttROI_reference=None
        self.ttProj=False
        self.runningRef=None
        self.refitData=False
        self.useProjection=False
        self.beamOff=[]
        self.laserOff=[]
        self.sb_convergence=1.
        self.ref_convergence=1.
        self.subtract_sideband=False
        self.ttCalib = [0.,1.]
        ttCfg=None
        self.evrdet=psana.Detector('NoDetector.0:Evr.0')
        if env is None:
            env = psana.Env
            #getting the env this way unfortunatly does not work. Find out how psana.DetNames() does it.
            return
        for cfgKey in env.configStore().keys():
            if cfgKey.type() == psana.TimeTool.ConfigV2:
                ttCfg = env.configStore().get(psana.TimeTool.ConfigV2, cfgKey.src())
                self.detname = cfgKey.alias()
                defaultDetector.__init__(self, self.detname, 'ttRaw')
        if ttCfg is not None:
            self.ttProj=ttCfg.write_projections()
            self.ttROI_signal = [[ttCfg.sig_roi_lo().row(),ttCfg.sig_roi_hi().row()],\
                                 [ttCfg.sig_roi_lo().column(),ttCfg.sig_roi_hi().column()]]
            self.ttROI_sideband = [[ttCfg.sb_roi_lo().row(),ttCfg.sb_roi_hi().row()],\
                                   [ttCfg.sb_roi_lo().column(),ttCfg.sb_roi_hi().column()]]
            if ttCfg.use_reference_roi()>0:
                self.ttROI_reference = [[ttCfg.ref_roi_lo().row(),ttCfg.ref_roi_hi().row()],\
                                        [ttCfg.ref_roi_lo().column(),ttCfg.ref_roi_hi().column()]]
            else:
                self.ttROI_reference = self.ttROI_signal

            self.weights = ttCfg.weights()
            self.sb_convergence=ttCfg.sb_convergence()
            self.ref_convergence=ttCfg.ref_convergence()
            self.subtract_sideband=ttCfg.subtract_sideband()
            self.ttCalib = ttCfg.calib_poly()
                
            for el in ttCfg.beam_logic():
                self.beamOff.append(el.event_code())
            for el in ttCfg.laser_logic():
                self.laserOff.append(el.event_code())
        else:
            defaultDetector.__init__(self, self.detname, 'ttRaw')
            
    def inRun(self):
        if self.detname=='':
            return False

    def setPars(self, ttPars):
        if 'ttProj' in ttPars.keys():
            self.ttProj=ttPars['ttProj']
        if 'ttROI_signal' in ttPars.keys():
            self.ttROI_signal=ttPars['ttROI_signal']
        if 'ttROI_reference' in ttPars.keys():
            self.ttROI_reference=ttPars['ttROI_reference']
        if 'ttROI_sideband' in ttPars.keys():
            self.ttROI_sideband=ttPars['ttROI_sideband']
        if 'weights' in ttPars.keys():
            self.weights=ttPars['weights']
        if 'runningRef' in ttPars.keys():
            self.runningRef=ttPars['runningRef']
        if 'refitData' in ttPars.keys():
            self.refitData=ttPars['refitData']
        if 'useProjection' in ttPars.keys():
            self.useProjection=ttPars['useProjection']
        if 'sb_convergence' in ttPars.keys():
            self.sb_convergence=ttPars['sb_convergence']
        if 'ref_convergence' in ttPars.keys():
            self.ref_convergence=ttPars['ref_convergence']
        if 'subtract_sideband' in ttPars.keys():
            self.subtract_sideband=ttPars['subtract_sideband']
        if 'ttCalib' in ttPars.keys():
            self.ttCalib=ttPars['ttCalib']

    def data(self, evt):
        retDict = self.getTraces(evt)
        if not self.refitData:
            return retDict
        data = self.prepareTrace(evt, retDict)
        fitDict = self.fitTraceData(data)
        for key in fitDict:
            retDict[key] = fitDict[key]
        return retDict

    def getTraces(self, evt):
        #check that we have reference, otherwise replace by none
        evtCodes = self.evrdet.eventCodes(evt)
        ttDet=self.det
        ttData={}
        ttData['tt_signal']=np.zeros(abs(self.ttROI_signal[1][1]-self.ttROI_signal[1][0]))
        if self.ttROI_sideband is not None:
            ttData['tt_sideband']=np.zeros(abs(self.ttROI_sideband[1][1]-self.ttROI_sideband[1][0]))
        #XXX
        if self.ttROI_reference is not None:
            ttData['tt_reference']=np.zeros(abs(self.ttROI_reference[1][1]-self.ttROI_reference[1][0]))
        for lOff in self.laserOff:
            if lOff in evtCodes:
                return ttData
            
        try:
            ttData['tt_signal_pj']=ttDet.projected_signal().astype(dtype='uint32').astype(float)
            ttData['tt_sideband_pj']=ttDet.projected_sideband().astype(dtype='uint32').astype(float)
            ttData['tt_reference_pj']=ttDet.projected_reference().astype(dtype='uint32').astype(float)
        except:
            pass
        ttImg = ttDet.raw(evt)
        if self.ttROI_signal is not None:
            ttData['tt_signal']=ttImg[self.ttROI_signal[0][0]:self.ttROI_signal[0][1],self.ttROI_signal[1][0]:self.ttROI_signal[1][1]].mean(axis=0)          
        if self.ttROI_sideband is not None:
            ttData['tt_sideband']=ttImg[self.ttROI_sideband[0][0]:self.ttROI_sideband[0][1],self.ttROI_sideband[1][0]:self.ttROI_sideband[1][1]].mean(axis=0)
        if self.ttROI_reference is not None:
            beamOff=False
            for bOff in self.beamOff:
                if bOff in evtCodes:
                    beamOff=True
            if beamOff:
                ttRef = ttImg[self.ttROI_reference[0][0]:self.ttROI_reference[0][1],self.ttROI_reference[1][0]:self.ttROI_reference[1][1]].mean(axis=0)
                if self.runningRef is None:
                    self.runningRef=ttRef
                else:
                    self.runningRef=ttRef*self.ref_convergence + self.runningRef*(1.-self.ref_convergence)
                #print 'update self.runningRef'
            ttData['tt_reference']=self.runningRef        
        return ttData

    def prepareTrace(self, evt, ttData=None):
        if ttData is None:
            ttData = self.getTraces(evt)
            if len(ttData.keys())==0:
                return None

        if self.useProjection:
            ttRef=ttData['tt_reference_pj']
            ttSignal=ttData['tt_signal_pj']
            if self.subtract_sideband>0:
                ttSignal-=ttData['tt_sideband_pj']
                ttRef-=ttData['tt_sideband_pj']
        else:
            ttRef=ttData['tt_reference'].copy()
            ttSignal=ttData['tt_signal']
            if self.subtract_sideband>0:
                ttSignal-=ttData['tt_sideband']
                ttRef-=ttData['tt_sideband']

        if ttData['tt_reference'].sum()==0:
            nanArray = np.ones(ttData['tt_reference'].shape[0])
            nanArray=nanArray*np.nan
            return nanArray

        return ttSignal/ttRef

    def fitTraceData(self, data):
        if data is None or len(data)<10:
            return
        lf = len(self.weights)
        halfrange = round(lf/10)
        retDict = {}
        retDict['pos']=0.
        retDict['amp']=0.
        retDict['fwhm']=0.
        retDict['pos_ps']=0.
        if np.isnan(data).sum()==data.shape[0]:
            return retDict
        
        f0 = np.convolve(np.array(self.weights).ravel(),data,'same')
        f = f0[lf/2:len(f0)-lf/2-1]
        retDict['f']=f
        if (self.kind=="stepDown"):
            mpr = f.argmin()
        else:
            mpr = f.argmax()
        # now do a parabolic fit around the max
        xd = np.arange(max(0,mpr-halfrange),min(mpr+halfrange,len(f)-1))
        yd = f[max(0,mpr-halfrange):min(mpr+halfrange,len(f)-1)]
        p2 = np.polyfit(xd,yd,2)
        tpos = -p2[1]/2./p2[0]
        tamp = np.polyval(p2,tpos)
        try:
            if self.kind == 'stepDown':
                beloh = (f>((f[-25:].mean()+tamp)/2.)).nonzero()[0]-mpr
            else:
                beloh = (f<tamp/2).nonzero()[0]-mpr            
            #print 'beloh ',len(beloh[beloh<0]),len(beloh[beloh>0])
            tfwhm = abs(beloh[beloh<0][-1]-beloh[beloh>0][0])
        except:
            tfwhm = 0.
        if self.kind == 'stepDown':
            tamp = abs(f[-25:].mean()-tamp)
        retDict['pos']=tpos + lf/2.
        retDict['amp']=tamp
        if np.isnan(tamp): 
            retDict['fwhm']=np.nan
        else:
            retDict['fwhm']=tfwhm 
        ttOrg = retDict['pos']
        ttCorr = self.ttCalib[0]+ ttOrg*self.ttCalib[1]
        if len(self.ttCalib)>2:
            ttCorr+=ttOrg*ttOrg*self.ttCalib[2]
        retDict['pos_ps']=ttCorr
        #for k in retDict.keys():
        #    print 'ret ',k,retDict[k]

        return retDict

class xtcavDetector(defaultDetector):
    def __init__(self, name='xtcav', detname='xtcav'):
        self.name = name
        self.detname = detname
        self.nb=1
        self.size=5000
        self.ShotToShotCharacterization = Xtcav.ShotToShotCharacterization()
        #this will NOT work. Have to ask Chris how to get access like psana does.
        self.ShotToShotCharacterization.SetEnv(psana.Env)
    def setPars(self, xtcavPars):
        self.nb=xtcavPars[0]
        self.size=xtcavPars[1]
    def data(self,evt):
        #check if detectors are in event
        dl={}
        xtcav_success=False
        self.ShotToShotCharacterization.SetCurrentEvent(evt)
        if self.ShotToShotCharacterization.SetCurrentEvent(evt):
            timeOrg,power,ok=self.ShotToShotCharacterization.XRayPower()
            if ok:
              arSize=timeOrg[0].shape[0]
              agreement,ok=self.ShotToShotCharacterization.ReconstructionAgreement()
              xtcav_success=True
              if arSize>=1 and arSize<=self.size:
                  time = np.append(timeOrg[0],np.array([np.nan] * (self.size-arSize)))
                  power = np.append(power[0],np.array([np.nan] * (self.size-arSize)))
              else:
                print 'Xtcav array is too small in run, please check littleData configuration',env.run()
                time = timeOrg[:self.size]
                power = power[:self.size]
              dl['agreement']=agreement
              dl['arSize']=arSize
              dl['time']=time
              dl['power']=power
                    
        return dl


