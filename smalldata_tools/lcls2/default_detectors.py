# importing generic python modules
import psana
import abc
import numpy as np
try:
    basestring
except NameError:
    basestring = str
from smalldata_tools.common.BaseDetector import DefaultDetector_base


defaultDetector = DefaultDetector_base


class genericDetector(defaultDetector):
    def __init__(self,  name=None, run=None, h5name=None):
        if name is None:
            self.name = 'anydet'
        else:
            self.name = name
        if h5name is None: h5name = self.name
        defaultDetector.__init__(self, self.name, h5name, run)

    def data(self, evt):
        dl = {}
        raw = getattr( self.det, 'raw')
        vetolist = ['TypeId', 'Version', 'config']
        if raw is not None:
           fields = [ field for field in dir(raw) if (field[0]!='_' and field not in vetolist) ]
           for field in fields:
               if getattr(raw, field)(evt) is None: continue
               dl[field] = getattr(raw, field)(evt)
               if isinstance(dl[field], list): dl[field]=np.array(dl[field])
        return dl


class ttDetector(defaultDetector):
    def __init__(self,  name=None, run=None, saveTraces=False):
        if name is None:
            self.name = 'anydet'
        else:
            self.name = name
        self.saveTraces = saveTraces
        defaultDetector.__init__(self, self.name, 'tt', run)

    def data(self, evt):
        dl={}
        fex=getattr( self.det, 'ttfex')
        veto_fields = ['TypeId', 'Version', 'calib', 'image', 'raw', 'config' ]
        if fex is not None:
           fields = [ field for field in dir(fex) if (field[0]!='_' and field not in veto_fields) ]
           for field in fields:
               if getattr(fex, field)(evt) is None: continue
               dl[field]=getattr(fex, field)(evt)
               if isinstance(dl[field], list): dl[field]=np.array(dl[field])

        if self.saveTraces:
            fex=getattr( self.det, 'ttproj')
            veto_fields = ['TypeId', 'Version', 'calib', 'image', 'raw', 'config' ]
            if fex is not None:
               fields = [ field for field in dir(fex) if (field[0]!='_' and field not in veto_fields) ]
               for field in fields:
                   if getattr(fex, field)(evt) is None: continue
                   dl[field]=getattr(fex, field)(evt)
                   if isinstance(dl[field], list): dl[field]=np.array(dl[field])        
        return dl


class fimfexDetector(defaultDetector):
    def __init__(self,  name=None, run=None):
        if name is None:
            self.name = 'anydet'
        else:
            self.name = name
        defaultDetector.__init__(self, self.name, self.name, run)

    def data(self, evt):
        dl={}
        fex=getattr( self.det, 'fex')
        veto_fields = ['TypeId', 'Version', 'calib', 'image', 'raw', 'config' ]
        if fex is not None:
           fields = [ field for field in dir(fex) if (field[0]!='_' and field not in veto_fields) ]
           for field in fields:
               if getattr(fex, field)(evt) is None: continue
               dl[field]=getattr(fex, field)(evt)
               if isinstance(dl[field], list): dl[field]=np.array(dl[field])
        return dl


class lightStatus(defaultDetector):
    def __init__(self, codes, run):
        defaultDetector.__init__(self, 'timing', 'lightStatus', run)
        self.xrayCodes_drop = [ c for c in codes[0] if c > 0]
        self.laserCodes_drop = [ c for c in codes[1] if c > 0]
        self.xrayCodes_req = [ -c for c in codes[0] if c < 0]
        self.laserCodes_req =  [ -c for c in codes[1] if c < 0]

    def data(self,evt):
        xfel_status, laser_status = (1,1) # default if no EVR code matches
        dl={}
        evtCodes = getattr(getattr( self.det, 'raw'), 'eventcodes')(evt)
        if evtCodes is not None:
            for xOff in self.xrayCodes_drop:
                if evtCodes[xOff]:
                    xfel_status = 0
            for lOff in self.laserCodes_drop:
                if evtCodes[lOff]:
                    laser_status = 0
            if len(self.xrayCodes_req)>0 and xfel_status==1:
                xfel_status = 0
                for xOff in self.xrayCodes_req:
                    if evtCodes[xOff]:
                        xfel_status = 1
            if len(self.laserCodes_req)>0 and laser_status==1:
                laser_status = 0
                for lOff in self.laserCodes_req:
                    if evtCodes[lOff]:
                        laser_status = 1

        else:
            xfel_status, laser_status = (-1,-1) # default if no EVR code matches
        dl['xray']=xfel_status
        dl['laser']=laser_status
        return dl


class epicsDetector(defaultDetector):
    def __init__(self, name='epics', PVlist=[],run=None):
        self.name = name
        self.detname='epics'
        self.PVlist = []
        self.missing = []
        self.missingPV = []
        self.pvs = []
        # run.epicsinfo is an alias or (alias, pvname) --> pvname dict.
        # Let's rearrange this somewhat to make it more useful.
        self.al2pv = {}
        self.pv2al = {}
        for (k, pv) in run.epicsinfo:
            if type(k) == tuple:
                al = k[0]
            else:
                al = k
            self.al2pv[al] = pv
            self.pv2al[pv] = al
        for p in PVlist:
            try:
                if type(p) == tuple:
                    # User specified PV and alias.
                    pv = p[0]
                    al = p[1]
                    self.pv2al[pv] = al
                    self.al2pv[al] = pv
                elif p in self.al2pv.keys():
                    # Known Alias
                    al = p
                    pv = self.al2pv[al]
                elif p in self.pv2al.keys():
                    # Known PV
                    pv = p
                    al = self.pv2al[pv]
                else:
                    # We don't know.  Assume it's a PV we'll find later.
                    pv = p
                    al = p
                    self.pv2al[pv] = al
                    self.al2pv[al] = pv
                self.pvs.append(run.Detector(pv))
                self.addPV(al, pv)
            except:
                print('could not find LCLS2 EPICS PV %s in data' % pv)
                self.missing.append(al)
                self.missingPV.append(pv)
        # Add these now so they don't interfere in the above loop.
        for (al, pv) in run.epicsinfo:
            self.al2pv[pv] = pv

    def addPV(self, al, pv):
        self.PVlist.append(al)

    def inRun(self):
        if len(self.pvs)>0:
            return True
        return False

    def data(self,evt):
        dl={}
        for pvname,pv in zip(self.PVlist,self.pvs):
            try:
                if pv(evt) is not None:
                    dl[pvname]=pv(evt)
                    if isinstance(dl[pvname], basestring):
                        dl[pvname]=np.nan
            except:
                #print('we have issues with %s in this event'%pvname)
                pass
        return dl

    def params_as_dict(self):
        """returns parameters as dictionary to be stored in the hdf5 file (once/file)"""
        parList =  {key:self.__dict__[key] for key in self.__dict__ if (key[0]!='_' and isinstance(getattr(self,key), (basestring, int, float, np.ndarray, tuple))) }
        PVlist = getattr(self,'PVlist')
        parList.update({'PV_%d'%ipv: pv for ipv,pv in enumerate(PVlist) if pv is not None})
        parList.update({'PVname_%d'%ipv: self.al2pv[pv] for ipv,pv in enumerate(PVlist) if pv is not None})
        return parList


class scanDetector(defaultDetector):
    def __init__(self, name='scan',run=None):
        self.name = name
        self.detname='scan'
        self.scans = []
        self.scanlist = []
        vetolist = ['step_docstring']
        try:
            scanlist = [k[0] for k in run.scaninfo if k[0] not in vetolist]
            for scan in scanlist:
                try:
                    self.scans.append(run.Detector(scan))
                    self.scanlist.append(scan)
                except:
                    print('could not find LCLS2 EPICS PV %s in data'%pv)
        except:
            pass

    def inRun(self):
        if len(self.scans)>0:
            return True
        return False

    def data(self,evt):
        dl={}
        for scanname,scan in zip(self.scanlist,self.scans):
            try:
                if scan(evt) is not None:
                    dl[scanname]=scan(evt)
                    if isinstance(dl[scanname], basestring):
                        dl[scanname]=np.nan
            except:
                #print('we have issues with %s in this event'%scanname)
                pass
        return dl
