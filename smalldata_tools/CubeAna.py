import glob

class CubeAna(object):
    def __init__(self, expname='', run=-1, dirname='', filename='', plotWith='matplotlib'):
        print 'DEBUG INIT; ',expname, run
        self._fields={}
        self.expname=expname
        self.run=run
        self.runLabel='Run%03d'%run
        self.plotWith=plotWith
        if len(expname)>3:
            self.hutch=self.expname[:3]
            if dirname=='':
                self.dirname='/reg/d/psdm/%s/%s/hdf5/smalldata'%(self.hutch,self.expname)
                self.plot_dirname='/reg/d/psdm/%s/%s/results/smalldata_plots/'%(self.hutch,self.expname)
                #run 13 and past.
                if not path.isdir(self.dirname):
                    self.dirname='/reg/d/psdm/%s/%s/ftc'%(self.hutch,self.expname)
                    self.plot_dirname='/reg/d/psdm/%s/%s/res/smalldata_plots/'%(self.hutch,self.expname)
            else:
                self.dirname=dirname
                self.plot_dirname = dirname+'/smalldata_plots'
            if not path.isdir(self.plot_dirname):
                makedirs(self.plot_dirname)

        self.fname=''
        allFiles = glob.glob('%s/Cube*_%s_Run%03d.h5'%(self.dirname,self.expname,self.run))
        if filename != '':
            for thisFile in allFiles:
                if thisFile.find(filename)>=0:
                    self.fname = thisFile
                    break
            if self.fname == '':
                print 'Could not find a file with that name'
        if filename == '' or self.fname=='':
            cubeNames=[]
            for thisFile in allFiles:
                cubeNames.append(thisFile.replace('%s/Cube'%dirname,'').replace('_%s_Run%03d.h5'%(self.expname,self.run),''))
            if len(cubeNames)==0:
                print 'no cube files found, return'
                return none
            if len(cubeNames)==1:
                self.fname = '%s/Cube%s_%s_Run%03d.h5'%(self.dirname,cubeName[0],self.expname,self.run))
            else:
                print 'Options for Cube names are: ',cubeNames
                selectCubeName = raw_input('Select one of these options: ')
                self.fname = '%s/Cube%s_%s_Run%03d.h5'%(self.dirname,selectCubeName,self.expname,self.run))

        if path.isfile(self.fname):
            self.fh5=tables.open_file(self.fname,'r')
        else: #if path.isfile(self.fname):
            print 'could not find file: ',self.fname
            return None


    def Keys(self):
        return self.fh5.root.keys()

    def plotCube(self, sig=None, i0='nEntries'):
        return None
