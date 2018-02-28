import numpy as np
from utilities_FitCenter import fitCircle

class fitCenter:
    def __init__(self, maskName='', threshold=90, name='', xArr=None, yArr=None, imgShape=None):
        if name!='':
            self.name=name
        else:
            self.name='center'
        self.threshold = threshold
        self.mask=None

        self.xArr = xArr
        self.yArr = yArr
        if self.xArr is None and (imgShape is None or len(imgShape)<2):
            print 'no geometry information is given, return'
            return
        if self.xArr is None:
            x = np.arange(0, imgShape[1])
            y = np.arange(0, imgShape[0])
            self.xArr,self.yArr = np.meshgrid(x,y)
        
        if maskName!='':
            self.mask=np.loadtxt(maskName)
            if self.xArr is not None and self.mask.shape != self.xArr.shape:
                if self.mask.shape[1] == 2:
                    self.mask = self.mask.transpose(1,0)
                self.mask = self.mask.reshape(self.xArr.shape)

    def getCenter(self,image):
        """
        select x/y array and fit circle to them. Return center & radius
    
        Parameters
        ----------
        image : np.ndarray
        The image or stack of images.
        Returns center & radius
        -------
        """
        if self.mask is not None:
            image = (image*self.mask)

        thresP = np.percentile(image, self.threshold)
        res = fitCircle(self.xArr.flatten()[image.flatten()>thresP],self.yArr.flatten()[image.flatten()>thresP])    
        res.pop('msg',None)
        res.pop('info',None)
        res.pop('success',None)
        res.pop('C',None)
        return res
