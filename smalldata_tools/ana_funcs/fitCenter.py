import numpy as np
from smalldata_tools.utilities_FitCenter import fitCircle
from smalldata_tools.DetObject import DetObjectFunc


class fitCenter(DetObjectFunc):
    def __init__(self, **kwargs):
        try:
            self._name = 'center_'+kwargs.get('_name')
        except:
            self._name = 'center'
        self.threshold = kwargs.get('threshold',90)

        self._xArr = kwargs.get('xArr',None)
        self._yArr = kwargs.get('yArr',None)
        imgShape = kwargs.get('imgShape', None)
        if self._xArr is None and imgShape is not None:
            x = np.arange(0, imgShape[1])
            y = np.arange(0, imgShape[0])
            self._xArr,self._yArr = np.meshgrid(x,y)

        userMask = kwargs.get('userMask',None) #fill these from detector object. Unless mask is set from userMask
        self._mask = None
        if userMask is not None:
            self._mask=np.loadtxt(userMask)
            if self._xArr is not None and self._mask.shape != self._xArr.shape:
                if self._mask.shape[1] == 2: #what is this for - seems to be turning a cs140 mask around?
                    self._mask = self._mask.transpose(1,0)
                self._mask = self._mask.reshape(self._xArr.shape)

    def setFromDet(self, det):
        super(fitCenter, self).setFromDet(det)
        if self._mask is None and det.mask is not None:
            setattr(self, '_mask', det.mask.astype(np.uint8))
        if det.x is not None:
            setattr(self, '_xArr', det.x)
        if det.y is not None:
            setattr(self, '_yArr', det.y)

    #def getCenter(self,image):
    def process(self,image):
        """
        select x/y array and fit circle to them. Return center & radius
    
        Parameters
        ----------
        image : np.ndarray
        The image or stack of images.
        Returns center & radius
        -------
        """
        if self._mask is not None:
            image = (image*self._mask)

        if self._xArr is None:
            x = np.arange(0, image.shape[1])
            y = np.arange(0, image.shape[0])
            self._xArr,self._yArr = np.meshgrid(x,y)

        thresP = np.percentile(image, self.threshold)
        res = fitCircle(self._xArr.flatten()[image.flatten()>thresP],self._yArr.flatten()[image.flatten()>thresP])    
        res.pop('msg',None)
        res.pop('info',None)
        res.pop('success',None)
        res.pop('C',None)
        return res
