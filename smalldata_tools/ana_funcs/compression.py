from smalldata_tools.common.detector_base import DetObjectFunc


class CompressDecompress(DetObjectFunc):

    def __init__(self, config, **kwargs):
        self._name = kwargs.get("name", "compression")
        super().__init__(**kwargs)

        self._config = config

        def setFromDet(self, det):
            super().setFromDet(det)

        def process(self, data):
            """
            Implement the compression / decompression here
            """
            data = self._compress(data)
            data = self._decompress(data)

            # self.dat will be used by any function chained after this one as the data input
            self.dat = data

            # Any returned dictionary will be saved to the h5 file under this detector
            return {}
        
        def _compress(self, data):
            return data
        
        def _decompress(self, data):
            return data

