import numpy as np

from smalldata_tools.common.detector_base import DetObjectFunc


class pressioCompressDecompress(DetObjectFunc):
    """Function to compress and decompress the data.

    This leaves the data "unchanged" (hopefully). It can be used to evaluate
    the effect of compression algorithms on the output data, and ultimately
    downstream scientific results.

    It makes use of the `libpressio` package.
    """

    def __init__(self, **kwargs):
        # We'll hide the import of libpressio here - it is not available in conda envs
        from libpressio import PressioCompressor

        self._name = kwargs.get("name", "pressio")
        super().__init__(**kwargs)

        # For now we will just transparently pass the entire config for libpressio
        # via the prod_config file.
        self._compressor_json = kwargs.get("pressio_config")
        if self._compressor_json is None:
            raise RuntimeError("Libpressio compressor configuration must be provided!")

        self._compressor = PressioCompressor.from_config(self._compressor_json)

    def process(self, data):
        decompressed_img = np.zeros_like(data)

        # First compress image - output is immediately decompressed
        compressed_img = self._compressor.encode(data)

        decompressed_img = self._compressor.decode(compressed_img, decompressed_img)

        data = decompressed_img
        return {}
