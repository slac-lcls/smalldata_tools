import numpy as np

from smalldata_tools.common.detector_base import DetObjectFunc


class PolynomialCurveCorrection(DetObjectFunc):
    """
    A class for applying polynomial curve corrections to detector images.
    This class applies a polynomial-based correction to shift columns (or rows) of an image
    according to a polynomial function. This is particularly useful for correcting geometric
    distortions in detector images where the distortion can be described by a polynomial.

    Parameters
    ----------
    polynomial_coefficients : list, optional
        Coefficients of the polynomial function that describes the correction curve.
        Coefficients are in descending order of degree (highest degree first).
    axis : int, default=0
        The axis along which to apply the correction:
        - 0: correction applied to columns (image is transposed first)
        - 1: correction applied to columns directly
    method : str, default='roll'
        The method used for applying the correction:
        - 'roll': basic implementation using np.roll (works for all cases)
        - 'vectorized': vectorized implementation (fastest, only for monotonic corrections)
        - 'ndimage': uses scipy.ndimage for non-integer shifts (slowest)
    name : str, optional
        Name of the correction function, used for identification in the UserDataCfg.

    Example
    -------
    >>> coeffs = [0.001, -0.1, 0]  # quadratic curve
    >>> corrector = PolynomialCurveCorrection(polynomial_coefficients=coeffs)
    >>> corrected_data = corrector.process(image_data)
    """

    def __init__(
        self,
        polynomial_coefficients: list = None,
        axis: int = 0,
        method="roll",
        **kwargs,
    ):
        self._name = kwargs.get("name", "curve_corr")
        super().__init__(**kwargs)
        self.polynomial_coefficients = polynomial_coefficients
        self.axis = axis
        self.method = method

    def setFromDet(self, det):
        """
        Just add some potentially useful attributes from the detector object
        to the function object, so they are saved under the UserDataCfg and can be used later.
        """
        if det.mask is not None:
            self.mask = self.process(det.mask)
        if det.cmask is not None:
            self.cmask = self.process(det.cmask)

    def process(self, data):
        if self.axis == 0:
            data = np.transpose(data)

        if self.method == "roll":
            shifted_image = self.polynomial_correction(
                data, self.polynomial_coefficients
            )
        elif self.method == "vectorized":
            shifted_image = self.polynomial_correction_vectorized(
                data, self.polynomial_coefficients
            )
        elif self.method == "ndimage":
            shifted_image = self.curve_correction_ndimage(
                data, self.polynomial_coefficients
            )

        if self.axis == 0:
            shifted_image = np.transpose(shifted_image)

        self.dat = shifted_image
        return {"corrected_image": shifted_image}

    @staticmethod
    def polynomial_correction(image, polynomial_coefficients):
        """
        Shifts columns of the image according to a polynomial function.
        """
        px = np.arange(0, image.shape[1])
        shift = np.polyval(polynomial_coefficients, px)
        shift = shift - np.min(shift)
        shift = np.round(shift).astype(np.int64)
        image_shifted = np.zeros_like(image)
        for i in range(0, image.shape[1]):
            image_shifted[:, i] = np.roll(image[:, i], -shift[i])
        return image_shifted

    @staticmethod
    def polynomial_correction_vectorized(data, poly_coefficients):
        """
        Most efficient but only works if the polynomial correction is monotonic.
        """
        rows, cols = image.shape
        px = np.arange(cols)
        shift = np.polyval(polynomial_coefficients, px)
        shift = shift - np.min(shift)
        shift = np.round(shift).astype(np.int64)

        # Create row indices for each column
        row_indices = np.arange(rows)[:, np.newaxis]

        # Calculate the shifted row indices for each column
        shifted_rows = (row_indices - shift[np.newaxis, :]) % rows

        # Create column indices
        col_indices = np.arange(cols)[np.newaxis, :]

        # Create the shifted image using advanced indexing
        return image[shifted_rows, col_indices]

    @staticmethod
    def curve_correction_ndimage(image, poly_coefficients):
        """
        Slowest, but might be useful if more generic correction are to be implemented
        in the future, or if we want to support non-integer shifts.
        """
        px = np.arange(0, image.shape[1])
        shift = np.polyval(polynomial_coefficients, px)
        shift = shift - np.min(shift)

        # Create output array
        output = np.zeros_like(image)

        # Apply shift to each column - still uses a loop but with optimized ndimage function
        for i in range(image.shape[1]):
            output[:, i] = ndimage_shift(
                image[:, i], -shift_values[i], mode="constant", cval=0, order=0
            )

        return output
