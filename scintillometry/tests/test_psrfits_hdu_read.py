# Licensed under the GPLv3 - see LICENSE
"""Full-package tests of psrfits reading routine sources."""

import pytest
import numpy as np
import os
import astropy.units as u
from astropy.time import Time

from ..io import open_read

test_data = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')


class TestRead:
    def setup(self):
        self.chime_data = os.path.join(test_data, "CHIME_B1937+21_grid01_beam_7"
                                            "_58464_79458.fits")
        self.crab_data = os.path.join(test_data, "pulse_169370884.rf")


class TestCHIMERead(TestRead):
    def setup(self):
        super().setup()
        self.readers = open_read(self.chime_data)

    def test_shape(self):
        # The header HUD does not have any data
        shape0 = self.readers[0].shape
        assert shape0 == (0,)
        # Subint shape
        shape1 = self.readers[1].shape
        assert shape1 == (1, 256, 1, 4), \
            "CHIME data shape did not read correctly."

    def test_start_time(self):
        # Header start time
        start_time0 = self.readers[0].start_time
        frac, mjd = np.modf(start_time0.mjd)
        sec_frac, sec_int = np.modf(frac * 86400.0)
        assert np.isclose(self.readers[0].header['STT_IMJD'], mjd), \
            "The header HDU start time's integer MJD is not reading correctly."
        assert np.isclose(self.readers[0].header['STT_SMJD'], sec_int), \
            ("The header HDU start time's integer second is not reading"
            "correctly.")
        assert np.isclose(self.readers[0].header['STT_OFFS'], sec_frac), \
            ("The header HDU start time's fractional second is not reading"
            "correctly.")
        # Subint start time
        start_time1 = self.readers[1].start_time
        dt  = start_time1 - start_time0
        # There is some problem about dt.
        # sec_farc, sec_int = np.modf(frac * 86400.0)
        # assert np.isclose(self.readers[0].header['STT_IMJD'], mjd), \
        #     "The header HDU start time's integer MJD is not reading correctly."
        # assert np.isclose(self.readers[0].header['STT_SMJD'], sec_int), \
        #     ("The header HDU start time's integer second is not reading"
        #     "correctly.")
        # assert np.isclose(self.readers[0].header['STT_OFFS'], sec_frac), \
        #     ("The header HDU start time's fractional second is not reading"
        #     "correctly.")

    def test_read(self):
        # Header HDU does not have any data.
        with pytest.raises(EOFError) as excinfo:
            data0 = self.readers[0].read(1)
        assert "cannot read from beyond end of input." in str(excinfo.value)

        data1 = self.readers[1].read(1)
        assert data1.shape == (1,) + self.readers[1].shape[1:], \
            "The result data shape does not match the header reported shape."

       # TODO, I think a better way to test read is to compare with PSRCHIVE
       # result.


class TestCrabRead(TestRead):
    def setup(self):
        super().setup()
        self.readers = open_read(self.crab_data)

    def test_shape(self):
        # The header HUD does not have any data
        shape0 = self.readers[0].shape
        assert shape0 == (0,)
        # Subint shape
        shape1 = self.readers[1].shape
        assert shape1 == (1, 8192, 20, 4), \
            "CHIME data shape did not read correctly."

    def test_start_time(self):
        # Header start time
        start_time0 = self.readers[0].start_time
        frac, mjd = np.modf(start_time0.mjd)
        sec_frac, sec_int = np.modf(frac * 86400.0)
        assert np.isclose(self.readers[0].header['STT_IMJD'], mjd), \
            "The header HDU start time's integer MJD is not reading correctly."
        assert np.isclose(self.readers[0].header['STT_SMJD'], sec_int), \
            ("The header HDU start time's integer second is not reading"
            "correctly.")
        assert np.isclose(self.readers[0].header['STT_OFFS'], sec_frac), \
            ("The header HDU start time's fractional second is not reading"
            "correctly.")
        # Subint start time
        start_time1 = self.readers[1].start_time
        dt  = start_time1 - start_time0
        # There is some problem about Subint start time.
        # sec_farc, sec_int = np.modf(frac * 86400.0)
        # assert np.isclose(self.readers[0].header['STT_IMJD'], mjd), \
        #     "The header HDU start time's integer MJD is not reading correctly."
        # assert np.isclose(self.readers[0].header['STT_SMJD'], sec_int), \
        #     ("The header HDU start time's integer second is not reading"
        #     "correctly.")
        # assert np.isclose(self.readers[0].header['STT_OFFS'], sec_frac), \
        #     ("The header HDU start time's fractional second is not reading"
        #     "correctly.")

    def test_read(self):
        # Header HDU does not have any data.
        with pytest.raises(EOFError) as excinfo:
            data0 = self.readers[0].read(1)
        assert "cannot read from beyond end of input." in str(excinfo.value)

        data1 = self.readers[1].read(1)
        assert data1.shape == (1,) + self.readers[1].shape[1:], \
            "The result data shape does not match the header reported shape."
