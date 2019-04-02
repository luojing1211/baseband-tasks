# Licensed under the GPLv3 - see LICENSE
"""Full-package tests of phaseutils sources."""
import numpy as np
from numpy.testing import assert_allclose

import astropy.units as u
from astropy.time import Time

from baseband_tasks import pfb
from baseband_tasks.convolution import ConvolveSamples
from baseband_tasks.generators import NoiseGenerator


class TestBasics:
    def setup(self):
        self.nh = NoiseGenerator(shape=(32, 2048),
                                 start_time=Time('2010-01-01'),
                                 sample_rate=1.*u.kHz, seed=12345,
                                 samples_per_frame=8, dtype='f8')
        self.nc = NoiseGenerator(shape=(32, 2048),
                                 start_time=Time('2010-01-01'),
                                 sample_rate=1.*u.kHz, seed=12345,
                                 samples_per_frame=8, dtype='c16')
        self.pfb = pfb.sinc_hamming(4, 2048).reshape(4, 2048)

    def test_understanding(self):
        """Stepping or frequency selection should give same answer.

        Normally think of PFB as multiplying with an array that is, e.g.,
        4 times larger, then FFTing and taking every 4th frequency.

        But this is equivalent to, after multiplication, summing 4
        subsequent arrays and taking the FT of that.
        """
        self.nh.seek(0)
        # First check for real data.
        d = self.nh.read(5)
        hd = self.pfb * d[:4]
        ft1_hd = np.fft.rfft(hd.ravel())[::4]
        ft2_hd = np.fft.rfft(hd.sum(0))
        assert_allclose(ft1_hd, ft2_hd)

        # Check convolution gives the same answer,
        # remembering that for convolution we have to flip the order.
        ch = ConvolveSamples(self.nh, self.pfb[::-1])
        cd = ch.read(1)[0]
        assert_allclose(cd, hd.sum(0))
        ft_cd = np.fft.rfft(cd)
        assert_allclose(ft_cd, ft2_hd)
        # Next item just for completeness
        cd2 = ch.read(1)[0]
        ft_cd2 = np.fft.rfft(cd2)
        assert_allclose(ft_cd2, np.fft.rfft((self.pfb * d[1:]).sum(0)))

    def test_understanding_complex(self):
        # Check above holds for complex too.
        self.nc.seek(0)
        c = self.nc.read(4)
        hc = self.pfb * c
        ft1_hc = np.fft.fft(hc.ravel())[::4]
        ft2_hc = np.fft.fft(hc.sum(0))
        assert_allclose(ft1_hc, ft2_hc)

        ch = ConvolveSamples(self.nc, self.pfb[::-1])
        cd = ch.read(1)[0]
        ft_cd = np.fft.fft(cd)
        assert_allclose(ft_cd, ft2_hc)
