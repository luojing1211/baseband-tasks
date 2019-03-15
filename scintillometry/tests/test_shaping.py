# Licensed under the GPLv3 - see LICENSE

import pytest
import numpy as np
from numpy.testing import assert_array_equal
import astropy.units as u

from ..shaping import (Reshape, Transpose, ReshapeAndTranspose,
                       ChangeSampleShape, GetItem)

from baseband import vdif, dada
from baseband.data import SAMPLE_VDIF, SAMPLE_DADA


def get_fh():
    """Get sample VDIF file with correct frequency, sideband, polarization."""
    fh = vdif.open(SAMPLE_VDIF)
    # Add frequency, sideband, and polarization information by hand.
    fh.frequency = 311.25 * u.MHz + (np.arange(8.) // 2) * 16. * u.MHz
    fh.sideband = 1
    fh.polarization = np.tile(['L', 'R'], 4)
    return fh


class TestReshape:
    @pytest.mark.parametrize('sample_shape', ((4, 2), (2, 4)))
    def test_reshape(self, sample_shape):
        fh = vdif.open(SAMPLE_VDIF)
        ref_data = fh.read().reshape((-1,) + sample_shape)

        rt = Reshape(fh, sample_shape=sample_shape)
        assert fh.sample_shape == (8,)
        assert rt.sample_shape == sample_shape
        assert rt.start_time == fh.start_time
        assert rt.sample_rate == fh.sample_rate

        data = rt.read()
        assert_array_equal(data, ref_data)
        rt.close()

    def test_frequency_sideband_polarization_propagation1(self):
        fh = get_fh()
        with Reshape(fh, (4, 2)) as rt:
            assert rt.frequency.shape == (4, 1)
            assert np.all(rt.frequency == fh.frequency[::2].reshape(4, 1))
            assert rt.sideband.shape == ()
            assert np.all(rt.sideband == 1)
            assert rt.polarization.shape == (2,)
            assert np.all(rt.polarization == fh.polarization[:2])

    def test_frequency_sideband_polarization_propagation2(self):
        fh = vdif.open(SAMPLE_VDIF)
        # Add different frequency, sideband, and polarization information.
        # (Note: these are incorrect; just for testing purposes.)
        fh.frequency = 311.25 * u.MHz + (np.arange(8.) // 4) * 16. * u.MHz
        fh.sideband = np.tile([-1, 1], 4)
        fh.polarization = np.tile(['L', 'L', 'R', 'R'], 2)
        with Reshape(fh, (2, 2, 2)) as rt:
            assert rt.frequency.shape == (2, 1, 1)
            assert np.all(rt.frequency == fh.frequency[::4].reshape(2, 1, 1))
            assert rt.sideband.shape == (2,)
            assert np.all(rt.sideband == fh.sideband[:2])
            assert rt.polarization.shape == (2, 1)
            assert np.all(rt.polarization == fh.polarization[:4:2].reshape(2, 1))

    def test_wrong_shape(self):
        with vdif.open(SAMPLE_VDIF) as fh:
            with pytest.raises(ValueError):
                Reshape(fh, (4, 4))
            with pytest.raises(ValueError):
                Reshape(fh, ())


class TestTranspose:
    """Test transpose on a reshaped stream."""
    @staticmethod
    def get_reshape_and_transpose(fh, sample_shape=(4, 2),
                                  sample_axes=(2, 1)):
        rt = Reshape(fh, sample_shape=sample_shape)
        return Transpose(rt, sample_axes=sample_axes)

    def test_basic(self):
        fh = vdif.open(SAMPLE_VDIF)
        ref_data = fh.read().reshape((-1, 4, 2)).transpose(0, 2, 1)
        tt = self.get_reshape_and_transpose(fh, (4, 2), (2, 1))
        assert tt.start_time == fh.start_time
        assert tt.sample_rate == fh.sample_rate
        data = tt.read()
        assert_array_equal(data, ref_data)
        tt.close()

    def test_frequency_sideband_polarization_propagation1(self):
        fh = get_fh()
        with self.get_reshape_and_transpose(fh, (4, 2), (2, 1)) as tt:
            assert tt.frequency.shape == (4,)
            assert np.all(tt.frequency == fh.frequency[::2])
            assert tt.sideband.shape == ()
            assert np.all(tt.sideband == 1)
            assert tt.polarization.shape == (2, 1)
            assert np.all(tt.polarization == fh.polarization[:2].reshape(2, 1))

    def test_frequency_sideband_polarization_propagation2(self):
        fh = vdif.open(SAMPLE_VDIF)
        # Add different frequency, sideband, and polarization information.
        # (Note: these are incorrect; just for testing purposes.)
        fh.frequency = 311.25 * u.MHz + (np.arange(8.) // 4) * 16. * u.MHz
        fh.sideband = np.tile([-1, 1], 4)
        fh.polarization = np.tile(['L', 'L', 'R', 'R'], 2)
        with self.get_reshape_and_transpose(fh, (2, 2, 2), (-1, -3, -2)) as tt:
            assert tt.frequency.shape == (2, 1)
            assert np.all(tt.frequency == fh.frequency[::4].reshape(2, 1))
            assert tt.sideband.shape == (2, 1, 1)
            assert np.all(tt.sideband == fh.sideband[:2].reshape(2, 1, 1))
            assert tt.polarization.shape == (2,)
            assert np.all(tt.polarization == fh.polarization[:4:2])

    def test_wrong_axes(self):
        with vdif.open(SAMPLE_VDIF) as fh:
            with pytest.raises(ValueError):
                self.get_reshape_and_transpose(fh, (4, 2), (1, 0))
            with pytest.raises(ValueError):
                self.get_reshape_and_transpose(fh, (4, 2), (2, 3))
            with pytest.raises(ValueError):
                self.get_reshape_and_transpose(fh, (4, 2), (1, 1))


class TestReshapeAndTranspose(TestTranspose):
    """Test reshape and transpose combination directly.

    Same tests as TestTranspose.
    """
    get_reshape_and_transpose = ReshapeAndTranspose


class TestChangeSampleShape(TestTranspose):
    """Test custom shaping using a reshape and transpose function.

    TestTranspose tests as well as some additional ones.
    """
    @staticmethod
    def get_reshape_and_transpose(fh, sample_shape=(4, 2),
                                  sample_axes=(2, 1)):
        new_shape = (-1,) + sample_shape
        new_axes = (0,) + sample_axes

        def task(data):
            return data.reshape(new_shape).transpose(new_axes)

        return ChangeSampleShape(fh, task)

    def test_swap_axes(self):
        fh = get_fh()
        st = ChangeSampleShape(
            fh, lambda data: data.reshape(-1, 4, 2).swapaxes(1, 2))
        assert st.frequency.shape == (4,)
        assert np.all(st.frequency == fh.frequency[::2])
        assert st.sideband.shape == ()
        assert np.all(st.sideband == 1)
        assert st.polarization.shape == (2, 1)
        assert np.all(st.polarization == fh.polarization[:2].reshape(2, 1))
        st.close()

    def test_get_item(self):
        """Selecting from both axes of two-dimensional samples."""
        fh = get_fh()
        with ChangeSampleShape(
                fh, lambda data: data.reshape(-1, 4, 2)[:, :2, 0]) as st:
            assert st.frequency.shape == (2,)
            assert_array_equal(st.frequency, fh.frequency[:4:2])
            assert st.polarization.shape == ()
            assert_array_equal(st.polarization, fh.polarization[0])

    def test_no_extra_arguments(self):
        with pytest.raises(TypeError):
            ChangeSampleShape(None, None, shape=())


class TestGetItem:

    @pytest.mark.parametrize('item', (0, slice(0, None, 2), [0, 1]))
    def test_basic(self, item):
        """Basic tests on one-dimensional samples."""
        fh = get_fh()
        ref_data = fh.read()[:, item]
        gih = GetItem(fh, item)
        data = gih.read()
        assert_array_equal(data, ref_data)
        ref_freq = fh.frequency[item]
        if ref_freq.size > 1 and np.all(ref_freq[0] == ref_freq):
            ref_freq = ref_freq[0]
        ref_sideband = fh.sideband
        ref_pol = fh.polarization[item]
        if ref_pol.size > 1 and np.all(ref_pol[0] == ref_pol):
            ref_pol = ref_pol[0, ...]

        assert_array_equal(gih.frequency, ref_freq)
        assert_array_equal(gih.sideband, ref_sideband)
        assert_array_equal(gih.polarization, ref_pol)
        gih.close()

    def test_specific1(self):
        """Selecting one item in first axis of two-dimensional samples."""
        rh = Reshape(get_fh(), (4, 2))  # freq, pol
        with GetItem(rh, 0) as gih:
            assert gih.frequency.shape == ()
            assert_array_equal(gih.frequency, rh.frequency[0, 0])
            assert gih.polarization.shape == (2,)
            assert_array_equal(gih.polarization, rh.polarization)

    def test_specific2(self):
        """Selecting one item in second axis of two-dimensional samples."""
        rh = Reshape(get_fh(), (4, 2))  # freq, pol
        with GetItem(rh, (slice(None), 1)) as gih:
            assert gih.frequency.shape == (4,)
            assert_array_equal(gih.frequency, rh.frequency.squeeze())
            assert gih.polarization.shape == ()
            assert_array_equal(gih.polarization, rh.polarization[1])

    def test_specific3(self):
        """Selecting from both axes of two-dimensional samples."""
        rh = Reshape(get_fh(), (4, 2))  # freq, pol
        with GetItem(rh, (slice(None, 2), 0)) as gih:
            assert gih.frequency.shape == (2,)
            assert_array_equal(gih.frequency, rh.frequency[:2].squeeze())
            assert gih.polarization.shape == ()
            assert_array_equal(gih.polarization, rh.polarization[0])

    def test_wrong_item(self):
        with vdif.open(SAMPLE_VDIF) as fh:
            with pytest.raises(IndexError):
                GetItem(fh, 10)
            with pytest.raises(IndexError):
                GetItem(fh, (1, 1))
