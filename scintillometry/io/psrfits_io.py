"""psrfits_translator.py defines the translator between baseband-style file
object and psrfits file hdus.
"""

from .psrfits_fields import *
from collections import namedtuple
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle, Latitude, Longitude
from astropy.io import fits
import numpy as np


__all__ = ["PsrfitsHearderHDU", "SubintHDU", "HDU_map"]


class PsrfitsHearderHDU(fits.PrimaryHDU):
    """HeaderHDU class provides the translator function between baseband-style
    file object and the PSRFITS main header HDU.

    Parameter
    ---------
    header_hdu : `Pdat.header_hdu`

    NOTE
    ----
    the frequency marker is on the center of the channels
    """
    _properties = ('start_time', 'observatory', 'frequency', 'ra', 'dec',
                   'shape', 'sample_rate')

    _header_defaults = main_header

    def __init__(self, header_hdu=None):
        if header_hdu is None:
            super(PsrfitsHearderHDU, self).__init__()
        else:
            super(PsrfitsHearderHDU, self).__init__(header=header_hdu.header,
                                                    data=header_hdu.data)
            self.verify()

    def verify(self):
        if 'SIMPLE' not in list(self.header.keys()):
            raise ValueError("The input HDU is not a fits headers HDU.")
        if self.header['FITSTYPE'] != "PSRFITS":
            raise ValueError("The input fits header is not a PSRFITS type.")

    @property
    def start_time(self):
        # TODO add location
        return (Time(self.header['STT_IMJD'], format='mjd', precision=9) +
                TimeDelta(self.header['STT_SMJD'], self.header['STT_OFFS'],
                          format='sec', scale='tai'))

    @start_time.setter
    def start_time(self, time):
        # TODO add location
        if not isinstance(time, Time):
            raise ValueError("Input time should be an `astropy.time.Time`"
                             " object.")
        if not time.isscalar:
            raise ValueError("Input time should be a scalar")

        mjd_int = int(time.mjd)
        mjd_int = Time(self.header['STT_IMJD'], scale=time.scale,
                       format=time.format)
        mjd_frac = (time - Time(mjd_int, scale=time.scale,
                                format=time.format)).jd
        if mjd_frac < 0:
            mjd_int -= 1
            mjd_frac += 1.
        mjd_frac = (mjd_frac * u.day).to(u.s)
        int_sec = int(mjd_frac.value)
        frac_sec = mjd_frac.value - int_sec
        self.header['STT_IMJD'] = '{0:05d}'.format(mjd_int)
        self.header['STT_SMJD'] = '{}'.format(int_sec)
        self.header['STT_OFFS'] = '{0:17.15f}'.format(frac_sec)
        self.header['DATE-OBS'] = time.fits

    @property
    def observatory(self):
        return self.header['TELESCOP']

    @observatory.setter
    def observatory(self, val):
        self.header['TELESCOP'] = val
        #TODO ADD ITRF coordinate.

    @property
    def frequency(self):
        try:
            n_chan = float(self.header['OBSNCHAN'])
            c_chan = float(self.header['OBSFREQ'])
            bw = float(self.header['OBSBW'])
        except Exception:
            return None
        chan_bw = bw / n_chan
        freq = np.arange(-n_chan / 2, n_chan / 2) * chan_bw + c_chan
        return u.Quantity(freq, u.MHz, copy=False)

    @frequency.setter
    def frequency(self, frequency):
        #NOTE this assumes the frequency marker is on the center of the channels
        frequency = u.Quantity(frequency, u.MHz, copy=False)
        n_chan = len(frequency)
        self.header['OBSNCHAN'] = n_chan
        self.header['OBSBW'] = (frequency[-1] - frequency[0]).value
        self.header['OBSFREQ'] = frequency[int(n_chan / 2)].value

    @property
    def ra(self):
        return Longitude(self.header['RA'], unit=u.hourangle)

    @ra.setter
    def ra(self, val):
        val = Longitude(val, unit=u.hourangle)
        self.header['RA'] = val.to_string(sep=':')

    @property
    def dec(self):
        return Latitude(self.header['DEC'], unit=u.deg)

    @ra.setter
    def dec(self, val):
        val = Latitude(val, unit=u.deg)
        self.header['DEC'] = val.to_string(sep=':')

    @property
    def obs_mode(self):
        return self.header['OBS_MODE']


class SubintHDUBase(fits.BinTableHDU):
    """SubintHDU class provides the translator functions between baseband-style
    file object and the PSRFITS SUBINT HDU.

    Parameter
    ---------
    header_hdu : PsrfitsHearderHDU
        The psrfits main header object
    subint_hdu : HDU object
        The psrfits data HDU.

    Note
    ----
    Right now we are assuming the data rows are continuous in time and the
    frequency are the same.
    """

    _properties = ('start_time', 'sample_rate', 'shape', 'samples_per_frame',
                   'polarization', 'frequency')
    _header_defaults = subint_header
    _req_columns = subint_columns

    def __init__(self, header_hdu, subint_hdu=None):
        self.header_hdu = header_hdu
        if subint_hdu is None:
            super(SubintHDU, self).__init__(name='SUBINT')
        else:
            super(SubintHDU, self).__init__(header=subint_hdu.header,
                                            data=subint_hdu.data)
        self.verify()
        self.offset = 0

    def verify(self):
        if self.header['EXTNAME'].strip() != "SUBINT":
            raise ValueError("Input HDU is not a SUBINT type.")
        # Check frequency
        if 'DAT_FREQ' in self.columns.names:
            freqs = u.Quantity(self.data['DAT_FREQ'],
                               u.MHz, copy=False)
            assert np.isclose(freqs, freqs[0],
                              atol=np.finfo(float).eps * u.MHz).all(), \
                "Frequencies are different within one subint rows."

    @property
    def mode(self):
        return self.header_hdu.obs_mode

    @property
    def sample_rate(self):
        # We assume the sample rate are the same for all the rows.
        sample_time = u.Quantity(self.data[0]['TSUBINT'] /
                                 self.samples_per_frame, u.s)
        return 1.0 / sample_time

    @sample_rate.setter
    def sample_rate(self, val):
        self.header['TBIN'] = (1.0 / val).to_value(u.s)

    @property
    def start_time(self):
        # NOTE should we get the start time for each raw, in case the time gaps
        # in between the rows
        file_start = self.header_hdu.start_time
        if "OFFS_SUB" in self.columns.names:
            subint_times = u.Quantity(self.data['OFFS_SUB'], u.s, copy=False)
            sample_time = 1.0 / self.sample_rate
            start_time = (file_start + subint_times[0] -
                          self.samples_per_frame / 2 * sample_time)
        else:
            start_time = file_start
        return start_time

    @start_time.setter
    def start_time(self, time):
        # NOTE this sets the start time of the HDU, not the file start time.
        dt = (time - self.header_hdu.start_time).to(u.s)
        center_off = dt + 1.0 / self.sample_rate * self.samples_per_frame / 2
        if "OFFS_SUB" in self.columns.names:
            first_off = self.data['OFFS_SUB'][0]
            self.data['OFFS_SUB'] = (self.data['OFFS_SUB'] - first_off +
                                     center_off)
        else:
            self._req_columns['OFFS_SUB']['value'] = center_off

    @property
    def nrow(self):
        return int(self.header['NAXIS2'])

    @property
    def nchan(self):
        return int(self.header['NCHAN'])

    @property
    def npol(self):
        return int(self.header['NPOL'])

    @property
    def nbin(self):
        return int(self.header['NBIN'])

    @property
    def samples_per_frame(self):
        try:
            return int(self.header['NSBLK'])
        except:
            return int(np.prod(self.data_shape)/(self.nrow * self.nchan *
                                                 self.npol * self.nbin))

    @samples_per_frame.setter
    def samples_per_frame(self, val):
        self.header['NSBLK'] = val

    @property
    def shape(self):
        raw_shape = self.raw_shape
        new_shape = namedtuple('shape', ['nsample','nbin', 'nchan', 'npol'])
        result = new_shape(raw_shape.nrow * raw_shape.samples_per_frame,
                           raw_shape.nbin, raw_shape.nchan, raw_shape.npol)
        return result

    @property
    def raw_shape(self):
        r_shape = namedtuple('shape', ['nrow', 'samples_per_frame','nbin',
                                       'nchan','npol'])
        result = r_shape(self.nrow, self.samples_per_frame, self.nbin,
                         self.nchan, self.npol)
        return result

    @property
    def data_shape(self):
        # Data are save in the fortran order. Reversed from the header label.
        d_shape_raw = self.data['DATA'].shape
        if self.mode == "SEARCH":
            d_shape_header = (self.nbin, self.nchan, self.npol,
                              self.samples_per_frame)
        else:
            d_shape_header = (self.nbin, self.nchan, self.npol)
        if d_shape_raw != (self.nrow, ) + d_shape_header[::-1]:
            raise ValueError("Data shape does not match with the header"
                             " information")
        if self.mode == "SEARCH":
            d_shape = namedtuple('d_shape', ['nrow', 'nsample', 'npol', 'nchan',
                                             'nbin'])
            result  = d_shape(self.nrow,  self.samples_per_frame, self.npol,
                              self.nchan, self.nbin)
        else:
            d_shape = namedtuple('d_shape', ['nrow', 'npol', 'nchan', 'nbin'])
            result  = d_shape(self.nrow,  self.npol, self.nchan, self.nbin)
        return result

    @property
    def polarization(self):
        return self.header['POL_TYPE']

    @polarization.setter
    def polarization(self, val):
        self.header['POL_TYPE'] = val

    @property
    def frequency(self):
        if 'DAT_FREQ' in self.columns.names:
            freqs = u.Quantity(self.data['DAT_FREQ'],
                               u.MHz, copy=False)[0]
        else:
            # Get the frequency from the header HDU
            try:
                freqs = self.header_hdu.frequency
            except:
                freqs = None
        return freqs

    @frequency.setter
    def frequency(self, val):
        if 'DAT_FREQ' in self.columns.names:
            # Check frequency length.
            assert self.shape.nchan == len(val), \
                "Input frequency does not match the number of channels in data."
            self.data['DAT_FREQ'] = val

        else:
            self._req_columns['DAT_FREQ']['value'] = val

    @property
    def dtype(self):
        return self.data['DATA'].dtype

    def read_data_row(self, row_index):
        if row_index >= self.shape[0]:
            raise EOFError("cannot read from beyond end of input SUBINT HDU.")

        row = self.data[row_index]
        # Reversed the header shape to match the data
        new_shape = self.raw_shape._replace(samples_per_frame=1,
                                            nbin=1)[-1:1:-1]
        data_scale = row['DAT_SCL'].reshape(new_shape)
        data_off_set = row['DAT_OFFS'].reshape(new_shape)
        zero_off = np.zeros((), dtype=self.dtype)
        try:
            zero_off = np.asarray(self.header['ZERO_OFF'], dtype=self.dtype)
        except:
            pass
        result = ((row['DATA'] - zero_off)* data_scale +
                   data_off_set)
        if 'DAT_WTS' in self.columns.names:
            data_wts = row['DAT_WTS']
            wts_shape = self.raw_shape._replace(samples_per_frame=1, nbin=1,
                                                npol=1)[-1:1:-1]
            result *= data_wts.reshape(wts_shape)
        return result


class SearchSubint(SubintHDUBase):
    """SearchSubint class is designed for handling the searching mode PSRFITS
    Subint HDU.

    Parameter
    ---------
    header_hdu : PsrfitsHearderHDU
        The psrfits main header object
    search_hdu : HDU object
        The psrfits HDU.
    """
    def __init__(self, header_hdu, search_hdu=None):
        super().__init__(header_hdu, search_hdu)
        self.verify()

    def verify(self):
        if self.header_hdu.obs_mode.upper() != 'SEARCH':
            raise ValueError("Header HDU is not in the searching mode.")
        try:
            # Check NSBLK field
            nsblk = int(self.header['NSBLK'])
        except:
            raise ValueError("Searching mode requires a 'NSBLK' field in the"
                             " header.")
        if nsblk <= 1:
            raise ValueError("Subint HDU is not in the searching mode.")

        try:
            # Check TBIN field
            nsblk = int(self.header['TBIN'])
        except:
            raise ValueError("Searching mode requires a 'TBIN' field in the"
                             " header.")

    @property
    def start_time(self):
        file_start = self.header_hdu.start_time
        if "OFFS_SUB" in self.columns.names:
            subint_times = u.Quantity(self.data['OFFS_SUB'], u.s, copy=False)
            sample_time = 1.0 / self.sample_rate
            start_time = (file_start + subint_times[0] -
                          self.samples_per_frame * sample_time / 2)
        else:
            start_time = file_start
        return start_time

    @start_time.setter
    def start_time(self, time):
        # NOTE this sets the start time of the HDU, not the file start time.
        dt = (time - self.header_hdu.start_time).to(u.s)
        center_off = dt + 1.0 / self.sample_rate * self.samples_per_frame / 2
        if "OFFS_SUB" in self.columns.names:
            first_off = self.data['OFFS_SUB'][0]
            self.data['OFFS_SUB'] = (self.data['OFFS_SUB'] - first_off +
                                     center_off)
        else:
            self._req_columns['OFFS_SUB']['value'] = center_off

    @property
    def samples_per_frame(self):
        return int(self.header['NSBLK'])

    @property
    def sample_rate(self):
        sample_time = u.Quantity(self.['TBIN'], u.s)
        return 1.0 / sample_time

    @property
    def data_shape(self):
        # Data are save in the fortran order. Reversed from the header label.
        d_shape_raw = self.data['DATA'].shape
        d_shape_header = (self.nbin, self.nchan, self.npol,
                          self.samples_per_frame)
        if d_shape_raw != (self.nrow, ) + d_shape_header[::-1]:
            raise ValueError("Data shape does not match with the header"
                             " information")
        d_shape = namedtuple('d_shape', ['nrow', 'nsample', 'npol', 'nchan',
                                         'nbin'])
        result  = d_shape(self.nrow,  self.samples_per_frame, self.npol,
                          self.nchan, self.nbin)
        return result


class PSRSubint(SubintHDUBase):
    """PSRSubint class is designed for handling the pulsar folding mode PSRFITS
    Subint HDU.

    Parameter
    ---------
    header_hdu : PsrfitsHearderHDU
        The psrfits main header object
    psr_subint : HDU object
        The psrfits subint HDU.
    """
    def __init__(self, header_hdu, psr_subint=None):
        super().__init__(header_hdu, psr_subint)
        self.verify()

    def verify(self):
        if self.header_hdu.obs_mode.upper() != 'PSR':
            raise ValueError("Header HDU is not in the folding mode.")
        try:
            # Check NSBLK field
            nbin = int(self.header['NBIN'])
        except:
            raise ValueError("Folding mode requires a 'NBIN' field in the"
                             " header.")
        if nbin <= 1:
            raise ValueError("Subint HDU is not in the folding mode.")

    @property
    def start_time(self):
        file_start = self.header_hdu.start_time
        if "OFFS_SUB" in self.columns.names:
            subint_times = u.Quantity(self.data['OFFS_SUB'], u.s, copy=False)
            sample_time = 1.0 / self.sample_rate
            start_time = (file_start + subint_times[0] -
                          self.samples_per_frame * sample_time / 2)
        else:
            start_time = file_start
        return start_time

    @start_time.setter
    def start_time(self, time):
        # NOTE this sets the start time of the HDU, not the file start time.
        dt = (time - self.header_hdu.start_time).to(u.s)
        center_off = dt + 1.0 / self.sample_rate * self.samples_per_frame / 2
        if "OFFS_SUB" in self.columns.names:
            first_off = self.data['OFFS_SUB'][0]
            self.data['OFFS_SUB'] = (self.data['OFFS_SUB'] - first_off +
                                     center_off)
        else:
            self._req_columns['OFFS_SUB']['value'] = center_off

    @property
    def samples_per_frame(self):
        return int(np.prod(self.data_shape)/(self.nrow * self.nchan *
                                             self.npol * self.nbin))

    @property
    def sample_rate(self):
        sample_time = u.Quantity(self.['TBIN'], u.s)
        return 1.0 / sample_time

    @property
    def data_shape(self):
        # Data are save in the fortran order. Reversed from the header label.
        d_shape_raw = self.data['DATA'].shape
        d_shape_header = (self.nbin, self.nchan, self.npol)
        if d_shape_raw != (self.nrow, ) + d_shape_header[::-1]:
            raise ValueError("Data shape does not match with the header"
                             " information")
        d_shape = namedtuple('d_shape', ['nrow', 'npol', 'nchan', 'nbin'])
        result  = d_shape(self.nrow,  self.npol, self.nchan, self.nbin)
        return result


HDU_map = {'PRIMARY': PsrfitsHearderHDU,
           'SUBINT': SubintHDU}
