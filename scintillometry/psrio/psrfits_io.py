"""psrfits_translator.py defines the translator between baseband-style file
object and psrfits file hdus.
"""

from .psrfits_fields import *
from collections import namedtuple
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle, Latitude, Longitude
from astropy.io import fits
import numpy as np


__all__ = ["PsrfitsHearderHDU", "SubintHDU"]


class PsrfitsHearderHDU(fits.PrimaryHDU):
    """HeaderHDU class provides the translator function between baseband-style
    file object and the PSRFITS main header HDU.

    Parameter
    ---------
    header_hdu : `Pdat.header_hdu`
    """
    _properties = ('start_time', 'observatory', 'frequency', 'source')

    _header_defaults = main_header

    def __init__(self, header_hdu=None):
        if header_hdu is None:
            super(PsrfitsHearderHDU, self).__init__()
            self._init_empty()
        else:
            super(PsrfitsHearderHDU, self).__init__(header=header_hdu.header,
                                                    data=header_hdu.data)
            self.verify()

    def verify(self):
        if 'SIMPLE' not in list(self.header.keys()):
            raise ValueError("The input HDU is not a fits headers HDU.")
        if self.header['FITSTYPE'] != "PSRFITS":
            raise ValueError("The input fits header is not a PSRFITS type.")

    def _init_empty(self):
        for k in self._defaults:
            self.header.set(k[0], k[1], k[2])
        self.header.add_comment('FITS (Flexible Image Transport System) format'
                                'is defined in Astronomy and Astrophysics, '
                                'volume 376, page 359; bibcode:'
                                ' 2001A&A...376..359H', after='EXTEND')

    @property
    def start_time(self):
        # TODO add location
        MJD_d_int = self.header['STT_IMJD'] * u.day
        MJD_s_int = self.header['STT_SMJD'] * u.s
        MJD_s_frac = self.header['STT_OFFS'] * u.s
        #NOTE I am assuming PSRFITS's time scale is UTC
        return Time(MJD_d_int, MJD_s_int + MJD_s_frac, format='mjd',
                    scale='utc', precision=9)

    @start_time.setter
    def start_time(self, time):
        # TODO add location
        if not isinstance(time, Time):
            raise ValueError("Input time should be an `astropy.time.Time`"
                             " object.")
        if not time.isscalar:
            raise ValueError("Input time should be a scalar")
        else:
            self.header['STT_IMJD'] = int(divmod(time.mjd, 1)[0])
            day_int = Time(self.header['STT_IMJD'], scale=time.scale,
                           format=time.format)
            dt = (time - day_int).to(u.s)
            int_sec, frac_sec = divmod(dt.value, 1)
            self.header['STT_SMJD'] = int(int_sec)
            self.header['STT_OFFS'] = frac_sec
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
        n_chan = self.header['OBSNCHAN']
        c_chan = self.header['OBSFREQ']
        bw = self.header['OBSBW']
        chan_bw = bw / n_chan
        #NOTE the frequency marker is on the center of the channels
        freq = np.arange(-n_chan / 2, n_chan / 2) * chan_bw + c_chan
        return u.Quantity(freq, u.MHz, copy=False)

    @frequency.setter
    def frequency(self, frequency):
        #NOTE this assumes the frequency marker is on the center of the channels
        if hasattr(frequency, 'unit'):
            frequency = frequency.to(u.MHz)
        else:
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


class SubintHDU(fits.BinTableHDU):
    """SubintHDU class provides the translator functions between baseband-style
    file object and the PSRFITS SUBINT HDU.

    Parameter
    ---------
    header_hdu : PsrfitsHearderHDU
        The psrfits main header object
    subint_hdu : HDU oject
        The psrfits data HDU.

    Note
    ----
    Right now we are assuming the data rows are continuous in time and the
    frequency are the same.
    """

    _properties = ('start_time', 'sample_rate', 'shape', 'samples_per_frame',
                   'polarization', 'frequency')
    _headr_defaults = subint_header
    _req_columns = subint_columns

    def __init__(self, header_hdu, subint_hdu=None):
        self.header_hdu = header_hdu
        if subint_hdu is None:
            super(SubintHDU, self).__init__(name='SUBINT')
            self._init_empty()
        else:
            super(SubintHDU, self).__init__(header=subint_hdu.header,
                                            data=subint_hdu.data)
        self.verify()
        self.offset = 0

    def _init_empty(self):
        for k in self._defaults:
            self.header.set(k[0], k[1], k[2])

    def make_columns(self, name):
        cols = []
        for cl in self._columns_defaults:
            col = fits.Column()

    def verify(self):
        if self.header['EXTNAME'].replace(' ', '') != "SUBINT":
            raise ValueError("Input HDU is not a SUBINT type.")

    @property
    def mode(self):
        return self.header_hdu.obs_mode

    @property
    def sample_rate(self):
        sample_time = u.Quantity(self.data[0]['TSUBINT'] / self.nrow /
                                 self.samples_per_frame, u.s)
        return 1.0 / sample_time

    @sample_rate.setter
    def sample_rate(self, val):
        if hasattr(val, 'unit'):
            val = val.to(u.Hz).value
        self.header['TBIN'] = 1.0 / val

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
                               u.MHz, copy=False)
            assert np.isclose(freqs, freqs[0],
                              atol=np.finfo(float).eps * u.MHz).all(), \
                "Frequencies are different within one subint rows."
            freqs = freqs[0]
        else:
            freqs = None
        return freqs

    @property
    def dtype(self):
        return self.data['DATA'].dtype

    @frequency.setter
    def frequency(self, val):
        if 'DAT_FREQ' in self.columns.names:
            self.data['DAT_FREQ'] = val
            ## Check frequency length.
            assert self.shape.nchan == len(val), \
                "Input frequency does not match the number of channels in data."
        else:
            self._req_columns['DAT_FREQ']['value'] = val

    def read_data_row(self, row_index):
        if row_index >= self.shape[0]:
            raise EOFError("cannot read from beyond end of input SUBINT HDU.")

        row = self.data[row_index]
        # Reversed the header shape to match the data
        new_shape = self.raw_shape._replace(samples_per_frame=1,
                                            nbin=1)[-1:1:-1]
        data_scale = row['DAT_SCL'].reshape(new_shape)
        data_off_set = row['DAT_OFFS'].reshape(new_shape)
        zero_off = np.asarray(0, dtype=self.dtype)
        if 'ZERO_OFF' in self.header.keys():
            try:
                zero_off = np.asarray(self.header['ZERO_OFF'], dtype=self.dtype)
            except:
                pass
        result = ((row['DATA'] - zero_off)* data_scale +
                   data_off_set)
        data_wts = np.ones(self.nchan)
        if 'DAT_WTS' in self.columns.names:
            data_wts = row['DAT_WTS']
        wts_shape = self.raw_shape._replace(samples_per_frame=1, nbin=1,
                                            npol=1)[-1:1:-1]
        result *= data_wts.reshape(wts_shape)
        return result


HDU_map = {'PRIMARY': PsrfitsHearderHDU,
           'SUBINT': SubintHDU}
