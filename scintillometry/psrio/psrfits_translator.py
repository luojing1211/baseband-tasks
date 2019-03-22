"""psrfits_translator.py defines the translator between baseband-style file
object and psrfits file hdus.
"""

from .translator import TranslatorBase
from collections import namedtuple
import astropy.units as u


__all__ = [HeaderHDU, SubintHDU]


class HeaderHDU(TranslatorBase):
    """HeaderHDU class provides the translator function between baseband-style
    file object and the PSRFITS main header HDU.

    Parameter
    ---------
    header_hdu : `Pdat.header_hdu`
    """
    _properties = ('start_time', 'observatory', 'filename')

    def __init__(self, header_hdu):
        self.header_hdu = header_hdu
        self.verify(self.header_hdu)
        super(HeaderHDU, self).__init__('baseband', 'psrfits_header')

    def verify(self, input_hdu):
        try:
            self.header =  input_hdu.read_header()
            assert self.header['SIMPLE']
        except:
            raise ValueError("Input HDU is not a main header HDU.")

    @property
    def start_time(self):
        MJD_d_int = self.header['STT_IMJD'] * u.day
        MJD_s_int = self.header['STT_SMJD'] * u.s
        MJD_s_frac = self.header['STT_OFFS'] * u.s
        #NOTE I am assuming PSRFITS's time scale is UTC
        return Time(MJD_d_int, MJD_s_int + MJD_s_frac, format='mjd',
                    scale='utc', precision=9)

    @start_time.setter
    def start_time(self, time):
        if not isinstance(time, Time):
            raise ValueError("Input time should be an `astropy.time.Time`"
                             " object.")
        if time.isisscalar:
            raise ValueError("Input time should be a scalar")
        else:
            self.header['STT_IMJD'] = int(divmod(time.mjd, 1)[0])
            day_int = Time(self.header['STT_IMJD'], scale=time.scale,
                           format=time.format)
            dt = (time - day_int).to(u.s)
            int_sec, frac_sec = divmod(dt.value, 1)
            self.header['STT_SMJD'] = int(int_sec)
            self.header['STT_OFFS'] = frac_sec


    def get_filename(self):
        return self.header_hdu.get_filename()


class SubintHDU(TranslatorBase):
    """SubintHDU class provides the translator functions between baseband-style
    file object and the PSRFITS SUBINT HDU.

    Parameter
    ---------
    name : str
        The name of the translator instance.
    header_hdu : object
        The psrfits main header object
    subint_hdu : HDU oject
        The psrfits data HDU.

    Note
    ----
    Right now we are assuming the data rows are continuous in time and the
    frequency are the same.
    """

    _properties = ('start_time', 'sample_rate', 'shape', 'pol', 'frequency')

    _defaults = [('BACKEND', 'GUPPI'),
                 ('BLOCSIZE', 0),
                 ('STT_OFFS', 0),
                 ('PKTIDX', 0),
                 ('OVERLAP', 0),
                 ('SRC_NAME', 'unset'),
                 ('TELESCOP', 'unset'),
                 ('PKTFMT', '1SFA'),
                 ('PKTSIZE', 8192),
                 ('NBITS', 8),
                 ('NPOL', 1),
                 ('OBSNCHAN', 1)]

    def __init__(self, header_hdu, subint_hdu):
        self.header_hdu = HeaderTranslator("file_header", header_hdu)
        self.data_hdu = data_hdu
        self.data_header = self.data_hdu.read_header()
        super(SubintHDU, self).__init__('baseband', 'SUBINT')
        self.verify(self.data_hdu)

    def verify(self, input_hdu):
        if input_hdu.get_filename() != self.header_hdu.get_filename():
            raise ValueError("Main header HDU and input HDU are not from the "
                             "same file.")
        else:
            if input_hdu.get_extname() != self.source:
                raise ValueError("Input HDU is not a SUBINT type.")
            else:
                return

    @property
    def sample_rate(self):
        return 1.0 / (self.data_header['TBIN'] * u.s)

    @sample_rate.setter
    def sample_rate(self, val):
        if hasattr(val, 'unit'):
            val = val.to(u.Hz).value
        self.data_hdu['TBIN'] = 1.0 / val

    @property
    def shape(self):
        nrows = self.data_hdu.get_nrows()
        samples_per_row = self.data_header['NSBLK']
        nchan = self.data_header['NCHAN']
        npol = self.data_header['NPOL']
        nbin = self.data_header['NBIN']
        return (nrows * samples_per_row, nbin, npol, nchan)

    def get_dim_label(self):
        return ('time', 'phase', 'pol', 'freq')

    @property
    def pol(self):
        return self.data_header['POL_TYPE']

    @pol.setter
    def pol(self, val):
        self.data_header['POL_TYPE'] = val

    @property
    def start_times(self):
        # NOTE should we get the start time for each raw, in case the time gaps
        # in between the rows
        file_start = self.header_hdu.get_start_time()
        subint_times = u.Quantity(self.header_hdu.read_column(col='OFFS_SUB'),
                                  u.s, copy=False)
        samples_per_row = self.data_header['NSBLK']
        sample_time = self.data_header['TBIN']
        start_time = (file_start + subint_times -
                      samples_per_row / 2 * sample_time)
        return start_time

    @start_time.setter
    def start_time(self):
        pass

    def get_freqs(self):
        # Those are the frequency for all the rows.
        freqs = u.Quantity(self.header_hdu.read_column(col='DAT_FREQ'),
                           u.MHz, copy=False)
        return freqs

    def get_sideband(self):
        # It is not clear for now.
        return 1

    def get_data(self, time_samples):
        # The seek is not working
        samples_per_row = self.data_header['NSBLK']
        num_rows = int(time_samples / samples_per_row)
        rest_samples = time_samples - num_rows * samples_per_row
        data = self.data_hdu.read(row=np.arange(num_rows))
        result = data['DATA'].reshape((num_rows * samples_per_row,
                                       self.get_shape[1::]))
        return result[0: time_samples]
