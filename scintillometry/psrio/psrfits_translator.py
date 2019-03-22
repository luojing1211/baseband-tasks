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

    _defaults = [('HDRVER', '6.1 ' 'Header version')
                 ('FITSTYPE', 'PSRFITS ', 'FITS definition for pulsar data'
                  'files'),
                 ('DATE', ' ', 'File creation UTC date (YYYY-MM-DDThh:mm:ss)'),
                 ('OBSERVER', ' ', 'Observer name(s)'),
                 ('PROJID', ' ', 'Project name'),
                 ('TELESCOP', ' ', 'Telescope name'),
                 ('ANT_X', 0  '[m] Antenna ITRF X-coordinate (D)'),
                 ('ANT_Y', 0  '[m] Antenna ITRF Y-coordinate (D)'),
                 ('ANT_Z', 0  '[m] Antenna ITRF Z-coordinate (D)'),
                 ('FRONTEND', ' ', 'Receiver ID'),
                 ('IBEAM', ' ', 'Beam ID for multibeam systems'),
                 ('NRCVR', 0, 'Number of receiver polarisation channels'),
                 ('FD_POLN', ' ', 'LIN or CIRC'),
                 ('FD_HAND', 0, '+/- 1. +1 is LIN:A=X,B=Y, CIRC:A=L,B=R (I)'),
                 ('FD_SANG', 0, '[deg] FA of E vect for equal sig in A&B (E)'),
                 ('FD_XYPH', 0, '[deg] Phase of A* B for injected cal (E)'),
                 ('BACKEND', ' ', 'Backend ID'),
                 ('BECONFIG', ' ', 'Backend configuration file name'),
                 ('BE_PHASE', 0, '0/+1/-1 BE cross-phase:0 unknown,+/-1 std/rev'),
                 ('BE_DCC' , 0, '0/1 BE downconversion conjugation corrected'),
                 ('BE_DELAY', 0, '[s] Backend propn delay from digitiser input'),
                 ('TCYCLE' , 0, '[s] On-line cycle time (D)'),
                 ('OBS_MODE', ' ', '(PSR, CAL, SEARCH)'),
                 ('DATE-OBS', ' ', 'UTC date of observation (YYYY-MM-DDThh:mm:ss)'),
                 ('OBSFREQ', 0, '[MHz] Centre frequency for observation'),
                 ('OBSBW', 0, '[MHz] Bandwidth for observation'),
                 ('OBSNCHAN', 0, 'Number of frequency channels (original)'),
                 ('CHAN_DM', 0, '[cm-3 pc] DM used for on-line dedispersion'),
                 ('PNT_ID', ' ', 'Name or ID for pointing ctr (multibeam feeds)'),
                 ('SRC_NAME', ' ', 'Source or scan ID'),
                 ('COORD_MD', ' ', 'Coordinate mode (J2000, GALACTIC, ECLIPTIC)'),
                 ('EQUINOX', 0, 'Equinox of coords (e.g. 2000.0)'),
                 ('RA', ' ', 'Right ascension (hh:mm:ss.ssss)'),
                 ('DEC', ' ', 'Declination (-dd:mm:ss.sss)'),
                 ('BMAJ', 0, '[deg] Beam major axis length'),
                 ('BMIN', 0, '[deg] Beam minor axis length'),
                 ('BPA' , 0, '[deg] Beam position angle'),
                 ('STT_CRD1', ' ', 'Start coord 1 (hh:mm:ss.sss or ddd.ddd)'),
                 ('STT_CRD2', ' ', 'Start coord 2 (-dd:mm:ss.sss or -dd.ddd)'),
                 ('TRK_MODE', ' ', 'Track mode (TRACK, SCANGC, SCANLAT)'),
                 ('STP_CRD1', ' ', 'Stop coord 1 (hh:mm:ss.sss or ddd.ddd)'),
                 ('STP_CRD2', ' ', 'Stop coord 2 (-dd:mm:ss.sss or -dd.ddd)'),
                 ('SCANLEN' ,0, '[s] Requested scan length (E)'),
                 ('FD_MODE', ' ', 'Feed track mode - FA, CPA, SPA, TPA'),
                 ('FA_REQ', 0, '[deg] Feed/Posn angle requested (E)'),
                 ('CAL_MODE', ' ', 'Cal mode (OFF, SYNC, EXT1, EXT2)'),
                 ('CAL_FREQ', 0, '[Hz] Cal modulation frequency (E)'),
                 ('CAL_DCYC', 0, 'Cal duty cycle (E)'),
                 ('CAL_PHS', 0, 'Cal phase (wrt start time) (E)'),
                 ('CAL_NPHS', 0, 'Number of states in cal pulse (I)'),
                 ('STT_IMJD', 0, 'Start MJD (UTC days) (J - long integer)'),
                 ('STT_SMJD', 0, '[s] Start time (sec past UTC 00h) (J)'),
                 ('STT_OFFS', 0, '[s] Start time offset (D)'),
                 ('STT_LST', 0, '[s] Start LST (D)')]

    def __init__(self, header_hdu, mode):
        self.header_hdu = header_hdu
        self.verify(self.header_hdu)
        super(HeaderHDU, self).__init__('psrfits', 'baseband',
                                        'psrfits_header', mode)

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


    _buffer = [('OFFS_SUB', 0),]

    def __init__(self, header_hdu, subint_hdu, mode):
        self.header_hdu = HeaderTranslator("file_header", header_hdu)
        self.data_hdu = data_hdu
        self.data_header = self.data_hdu.read_header()
        super(SubintHDU, self).__init__('psrfits','baseband', 'SUBINT', mode)
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
    def samples_per_frame(self):
        return self.data_header['NSBLK']

    @samples_per_frame.setter
    def samples_per_frame(self, val):
        self.data_header['NSBLK'] = val

    @property
    def start_times(self):
        # NOTE should we get the start time for each raw, in case the time gaps
        # in between the rows
        file_start = self.header_hdu.start_time
        subint_times = u.Quantity(self.data_hdu.read_column(col='OFFS_SUB'),
                                  u.s, copy=False)
        samples_per_row = self.data_header['NSBLK']
        sample_time = 1.0 / self.sample_rate
        start_time = (file_start + subint_times[0] -
                      self.samples_per_frame / 2 * sample_time)
        return start_time

    @start_time.setter
    def start_time(self, time):
        # NOTE this sets the start time of the HDU, not the file start time.
        dt = (time - self.header_hdu.start_time.to(u.s)
        center_off = dt + 1.0 / self.sample_rate * self.samples_per_frame / 2
        self._buffer['OFFS_SUB'] = center_off


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
