"""core.py defines the classes for reading pulsar data from a non-baseband
format.
"""


from ..generators import StreamGenerator
from ..base import BaseTaskBase
from astropy.io import fits
from .psrfits_io import HDU_map
from astropy import log
from collections import defaultdict

__all__ = ['open_read', 'Reader', 'HDUReader']


def open_read(filename, memmap=None):
    """ Function to read one PSRFITS file into a list of HDU Readers.

    Parameter
    ---------
    filename : str
        File name of the input PSRFITS file
    memmap : bool, optional
        Is memory mapping to be used? This value is obtained from the
        configuration item astropy.io.fits.Conf.use_memmap. Default is True.

    Return
    ------
    A list of the HDU readers.
    """
    hdus = fits.open(filename, 'readonly', memmap=memmap)
    buffer = defaultdict(list)
    for ii, hdu in enumerate(hdus):
        if hdu.name in HDU_map.keys():
            buffer[hdu.name].append(hdu)
        else:
            log.warn("Skipping HDU {} ({}), as it is not a known PSRFITs"
                      " HDU.".format(ii, hdu.name))

    primary = buffer.pop('PRIMARY', [])
    if len(primary) != 1:
        raise ValueError("File `{}` does not have a header"
                         " HDU or have more than one header"
                         " HDU.".format(filename))
    primary_hdu = HDU_map['PRIMARY'](primary[0])
    psrfits_hdus = []
    psrfits_hdus.append(primary_hdu)
    for k, v in buffer.items():
        for hdu in v:
            psrfits_hdus.append(HDU_map[k](primary_hdu, hdu))

    # Build reader on the HDUs
    readers = []
    for hdus in psrfits_hdus:
        readers.append(HDUReader(hdus))
    return readers


class Reader(StreamGenerator):
    """Reader class defines the common API for the Read sub_class.

    Parameter
    ---------
    scource : object
        The source object for the input file.
    """
    _req_args = {'shape': None, 'start_time': None, 'sample_rate': None}
    _opt_args = {'samples_per_frame': 1, 'frequency': None, 'sideband': None,
                 'polarization': None, 'dtype':None}
    def __init__(self, source, function, **kwargs):
        self.source = source
        self.input_args = kwargs
        self._prepare_args()
        self._setup_args()

        super().__init__(*(function, self._req_args['shape'],
                                       self._req_args['start_time'],
                                       self._req_args['sample_rate']),
                                     **self._opt_args)

    def _prepare_args(self):
        """This setup function setups up the argrument for initializing the
        StreamGenerator.
        """
        input_args_keys = self.input_args.keys()
        source_properties = self.source._properties
        for rg in self._req_args.keys():
            try:
                self._req_args[rg] = getattr(self.source, rg)
            except AttributeError as exc:
                exc.args += ("souce file should define '{}'.".format(rg),)
                raise exc

        for og in self._opt_args.keys():
            if og in input_args_keys:
                self._opt_args[og] = self.input_args[og]
            else:
                try:
                    self._opt_args[og] = getattr(self.source, og)
                except AttributeError:
                    pass
        self._setup_args()
        return

    def _setup_args(self):
        pass


class HDUReader(Reader):
    """ Wrapper to treat PSRFITS HDUs as streams.

    Parameter
    ---------
    psrfits_hdus: hdu
        psrfits HDUs
    """
    def __init__(self, psrfits_hdu):
        super().__init__(psrfits_hdu, None)

    def _read_frame(self, frame_index):
        res = self.source.read_data_row(frame_index).T
        return res.reshape((self.samples_per_frame, ) + self.sample_shape)

    def _setup_args(self):
        # Reshape frequency.
        shape = self._req_args['shape']
        if shape != (0, ):
            #TODO Need to be more generic here.
            freq_shape = shape._replace(npol=1, nbin=1)[1:]
            self._opt_args['frequency'] = self._opt_args['frequency'].reshape(freq_shape)
        else:
            self._opt_args['frequency']= None

    @property
    def header(self):
        return self.source.header
