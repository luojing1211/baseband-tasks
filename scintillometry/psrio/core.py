"""psr_reader.py defines the classes for reading pulsar data from a non-baseband
format.
"""


from ..generators import StreamGenerator
from astropy.io import fits
from .psrfits_io import HDU_map
from astropy import log


__all__ = ['Reader', 'PsrfitsReader']

def open_read(filename, memmap=None):
    hdus = fits.open(filename, 'readonly', memmap=memmap)
    buffer = {'PRIMARY':[]}
    for ii, hdu in enumerate(hdus):
        if hdu.name in HDU_map.keys():
            if hdu.name in buffer.keys():
                buffer[hdu.name].append(hdu)
            else:
                buffer[hdu.name] = [hdu,]
        else:
             log.warn("The {}th HDU '{}'' is not a known PSRFITs"
                      " HDU.".format(ii, hdu.name))
    if len(buffer['PRIMARY']) > 1 or len(buffer['PRIMARY']) < 1:
        raise ValueError("File `{}` does not have a header"
                         " HDU or have more than one header"
                         " HDU.".format(filename))
    header_hdu = buffer['PRIMARY'][0]

    psrfits_hdus = []
    psrfits_hdus.append(HDU_map['PRIMARY'](header_hdu))
    buffer.pop('PRIMARY')
    for k, v in buffer.items():
        for hdu in v:
            psrfits_hdus.append(HDU_map[k](psrfits_hdus[0], hdu))
    return psrfits_hdus


class Reader(StreamGenerator):
    """Reader class defines the common API for the Read sub_class.

    Parameter
    ---------
    translator : `Translator` object or its sub_class
        The class for translating head and data from the other format
    """


    def __init__(self, source, function, **kwargs):
        self.source = source
        self.input_args = kwargs
        # The required argument will come from the source.
        self.req_args = {'shape': None, 'start_time': None,
                         'sample_rate': None}
        self.opt_args = {'samples_per_frame': 1, 'frequency': None,
                         'sideband': None, 'polarization': None, 'dtype':None}

        self._prepare_args()
        self._setup_args()

        super(Reader, self).__init__(*(function, self.req_args['shape'],
                                       self.req_args['start_time'],
                                       self.req_args['sample_rate']),
                                     **self.opt_args)

    def _prepare_args(self):
        """This setup function setups up the argrument for initializing the
        StreamGenerator.
        """
        input_args_keys = self.input_args.keys()
        source_properties = self.source._properties
        for rg in self.req_args.keys():
            if rg not in source_properties:
                raise ValueError("'{}' is required.".format(rg))

            self.req_args[rg] = getattr(self.source, rg)

        for og in self.opt_args.keys():
            if og in source_properties:
                self.opt_args[og] = getattr(self.source, og)
            elif og in input_args_keys:
                self.opt_args[og] = self.input_args[og]
            else:
                continue
        self._setup_args()
        return

    def _setup_args(self):
        pass

class HDUReader(Reader):
    """ This is a class for reading PSRFITS HDUs to scintillometry
    StreamGenerator style of file handleself.

    Parameter
    ---------
    psrfits_hdus: hdu object
        psrfits HDUs
    """
    def __init__(self, psrfits_hdu):
        super(HDUReader, self).__init__(psrfits_hdu, None)

    def _read_frame(self, frame_index):
        res = self.source.read_data_row(frame_index).T
        return res.reshape((self.samples_per_frame, ) + self.shape[1:])

    def _setup_args(self):
        # Reshape frequency.
        shape = self.req_args['shape']
        freq_shape = shape._replace(npol=1, nbin=1)[1:]
        self.opt_args['frequency'] = self.opt_args['frequency'].reshape(freq_shape)
