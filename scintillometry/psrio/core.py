"""psr_reader.py defines the classes for reading pulsar data from a non-baseband
format.
"""


from ..generators import StreamGenerator
from astropy.io import fits
from .psrfits_io import HDU_map


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
             warn("HDU {} is not a PSRFITs HDUs.".format(ii))
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
        self.args = {'function': function}
        self.args.update(kwargs)
        self.required_args = ['shape', 'start_time', 'sample_rate']
        self.optional_args = ['samples_per_frame', 'frequency', 'sideband',
                              'polarization', 'dtype']
        self._prepare_args()
        super(StreamGenerator, self).__init__(*self.args)

    def _prepare_args(self):
        """This setup function setups up the argrument for initializing the
        StreamGenerator.
        """
        input_args_keys = self.args.keys()
        source_properties = self.source._properties
        for rg in self.required_args:
            if rg not in source_properties and rg not in input_args_keys:
                raise ValueError("'{}' is required. You can input it while "
                                 "initialization or give a function in the "
                                 "translator.")
            self.args[rg] = getattr(self.source, rg)
            print(rg, self.args[rg])

        for og in self.optional_args:
            if og in source_properties and og not in input_args_keys:
                self.args[og] = getattr(self.source, og)
                print(og, self.args[og])

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
        return self.source.read_data_row(frame_index)
