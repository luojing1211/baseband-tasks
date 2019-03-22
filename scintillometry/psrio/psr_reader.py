"""psr_reader.py defines the classes for reading pulsar data from a non-baseband
format.
"""


from ..Generator import StreamGenerator


__all__ = ['Reader', 'PsrfitsReader']

class Reader(StreamGenerator):
    """Reader class defines the common API for the Read sub_class.

    Parameter
    ---------
    translator : `Translator` object or its sub_class
        The class for translating head and data from the other format
    """


    def __init__(self, translator, **kwargs):
        self.translator = translator
        self.args = {'function': self.read_format_data}
        self.args.update(kwargs)
        self.required_args = ['shape', 'start_time', 'sample_rate']
        self.optional_args = ['samples_per_frame', 'frequency', 'sideband',
                              'polarization', 'dtype']
        self._prepare_args()
        super(Reader, self).__init__(**self.args)

    def _prepare_args(self):
        """This setup function setups up the argrument for initializing the
        StreamGenerator.
        """
        input_args_keys = self.args.keys()
        translatort_keys = self.translator.keys()
        for rg in self.required_args:
            if rg not in translatort_keys and rg not in input_args_keys:
                raise ValueError("'{}' is required. You can input it while "
                                 "initialization or give a function in the "
                                 "translator.")
            self.args[rg] = self.translator[rg]()

        for og in self.required_args:
            if og in translatort_keys and og not in input_args_keys:
                self.args[og] = self.translator[og]()
