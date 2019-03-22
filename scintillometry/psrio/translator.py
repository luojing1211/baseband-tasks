"""translator.py defines a the base class for translating data between two file
format.
"""

class TranslatorBase:
    """The base translator class which defines set of common API functions for
    translatort container subclasses.
    Parameter
    ---------
    name : str
        The name of translator instance.
    source : str
        The name of the input format.
    target : str
        The name of the target format.
    """
    def __init__(self, file_type, source_name, client_name, mode):
        self.file_type = file_type
        self.source_name = source_name
        self.client_name = client_name
        self.mode = mode

    def verify(self, input_object):
        """ This function checks if the input file object matches the translator
        source.
        """
        raise NotImplementError
