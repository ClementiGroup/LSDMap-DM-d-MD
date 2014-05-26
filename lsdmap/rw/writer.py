import os

class WriterError(Exception):
    pass

class GroFormat(object):
    """
    class used to write .gro file(s)
    """
    def open(self, *args):
        import gro
        pattern = args[0]
        return gro.Writer(pattern)


class EvFormat(object):
    """
    class used to write .ev files
    """
    def open(self, *args):
        import ev
        return ev.Writer()

known_formats={'.gro': GroFormat(),
               '.ev': EvFormat()
}

class WriterFormatError(WriterError):

    def __init__(self, format, known_formats):
        self.format = format
        self.known_formats = known_formats

    def __str__(self):
        global known_formats
        return 'Unknown molecule file format "%s"\n Available formats are %s.' % (self.format, self.known_formats.keys())


def open(*args, **kargs):
    """
    *args corresponds to the file that will be used as a pattern to write the new data
    """

    if 'format' not in kargs:
        filename = args[0]
        format = os.path.splitext(filename)[1]
    else:
        filename = None
        format = kargs['format']

    if format not in known_formats:
        raise WriterFormatError(format, known_formats)

    return known_formats[format].open(filename)
