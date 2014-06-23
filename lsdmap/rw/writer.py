import os

class WriterError(Exception):
    pass

class GroFormat(object):
    """
    class used to write .gro file(s)
    """
    def open(self, kwargs):
        import gro
        return gro.Writer(**kwargs)


class EvFormat(object):
    """
    class used to write .ev files
    """
    def open(self, **kwargs):
        import ev
        return ev.Writer(kwargs)

class SlFormat(object):
    """
    class used to write files with a single column
    """
    def open(self, *kwargs):
        import sl
        return sl.Writer()


known_formats={'.gro': GroFormat(),
               '.ev': EvFormat(),
               '.eps': SlFormat(),
               '.w': SlFormat()
}

class WriterFormatError(WriterError):

    def __init__(self, format, known_formats):
        self.format = format
        self.known_formats = known_formats

    def __str__(self):
        global known_formats
        return 'Unknown molecule file format "%s"\n Available formats are %s.' % (self.format, self.known_formats.keys())


def open(format, **kwargs):
    """
    *args corresponds to the file that will be used as a pattern to write the new data

    Examples
    --------

    writer.open('.eps') will create a sl.writer instance
    writer.open('.gro', pattern='aladip.gro') will create a gro.writer instance with file pattern aladip.gro

    """

    if format not in known_formats:
        raise WriterFormatError(format, known_formats)

    return known_formats[format].open(kwargs)
