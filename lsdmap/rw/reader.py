import os

class ReaderError(Exception):
    pass

class GroFormat(object):
    """
    class used to read .gro file(s)
    """
    def open(self, filename, kwargs):
        import gro
        return gro.Reader(filename)

class XvgFormat(object):
    """
    class used to read .xvg file(s)
    """
    def open(self, filename, kwargs):
        import xvg
        return xvg.Reader(filename, **kwargs)

class SlFormat(object):
    """
    class used to read files containing one single line (.eps, .w)
    """
    def open(self, filename, kwargs):
        import sl
        return sl.Reader(filename)

class EvFormat(object):
    """
    class used to read .ev files
    """
    def open(self, filename, kwargs):
        import ev
        return ev.Reader(filename)

known_formats={'.gro': GroFormat(),
               '.xvg': XvgFormat(),
               '.eps': SlFormat(),
               '.w': SlFormat(),
               '.ev': EvFormat() 
}

class ReaderFormatError(ReaderError):

    def __init__(self, format, known_formats):
        self.format = format
        self.known_formats = known_formats

    def __str__(self):
        global known_formats
        return 'Unknown molecule file format "%s"\n Available formats are %s.\n' % (self.format, self.known_formats.keys())

def open(filename, **kargs):

    if isinstance(filename, (list, tuple)):
        nfiles = len(filename)
        if nfiles == 1:
            filename = filename[0]
    elif isinstance(filename, basestring):
        nfiles = 1

    if nfiles > 1:
        formats = [os.path.splitext(name)[1] for name in filename]
        format = formats[0]
        if not all(format == format_q for format_q in formats):
            raise ReaderError('filenames provided should have the same format')
    elif nfiles == 1:
         format = os.path.splitext(filename)[1]

    if format not in known_formats:
        raise ReaderFormatError(format, known_formats)

    return known_formats[format].open(filename, kargs)
