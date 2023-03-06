import warnings

# https://stackoverflow.com/a/26433913
def _warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return "%s:%s: %s: %s\n" % (filename, lineno, category.__name__, message)


def make_warnings_nicer():
    """
    Get rid of the annoying second line of warnings.
    """
    warnings.formatwarning = _warning_on_one_line


def remove_string_end(s: str, end: str):
    assert s.endswith(end)
    return s[: -len(end)]

class UnsupportedError(Exception):
    pass
