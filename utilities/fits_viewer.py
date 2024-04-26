import astrodata
import numpy as np
import matplotlib.pyplot as plt
from gempy.gemini import gemini_tools as gt
import sys

def view_extensions(adinput):
    """
    View the extensions of a FITS file.

    Parameters
    ----------
    filename : str
        The name of the FITS file to view.
    """
    log = gt.logutils.get_logger(__name__)
    adobj = astrodata.open(adinput)
    log.stdinfo("Filename: {}".format(adobj.filename))
    log.stdinfo(adobj.info())
    print(adobj.info())

def peaks(adinput):
    log = gt.logutils.get_logger(__name__)
    log.stdinfo("Opening peaks extension")
    adobj = astrodata.open(adinput)
    peaks = adobj[0].PEAKS
    print(peaks)


def main(args):
    print(args)
    view_extensions(args[1])
    peaks(args[1])

if __name__ == "__main__":
    main(args=sys.argv)

