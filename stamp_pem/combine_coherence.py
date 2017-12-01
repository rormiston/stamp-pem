#!/usr/bin/python
import h5py
from gwpy.spectrum import Spectrum
import optparse
from stamp_pem import coh_io
from stamp_pem.coherence_segment import PEMSubsystem
from gwpy.segments import (DataQualityFlag, Segment, SegmentList)
import numpy as np


def combine_coherence(darm_channel, subsystem, st, et, directory, jobdur,
        verbose=True):
    """
    combine coherence between two times

    Parameters
    ----------
    darm_channel : `str`
        differential arm channel
    subsystem : `str`
        subsystem to combine data for
    st : `int`
        start time for combining data
    et : `int`
        end time for combining data
    jobdur : `int`
        duration of analyzed data saved in files

    Returns
    -------
    subsys : `stamp_pem.coherence_segment.PEMCoherenceSubsystem`
        combined data for subsystem
    """
    segs = SegmentList()
    info_str = ""
    ii = 0
    while st + (ii+1)*jobdur <= et:
        segs.append(Segment([st + (ii*params.jobdur), st + ((ii+1)*params.jobdur)]))
        ii + 1
    First = 1
    subsys=None
    for seg in segs:
        seg_st = seg[0]
        seg_et = seg[1]
        cohdir = coh_io.get_directory_structure(subsystem, seg_st, directory=directory)
        cohfile = coh_io.create_coherence_data_filename(darm_channel, subsystem, seg_st,
                                                        seg_et, directory=cohdir)
        if First:
            try:
                subsys = PEMSubsystem.read(subsystem, cohfile)
                First = 0
            except IOError:
                info_str+=("Couldn't load %s\n" % cohfile)
                continue
        else:
            try:
                temp = PEMSubsystem.read(subsystem, cohfile)
            except IOError:
                info_str+=("Couldn't load %s\n" % cohfile)
                continue
            subsys.update(temp)
    if verbose:
        print info_str
    return subsys


