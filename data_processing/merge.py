__author__ = "davidmurphy"


"""
A function used collapse a set of overlapping segments and return a new set of into non-overlapping merged segments
"""


def merge(segments):
    """
    A generator that yields new merged segments by collapsing overlapping segments (the union).
    :param segments: segment list 1 in format [(start1, end1) ... (startN, endN)]
    Example:
    --------
    segments = [(1, 10), (7, 12), (15, 30), (25, 55), (45, 51)]
    result   = list(merge(segments))
    result   = [(1, 12), (15, 55)]
    """

    # store number of segments
    nsegs = len(segments)

    # segments need to be pre-sorted by segment starts but they may overlap
    assert all([segments[i][0] <= segments[i+1][0] for i in xrange(nsegs-1)])

    # initialize the current start and stop with the first segment
    (current_start, current_stop) = segments[0]

    # process the remaining segments merging where necessary
    for (start, stop) in segments[1:]:

        # new segment disjoint with previous: yield previous segment
        if start > current_stop:
            yield (current_start, current_stop)
            (current_start, current_stop) = (start, stop)

        # new segment overlaps previous segment: merge previous+new
        else:
            current_stop = max(stop, current_stop)

    # yield the final segment and stop iteration
    yield (current_start, current_stop)


def merge_weighted(segments, values):
    """
    A generator that yields new merged segments by collapsing overlapping
    segments (the union) and taking the average of values corresponding to
    segments weighted by segment length.
    :param segments: segment list 1 in format [(start1, end1) ... (startN, endN)]
    :param values: a list of values corresponding to segments
    Example:
    --------
    segments = [(1, 10), (7, 12), (15, 30), (25, 55), (45, 51)]
    result   = list(merge(segments))
    result   = [(1, 12), (15, 55)]
    """

    # store number of segments
    nsegs = len(segments)

    # segments need to be pre-sorted by segment starts but they may overlap
    assert all([segments[i][0] <= segments[i+1][0] for i in xrange(nsegs-1)])

    # initialize the current start and stop with the first segment
    current_start, current_stop = segments[0]
    current_value = values[0]

    # process the remaining segments merging where necessary
    for idx in xrange(nsegs-1):
        start, stop = segments[idx+1]
        value = values[idx+1]

        # new segment disjoint with previous: yield previous segment
        if start > current_stop:
            yield current_start, current_stop, current_value
            current_start, current_stop, current_value = start, stop, value

        # new segment overlaps previous segment: merge previous+new
        else:
            # get length and value of each segment
            lcur = float(current_stop - current_start)
            lnew = float(stop - start)
            ltot = lcur + lnew

            # set current value to length-weighted mean val to nearest int
            wavg = (lcur*current_value + lnew*value) / ltot
            current_value = int(0.5 + wavg)
            current_stop = max(stop, current_stop)

    # yield the final segment and stop iteration
    yield current_start, current_stop, current_value
