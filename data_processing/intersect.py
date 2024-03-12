__author__ = "davidmurphy"


"""
A function used to compare 2 (potentially more if script is modified) sets of segments and return a new set of
segments of overlapping regions of the input segments
"""


def intersect(segments1, segments2):
    """
    A generator that yields all intersections of segments1 and segments2.
    ** NOTE: In principle this function can be modified to take more than 2
    sets of segments **
    :param segments1: segment list 1 in format [(start1, end1)..(startN, endN)]
    :param segments2: segment list 2 in format [(start1, end1)..(startM, endM)]
    Example:
    --------
    segments1 = [(1, 10), (25, 55)]
    segments2 = [(7, 12), (15, 30), (45, 51)]
    result    = [(7, 10), (25, 30), (45, 51)]
    """
    # segments must be sorted and consecutive segments cannot overlap.
    # lengths must be non-negative
    for s in [segments1, segments2]:
        assert all([s[i][1] < s[i + 1][0] for i in range(len(s) - 1)])
        assert all([(start <= end) for (start, end) in s])

    # turn each list of segments into an iterator
    segments1, segments2 = map(iter, [segments1, segments2])

    # initialize each start and end point from the segment lists
    (s1, e1) = next(segments1)
    (s2, e2) = next(segments2)
    (start, end) = (0, 0)

    # loop through the segments and find intersections until all segments are
    # exhausted
    while True:

        try:
            # find the max starting point of the two segments
            start = max(s1, s2)

            # find the min end point and advance to next for the segment list
            # with the min end point
            if e1 < e2:
                end = e1
                (s1, e1) = next(segments1)
            else:
                end = e2
                (s2, e2) = next(segments2)

            # if the max start to min end has positive length, yield the
            # segment as an intersection
            if start <= end:
                yield (start, end)

        except StopIteration:
            # if the last intersection had positive length, yield
            if start <= end:
                yield (start, end)
            break
