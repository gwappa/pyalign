#
# MIT License
#
# Copyright (c) 2019 Keisuke Sehara
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

from math import inf as _INFTY
from collections import namedtuple as _namedtuple

VERSION_STR = '1.0.0a1'

DEBUG = True

def debug(msg, end='\n'):
    if DEBUG is True:
        print(f"... {msg}", end=end, flush=True)

Aligned = _namedtuple('Aligned', ('first', 'second'))

class AlignmentResult(_namedtuple('_AlignmentResult', ('aligned', 'loss'))):
    @staticmethod
    def from_pair(item1, item2, distance):
        try:
            return AlignmentResult([ Aligned(item1, item2) ], distance(item1, item2))
        except ValueError:
            return AlignmentResult([], _INFTY)

    def __gt__(self, other):
        if not isinstance(other, AlignmentResult):
            raise ValueError(f"expected AlignmentResult, got {other.__class__.__name__}")
        return self.loss > other.loss

    def __ge__(self, other):
        if not isinstance(other, AlignmentResult):
            raise ValueError(f"expected AlignmentResult, got {other.__class__.__name__}")
        return self.loss >= other.loss

    def __lt__(self, other):
        if not isinstance(other, AlignmentResult):
            raise ValueError(f"expected AlignmentResult, got {other.__class__.__name__}")
        return self.loss < other.loss

    def __le__(self, other):
        if not isinstance(other, AlignmentResult):
            raise ValueError(f"expected AlignmentResult, got {other.__class__.__name__}")
        return self.loss <= other.loss

    def __add__(self, other):
        if not isinstance(other, AlignmentResult):
            raise ValueError(f"expected AlignmentResult, got {other.__class__.__name__}")
        return self.__class__(self.aligned + other.aligned, self.loss + other.loss)

def difference(item1, item2):
    return abs(item1 - item2)

def SmithWaterman(first_seq, second_seq, distance=None):
    """perform alignment between `first_seq` and `second_seq`
    based on the Smith-Waterman method.

    `distance` must be the function that takes the form of `distance(first_item, second_item)`
    and returns the distance (loss function) between the two. If either is a skip, None will
    be assigned.

    Returns the AlignmentResult named tuple.
    """
    first_size  = len(first_seq)
    second_size = len(second_seq)
    if first_size > 900 or second_size > 900:
        raise ValueError("size of the sequence may be too long; consider splitting in pieces")

    if distance is None:
        distance = difference

    # prepare memo
    memo = {}

    def _compute_alignment(offset1, offset2):
        debug(f"--> {offset1}, {offset2}")
        if (offset1, offset2) not in memo.keys():
            if offset1 == first_size:
                if offset2 == second_size:
                    memo[offset1, offset2] = AlignmentResult([], 0)
                else:
                    memo[offset1, offset2] = AlignmentResult.from_pair(None,
                                                                second_seq[offset2],
                                                                distance) \
                                                + _compute_alignment(offset1, offset2+1)
            elif len(_seq2) == 0:
                memo[offset1, offset2] = AlignmentResult.from_pair(first_seq[offset1],
                                                                None,
                                                                distance) \
                                            + _compute_alignment(offset1+1, offset2, base_loss=distance(first_seq[offset1], None))
            else:
                # compute recursively
                item1      = first_seq[offset1]
                item2      = second_seq[offset2]

                match      = AlignmentResult.from_pair(item1, item2, distance) \
                                + _compute_alignment(offset1+1, offset2+1)
                skip_item1 = AlignmentResult.from_pair(None, item2, distance) \
                                + _compute_alignment(offset1,   offset2+1)
                skip_item2 = AlignmentResult.from_pair(item1, None, distance) \
                                + _compute_alignment(offset1+1, offset2)
                memo[offset1, offset2] = min(match, skip_item1, skip_item2)

        # must have already computed
        debug(f"<-- {offset1}, {offset2} = {memo[offset1, offset2].loss}")
        return memo[offset1, offset2]

    return _compute_alignment(0, 0)
