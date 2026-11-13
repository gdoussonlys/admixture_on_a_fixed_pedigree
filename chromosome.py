from typing import Optional, Tuple, List

from bisect import bisect_right, insort



class Chromosome:
    """
        This class represents a chromosome (or a chromosome segment) using an
        interval along with origin-based coloring (known as chromosome painting
        or local ancestries).
    """

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Initialization
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def __init__(self, length: float):
        """
            Parameters
            ----------
            length
                The chromosome length.
        """
        if length <= 0:
            raise ValueError("The length of a chromosome must be positive.")
        self.length = float(length)
        self.ibd_segments: List[IBD_Segment] = []
        self._starts: List[float] = []
    
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Auxiliary methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    
    def _check_bounds(self, a: float, b: float):
        if not (a < b):
            raise ValueError("The bounds of a non-empty interval must be increasing.")
        if not (0 <= a and b <= self.length):
            raise ValueError("The interval must be included in the chromosome.")

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Observation methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def get_length(self) -> float:
        return self.length

    def get_origin_at(self, x: float) -> int:
        """
            This method returns the origin of the allele at locus `x`.
        """
        if not (0 <= x < self.length):
            raise ValueError("The locus must belong to the chromosome.")

        # Find the first IBD segment next to the locus.
        i = bisect_right(self._starts, x)
        if i == 0:
            return 0
        
        # Test whether the previous IBD segment covers the locus.
        a, b, origin = self.ibd_segments[i - 1]
        if x < b:
            return origin
        return 0

    def get_origins_between(self, a: Optional[float] = 0.0, b: Optional[float] = None, intersect: Optional[bool] = True) -> List[Tuple[float, float, int]]:
        """
            This method returns the origin of the alleles between loci `a` and
            `b`.
        """
        if b is None:
            b = self.length
        self._check_bounds(a, b)
        ibd_segments = []
        i = bisect_right(self._starts, a)

        # Check whether the previous segment covers a.
        if i > 0:
            left = self.ibd_segments[i-1]
            if left[1] > a:
                if intersect:
                    seg_a = a
                    seg_b = min(left[1], b)
                    ibd_segments.append((seg_a, seg_b, left[2]))
                else:
                    ibd_segments.append(left)

        # Iterate over the next segments until exceeding b.
        j = i
        while j < len(self.ibd_segments) and self.ibd_segments[j][0] < b:
            right = self.ibd_segments[j]
            if intersect:
                seg_a = max(right[0], a)
                seg_b = min(right[1], b)
                ibd_segments.append((seg_a, seg_b, right[2]))
            else:
                ibd_segments.append(right)
            j += 1

        return ibd_segments
    
    def get_unknown_origins_between(self, a: Optional[float] = 0.0, b: Optional[float] = None) -> List[Tuple[float, float]]:
        """
            This method returns interval bounds between which local ancestries
            are unknown.
        """
        if b is None:
            b = self.length
        self._check_bounds(a, b)
        unknown_regions = []

        i = bisect_right(self._starts, a)
        seg_a = a

        # Check whether the previous segment covers a.
        if i > 0 and self.ibd_segments[i-1][1] > a:
            seg_a = min(self.ibd_segments[i-1][1], b)

        # Iterate over the next segments until exceeding b.
        j = i
        while j < len(self.ibd_segments) and self.ibd_segments[j][0] < b:
            right = self.ibd_segments[j]
            if right[0] > seg_a:
                seg_b = min(right[0], b)
                unknown_regions.append((seg_a, seg_b))
            if right[1] > seg_a:
                seg_a = min(right[1], b)
            j += 1
        
        if seg_a < b:
            unknown_regions.append((seg_a, b))

        return unknown_regions
    
    def get_admixture_fraction(self, p: int) -> float:
        """
            This method computes the admixture fraction with respect to source
            population `p` (ignoring regions of unknown ancestry).
        """
        if p <= 0:
            raise ValueError("A source population must be indexed by a positive integer.")
        h = 0.0
        for a, b, origin in self.ibd_segments:
            if origin == p:
                h += (b - a)
        return h / self.length
    
    def to_text(self, char_nb: Optional[int] = 100) -> str:
        """
            This method returns a colored string representation of the
            chromosome painting.
        """
        text = ""
        colors = ["\033[30m", "\033[31m", "\033[34m", "\033[32m", "\033[35m", "\033[36m", "\033[33m", "\033[37m"]
        morgan_per_char = self.length / char_nb
        locus = 0
        for _ in range(char_nb):
            p = self.get_origin_at(locus)
            text += colors[p] + '█'
            locus += morgan_per_char
        text += "\033[30m"
        return text

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Editing methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def add_ibd_segment(self, a: float, b: float, origin: int):
        """
            This method adds an IBD segment, ranging from `a` to `b` along with
            its `origin`, merging contiguous or overlapping segments as needed.
        """
        if origin <= 0:
            raise ValueError("The origin of an IBD segment must be indexed by a positive integer.")
        a = float(a); b = float(b)
        self._check_bounds(a, b)

        new_a, new_b = a, b
        i = bisect_right(self._starts, a)

        to_delete = []

        # Check the segment to the left that may overlap.
        if i > 0:
            left = self.ibd_segments[i-1]
            if left[1] > a and left[2] != origin:
                    raise ValueError(f"The IBD segments [{a}, {b}) and [{left[0]}, {left[1]}) overlap, but their origins differ.")
            if left[1] > a or (left[1] == a and left[2] == origin):
                new_a = min(new_a, left[0])
                new_b = max(new_b, left[1])
                to_delete.append(i - 1)

        # Check the segments to the right that may overlap.
        while i < len(self.ibd_segments):
            right = self.ibd_segments[i]
            if right[0] > b or (right[0] == b and right[2] != origin):
                break
            if right[0] < b and right[2] != origin:
                raise ValueError(f"The IBD segments [{a}, {b}) and [{right[0]}, {right[1]}) overlap, but their origins differ.")
            new_a = min(new_a, right[0])
            new_b = max(new_b, right[1])
            to_delete.append(i)
            i += 1

        # Remove all merged segments
        for idx in reversed(to_delete):
            del self.ibd_segments[idx]
            del self._starts[idx]

        # Add the complete segment
        insort(self.ibd_segments, (new_a, new_b, origin))
        insort(self._starts, new_a)
    


    def add_ibd_segments(self, segments: List[Tuple[float, float, int]]):
        """
            This method adds multiple IBD segments using `add_ibd_segment`.
        """
        segments = sorted(segments, key=lambda s: s[0])
        for a, b, origin in sorted(segments):
            self.add_ibd_segment(a, b, origin)
