from typing import Optional, Tuple, List

from chromosome import *
from pedigree import *



class AdmixtureFraction:
    """
        This class represents an individual using a list of admixture
        fractions (one for each chromosome pair).

        This representation is intended to display the results of an admixture
        history in a way that facilitates the computation of summary
        statistics. It is not used in simulations.
    """

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Initialization
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def __init__(self, adm_frac: List[float]):
        if not (abs( sum(adm_frac) - 1.0 ) < 10e-9):
            raise ValueError("The admixture fractions of an individual must sum to 1.")
        
        self.adm_frac = adm_frac
    
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Observation methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    
    def get_admixture_fraction(self, p: int) -> float:
        return self.adm_frac[p-1]



class Individual:
    """
        This class represents an individual as a set of chromosome pairs.
        
        This representation is intended to display the results of an admixture
        history in a way that facilitates the computation of summary
        statistics. It is not used in simulations.
    """

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Initialization
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def __init__(self, chr1s: List[Chromosome], chr2s: List[Chromosome]):
        """
            Parameters
            ----------
            chr1
                The list of chromosomes (from parent 1).
            chr2
                The list of chromosomes (from parent 2).
        """
        chr_nb = len(chr1s)
        if not (chr_nb == len(chr2s)):
            raise ValueError("The numbers of maternal and paternal chromosomes must be the same.")
        
        self.chrs = []
        for chr1, chr2 in zip(chr1s, chr2s):
            if not (chr1.get_length() == chr2.get_length()):
                raise ValueError("The lengths of two homologous chromosomes must be the sames.")
            self.chrs.append((chr1, chr2))
    
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Observation methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def get_chromosome(self, p: int, c: int) -> Chromosome:
        """
            This method returns the chromosome `c` in the pair `p`.
        """
        return self.chrs[p][c-1]
    
    def get_admixture_fraction(self, p: int) -> float:
        """
            This method returns the admixture fraction with respect to source
            population `p`.
        """
        H, L = 0, 0
        for chr1, chr2 in self.chrs:
            length = chr1.get_length()
            L += 2 * length
            h = chr1.get_admixture_fraction(p) + chr2.get_admixture_fraction(p)
            H += length * h
        return H / L



class AdmixtureFractionHistory:
    """
        This class represents an admixture history with every admixture
        fraction that existed and all kinship links between them.
        
        This representation is intended to compute the genetic diversity of
        finitely many infinitely long chromosomes (equivalently, infinite
        chromosomal crossover intensity).
    """

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Initialization
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def __init__(self, contrib: List[List[float]], pop_sizes: List[int], endogamy: Optional[float] = 0.0):
        """
            Parameters
            ----------
            contrib
                The contributions of each source populations at each
                generation.
            pop_sizes
                The sizes of the admixed population at each generation.
            endogamy (optional)
                A coefficient controlling the closeness of relatives in mating
                patterns: large and positive for endogamy; small or negative
                for exogamy (can be infinite, default 0.0).
        """
        source_nb = len(contrib)
        gen_nb = len(pop_sizes)

        self.pedigree = Pedigree(contrib, pop_sizes, endogamy)
        self.adm_frac: List[List[List[float]]] = [[[None for _ in range(source_nb)] for _ in range(pop_sizes[t])] for t in range(gen_nb)]
    
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Auxiliary methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def _compute_adm_frac(self, t: int, i: int):
        """
            This method computes the admixture fractions of the individual `i`
            at generation `t`.
        """
        if self.adm_frac[t][i][0] is None:
            source_nb = len(self.adm_frac[t][i])
            parent1, parent2 = self.pedigree.get_kinship_of(t, i)

            if parent1 < 0:
                adm_frac1 = [0.0 for _ in range(source_nb)]
                adm_frac1[-parent1-1] = 1.0
            else:
                self._compute_adm_frac(t-1, parent1)
                adm_frac1 = self.adm_frac[t-1][parent1]
            
            if parent2 < 0:
                adm_frac2 = [0.0 for _ in range(source_nb)]
                adm_frac2[-parent2-1] = 1.0
            else:
                self._compute_adm_frac(t-1, parent2)
                adm_frac2 = self.adm_frac[t-1][parent2]
            
            for p in range(source_nb):
                self.adm_frac[t][i][p] = (adm_frac1[p] + adm_frac2[p]) / 2
            
    
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Observation methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def sample(self, size: Optional[int] = 1, t: Optional[int] = -1) -> List[AdmixtureFraction]:
        """
            This method computes and returns admixture fractions of `size`
            different individuals in generation `t`.
        """
        if size > self.pedigree.pop_sizes[t]:
            raise ValueError("The sample size must be smaller than the sampled population size.")
        
        indivs = []
        for i in range(size):
            self._compute_adm_frac(t, i)
            indivs.append(AdmixtureFraction(self.adm_frac[t][i]))
        
        return indivs



class AdmixtureHistory:
    """
        This class represents an admixture history with every admixed
        chromosome that existed and all kinship links between them.
    """

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Initialization
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def __init__(self, contrib: List[List[float]], pop_sizes: List[int], chr_lengths: List[float], endogamy: Optional[float] = 0.0, intensity: Optional[float] = 1.0, fitness: Optional[float] = 0.0, selected_locus: Optional[float] = 0.0):
        """
            Parameters
            ----------
            contrib
                The contributions of each source populations at each
                generation.
            pop_sizes
                The sizes of the admixed population at each generation.
            chr_lengths
                The list of each chromosome pair lengths.
            endogamy (optional)
                A coefficient controlling the closeness of relatives in mating
                patterns: large and positive for endogamy; small or negative
                for exogamy (can be infinite, default 0.0).
            intensity (optional)
                The intensity for the Poisson point process modeling
                recombination (can be infinite, default 1.0).
            fitness (optional)
                The additional odd to one of mating for selected individuals
                (must be at least -1, can be infinite, default 0.0).
            selected_locus (optional)
                The relative position on the first chromosome pair on which
                alleles originating from source population 1 are selected
                (between 0.0 and 1.0, default 0.0).
        """
        gen_nb = len(pop_sizes)
        self.chr_nb = len(chr_lengths)

        self.pedigree = GeneticPedigree(contrib, pop_sizes, chr_lengths, endogamy, intensity, fitness, selected_locus)
        self.chromosomes: List[List[List[Chromosome]]] = [[[Chromosome(length) for length in chr_lengths] for _ in range(2 * pop_sizes[t])] for t in range(gen_nb)]
    
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Auxiliary methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def _compute_chromosome(self, t: int, i: int, p: int, c: int, a: Optional[float] = 0.0, b: Optional[float] = None):
        """
            This method computes the local ancestries between loci `a` and `b`
            on the chromosome `c` in pair `p` of the individual `i` at
            generation `t`.
        """
        unknown_regions = self.chromosomes[t][2*i+c][p].get_unknown_origins_between(a, b)
        parent = self.pedigree.get_kinship_of(t, i)[c]
        inter_recomb, intra_recomb = self.pedigree.get_recombination_event_of(t, i, p, c)

        # Cut regions of unknown origin using recombination event and associate them to a chromosome in the parental pair.
        parental_regions = []
        j = 0
        for a, b in unknown_regions:
            while j < len(intra_recomb) and intra_recomb[j] < a:
                j += 1
            seg_a = a
            while j < len(intra_recomb) and intra_recomb[j] < b:
                parent_c = (j + inter_recomb) % 2
                parental_regions.append((seg_a, intra_recomb[j], parent_c))
                seg_a = intra_recomb[j]
                j += 1
            parent_c = (j + inter_recomb) % 2
            parental_regions.append((seg_a, b, parent_c))

        # Compute the origin based coloring of each parental region and retrieve it.
        for parent_a, parent_b, parent_c in parental_regions:
            if parent < 0:
                self.chromosomes[t][2*i+c][p].add_ibd_segment(parent_a, parent_b, -parent)
            else:
                self._compute_chromosome(t-1, parent, p, parent_c, parent_a, parent_b)
                ibd_segments = self.chromosomes[t-1][2*parent+parent_c][p].get_origins_between(parent_a, parent_b)
                self.chromosomes[t][2*i+c][p].add_ibd_segments(ibd_segments)
    
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Observation methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def sample(self, size: Optional[int] = 1, t: Optional[int] = -1) -> List[Individual]:
        """
            This method computes and returns painted chromosomes of `size`
            different individuals in generation `t`.
        """
        if size > self.pedigree.pop_sizes[t]:
            raise ValueError("The sample size must be smaller than the sampled population size.")
        
        indivs = []
        for i in range(size):
            chr1s, chr2s = [], []
            for p in range(self.chr_nb):
                self._compute_chromosome(t, i, p, 0)
                chr1s.append(self.chromosomes[t][2*i][p])
                self._compute_chromosome(t, i, p, 1)
                chr2s.append(self.chromosomes[t][2*i+1][p])
            indivs.append(Individual(chr1s, chr2s))
        
        return indivs
