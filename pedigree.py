from typing import Optional, Tuple, List

from numpy.random import randint, random, exponential, choice
from math import exp



class Pedigree:
    """
        This class represents a pedigree of an admixed population using a
        `kinship` matrix associating each admixed individual to two distinct
        parents in the previous generation.
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

        # Check whether at least one source population exists.
        if source_nb == 0:
            raise ValueError("There must be at least one source population.")
        
        # Check whether the number of generations matches across parameters.
        for p in range(source_nb):
            if len(contrib[p]) != gen_nb:
                raise ValueError("The number of generations must be the same throughout all parameters.")

        # Check whether the contributions are non-negative and sum to 1.
        for t in range(gen_nb):
            h = 1
            for p in range(source_nb):
                s_p = contrib[p][t]
                h -= s_p
                if s_p < 0:
                    raise ValueError("The contributions of the source populations must be non-negative.")
            if t == 0 and h > 0:
                raise ValueError("The contributions of the source populations must sum to 1 at the foundation.")
            if h < 0:
                raise ValueError("The contributions of the source populations must not sum to more than 1.")

        self.kinship: List[List[Tuple[int, int]]] = [[None for _ in range(2 * pop_sizes[t])] for t in range(gen_nb)]
        self.contrib: List[List[float]] = contrib
        self.pop_sizes: List[int] = pop_sizes

        self.endogamy: float = endogamy
        if endogamy != 0.0:
            self.tlca: List[List[List[int]]] = [[[None for _ in range(pop_sizes[t])] for _ in range(pop_sizes[t])] for t in range(gen_nb)]
    
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Auxiliary methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    @staticmethod
    def normalize(weights: List[float]) -> List[float]:
        """
            This method converts a list of `weights` into a distribution
            through normalization.
        """
        if max(weights) == float("inf"):
            for k in range(len(weights)):
                if weights[k] == float("inf"):
                    weights[k] = 1
                else:
                    weights[k] = 0
        Z = sum(weights)
        for k in range(len(weights)):
            weights[k] /= Z
        return weights
    
    @staticmethod
    def boltzmann_gibbs_distribution(energies: List[float], beta: Optional[float] = 0.0) -> List[float]:
        """
            This method converts a list of `energies` into a Boltzmann-Gibbs
            distribution.
        """
        weights = []
        if abs(beta) == float("inf"):
            if beta == float("inf"):
                extremum = min(energies)
            else:
                extremum = max(energies)
            for E in energies:
                if E == extremum:
                    weights.append(1.0)
                else:
                    weights.append(0.0)
        else:
            for E in energies:
                weights.append(exp(- beta * E))
        Z = sum(weights)
        for k in range(len(weights)):
            weights[k] /= Z
        return weights

    def _draw_population(self, t: int) -> int:
        """
            This method randomly samples a population according to the
            contributions at generation `t`.
        """
        source_nb = len(self.contrib)

        r = random()
        for p in range(source_nb):
            r -= self.contrib[p][t]
            if r < 0:
                return p + 1
        return 0
    
    def _get_tlca_of(self, t: int, i1: int, i2: int) -> int:
        """
            This method retrieves the time to the last common ancestor (TLCA)
            of two individuals `i1` and `i2` of a same generation `t`.
        """
        if self.endogamy == 0.0:
            raise ValueError("The admixed population must be endo/exogamous to be interested in TLCAs.")
        
        if i1 < 0 or i2 < 0:
            return float("inf")
        
        if self.tlca[t][i1][i2] is None:
            if i1 == i2:
                tlca = 0
            else:
                p11, p12 = self.get_kinship_of(t, i1)
                p21, p22 = self.get_kinship_of(t, i2)
                tlca11 = self._get_tlca_of(t-1, p11, p21)
                tlca12 = self._get_tlca_of(t-1, p11, p22)
                tlca21 = self._get_tlca_of(t-1, p12, p21)
                tlca22 = self._get_tlca_of(t-1, p12, p22)
                tlca = min(tlca11, tlca12, tlca21, tlca22) + 1
            self.tlca[t][i1][i2] = tlca
            self.tlca[t][i2][i1] = tlca
        
        return self.tlca[t][i1][i2]
    
    def _draw_relative_of(self, t: int, i: int) -> int:
        """
            This method randomly samples an individual in generation `t`
            according to the Boltzmann-Gibbs distribution, with TLCA regarding
            individual `i` as energies.
        """
        js = [j for j in range(i)] + [j for j in range(i+1, self.pop_sizes[t])]
        tlcas = [min(self._get_tlca_of(t, i, j), 5) for j in js]
        distrib = self.boltzmann_gibbs_distribution(tlcas, self.endogamy)
        return choice(js, p=distrib)

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Observation methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def get_kinship_of(self, t: int, i: int) -> Tuple[int, int]:
        """
            This method returns the two parents of the individual `i` at
            generation `t`. Source populations are represented by negative
            integers.
        """
        possible_selfing = False
        
        # Draw a source population or an admixed parent for the maternal chromosome.
        if self.kinship[t][2*i] is None:
            p = self._draw_population(t)
            if p == 0:
                possible_selfing = True
                self.kinship[t][2*i] = randint(0, self.pop_sizes[t-1])
            else:
                self.kinship[t][2*i] = -p

        # Draw a source population or an admixed parent for the paternal chromosome.
        if self.kinship[t][2*i+1] is None:
            p = self._draw_population(t)
            if p == 0:
                possible_selfing = True
                self.kinship[t][2*i+1] = randint(0, self.pop_sizes[t-1])
            else:
                self.kinship[t][2*i+1] = -p
        
        # Change the admixed parent as long as they are the same (to avoid selfing).
        if possible_selfing:
            if self.endogamy == 0:
                js = [j for j in range(self.kinship[t][2*i])] + [j for j in range(self.kinship[t][2*i]+1, self.pop_sizes[t-1])]
                distrib = self.normalize([1 for _ in js])
                self.kinship[t][2*i+1] = choice(js, p=distrib)
            else:
                self.kinship[t][2*i+1] = self._draw_relative_of(t-1, self.kinship[t][2*i])
        
        return self.kinship[t][2*i], self.kinship[t][2*i+1]



class GeneticPedigree(Pedigree):
    """
        This class extends the `Pedigree` class with recombination events,
        represented by a Poisson point process for crossover points and by
        a Bernoulli variable (with probability 1/2) for segregation.
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
        # Check whether endogamy and fitness are not both simulated.
        if not (endogamy == 0.0 or fitness == 0.0):
            raise ValueError("The admixed population cannot be both endo/exogamous and under selective pressure.")
        
        super().__init__(contrib, pop_sizes, endogamy)

        gen_nb = len(pop_sizes)
        chr_nb = len(chr_lengths)
        self.recombination: List[List[List[Tuple[int, List[float]]]]] = [[[None for _ in range(chr_nb)] for _ in range(2 * pop_sizes[t])] for t in range(gen_nb)]
        self.chr_lengths: List[float] = chr_lengths
        self.intensity: float = intensity

        self.fitness: float = fitness
        self.selected_locus: float = selected_locus
        if fitness != 0.0:
            self.mutant: List[List[bool]] = [[None for _ in range(2 * pop_sizes[t])] for t in range(gen_nb)]

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Auxiliary methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    @staticmethod
    def _poisson_point_process(length: float, intensity: Optional[float] = 1.0) -> List[float]:
        """
            This method randomly samples a Poisson point process on a finite
            interval ranging from 0 to `length`.
        """
        points = []
        if intensity > 0:
            point = exponential(1 / intensity)
            while point < length:
                points.append(point)
                point += exponential(1 / intensity)
        return points
    
    def _recombinate(self, length: float) -> Tuple[int, List[float]]:
        """
            This method draws crossover points using `_poisson_point_process`
            and segregation as a Bernoulli variable with probability 1/2.
        """
        segregation = randint(0, 2)
        crossovers = self._poisson_point_process(length, self.intensity)
        return segregation, crossovers
    
    def _get_type_of(self, t: int, i: int, c: int) -> bool:
        """
            This method checks whether the individual `i` at generation `t` is
            a mutant or not (i.e, carries at least one selected allele).
        """
        if self.fitness == 0.0:
            raise ValueError("The admixed population must be under selective pressure to be interested in wild/mutant types.")

        if self.mutant[t][2*i+c] is None:
            parent = self.get_kinship_of(t, i)[c]
            inter_recomb, intra_recomb = self.get_recombination_event_of(t, i, 0, c)
            parent_c = 0
            while parent_c < len(intra_recomb) and intra_recomb[parent_c] < self.selected_locus * self.chr_lengths[0]:
                parent_c += 1
            parent_c = (parent_c + inter_recomb) % 2
            
            if parent < 0:
                self.mutant[t][2*i+c] = (parent == -1)
            else:
                self.mutant[t][2*i+c] = self._get_type_of(t-1, parent, parent_c)
        
        return self.mutant[t][2*i+c]
    
    def _get_fitness_of(self, t: int, i: int) -> float:
        """
            This method returns the fitness of the individual `i` at generation
            `t`.
        """
        if self.fitness == 0.0:
            raise ValueError("The admixed population must be under selective pressure to be interested in fitnesses.")
        
        if self._get_type_of(t, i, 0) or self._get_type_of(t, i, 1):
            return self.fitness
        else:
            return 0.0
    
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Observation methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    
    def get_recombination_event_of(self, t: int, i: int, p: int, c: int) -> Tuple[int, List[float]]:
        """
            This method returns the recombination event that produced the
            chromosome `c` in pair `p` of the individual `i` at generation `t`.
        """
        if self.recombination[t][2*i+c][p] is None:
            self.recombination[t][2*i+c][p] = self._recombinate(self.chr_lengths[p])
        return self.recombination[t][2*i+c][p]
    
    def get_kinship_of(self, t: int, i: int) -> Tuple[int, int]:
        """
            This method returns the two parents of the individual `i` at
            generation `t`. Source populations are represented by negative
            integers.
        """
        if self.fitness == 0.0:
            return super().get_kinship_of(t, i)
        
        possible_selfing = False
        
        # Draw a source population or an admixed parent for the maternal chromosome.
        if self.kinship[t][2*i] is None:
            p = self._draw_population(t)
            if p == 0:
                possible_selfing = True
                js = [j for j in range(self.pop_sizes[t-1])]
                distrib = self.normalize([1.0 + self._get_fitness_of(t-1, j) for j in js])
                self.kinship[t][2*i] = choice(js, p=distrib)
            else:
                self.kinship[t][2*i] = -p

        # Draw a source population or an admixed parent for the paternal chromosome.
        if self.kinship[t][2*i+1] is None:
            p = self._draw_population(t)
            if p == 0:
                possible_selfing = True
                js = [j for j in range(self.pop_sizes[t-1])]
                distrib = self.normalize([1.0 + self._get_fitness_of(t-1, j) for j in js])
                self.kinship[t][2*i+1] = choice(js, p=distrib)
            else:
                self.kinship[t][2*i+1] = -p
        
        # Change the admixed father in case of possible selfing.
        if possible_selfing:
            js = [j for j in range(self.kinship[t][2*i])] + [j for j in range(self.kinship[t][2*i]+1, self.pop_sizes[t-1])]
            distrib = self.normalize([1.0 + self._get_fitness_of(t-1, j) for j in js])
            self.kinship[t][2*i+1] = choice(js, p=distrib)
        
        return self.kinship[t][2*i], self.kinship[t][2*i+1]
