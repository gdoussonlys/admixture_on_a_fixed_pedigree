from admixture_history import *
from numpy import random

from multiprocessing import Pool
from functools import partial

def dirac(contrib, pop_sizes, endogamy, fitness, selected_locus):
    return contrib, pop_sizes, endogamy, fitness, selected_locus

class SimulationFromPrior:
    """
        This class is for managing simulation of multiple independent and
        identically distributed admixture histories.
    """

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Initialization
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def __init__(self, prior, chr_lengths, intensity=1):
        """
            Parameters
            ----------
            prior
                A function taking no arguments and returning contributions,
                population sizes, an endogamy coefficient, a fitness
                coefficient and a selected locus (when alleles originates from
                source population 1).
            chr_lengths
                The list of each chromosome pair lengths.
            intensity (optional)
                The intensity for the Poisson point process modeling
                recombination (can be infinite, default 1.0).
        """
        self.prior = prior

        self.intensity = intensity
        
        self.chr_nb = len(chr_lengths)
        no_recomb, infinite_recomb = True, False
        for chr_length in chr_lengths:
            no_recomb = no_recomb and (chr_length == 0)
            infinite_recomb = infinite_recomb or (chr_length == float("inf"))
        if no_recomb:
            self.chr_lengths = []
            self.intensity = 0
            for _ in range(self.chr_nb):
                self.chr_lengths.append(1)
        elif infinite_recomb:
            self.chr_lengths = []
            self.intensity = float("inf")
            for chr_length in chr_lengths:
                if chr_lengths == float("inf"):
                    self.chr_lengths.append(1)
        else:
            self.chr_lengths = chr_lengths

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Auxiliary methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def _sample_one_history(self, rep, size=1, gen=-1):
        """
            This method simulates an admixture history (with seed `rep`) and
            sample `size` different individuals in generation `gen`.
        """
        random.seed(rep)
        params = self.prior()
        contrib, pop_sizes, endogamy, fitness, selected_locus = params
        
        if self.intensity == float("inf"):
            if fitness == 0.0:
                history = AdmixtureFractionHistory(contrib, pop_sizes, endogamy)
            else:
                raise ValueError("The admixed population cannot be under selective pressure if recombination intensity is infinite.")
        else:
            history = AdmixtureHistory(contrib, pop_sizes, self.chr_lengths, endogamy, self.intensity, fitness, selected_locus)
        
        sample = history.sample(size, gen)
        return sample, params
    
    def _measure_one_sample(self, rep, sample, params, param=True, adm_frac=False, ibd_seg=False):
        """
            This method takes a `sample` of individuals and the corresponding
            seed `rep` and parameters `params`, compute statistics (`adm_frac`
            or `ibd_seg`), and returns them in long data format.
        """
        measurements = []

        contrib, pop_sizes, endogamy, fitness, selected_locus = params
        source_nb = len(contrib)

        rep_measurement = [rep]

        if adm_frac or ibd_seg:
            for i in range(len(sample)):
                indiv_measurement = [i]
                individual = sample[i]

                param_measurement = []
                if param:
                    param_measurement.append(endogamy)
                    param_measurement.append(fitness)
                    param_measurement.append(selected_locus)

                if adm_frac:
                    for p in range(1, source_nb + 1):
                        source_measurement = [p]
                        source_measurement.append(individual.get_admixture_fraction(p))
                            
                        measurement = rep_measurement + param_measurement + indiv_measurement + source_measurement
                        measurements.append(measurement)

                if ibd_seg:
                    for p in range(self.chr_nb):
                        pair_measurement = [p + 1]
                        for c in [1, 2]:
                            chrm_measurement = [c]
                            chrm = individual.get_chromosome(p, c)

                            ibd_segments = chrm.get_origins_between()
                            s = 0
                            for seg_a, seg_b, origin in ibd_segments:
                                ibd_seg_measurement = [s, origin, seg_b - seg_a]
                                measurement = rep_measurement + param_measurement + indiv_measurement + pair_measurement + chrm_measurement + ibd_seg_measurement
                                measurements.append(measurement)
                                s += 1

        else:
            measurement = rep_measurement
            measurements.append(measurement)

        return measurements
    
    def _measure_one_history(self, rep, size=1, gen=-1, param=True, adm_frac=False, ibd_seg=False):
        """
            This method combines `_sample_one_history` and
            `_measure_one_sample`.
        """
        sample, params = self._sample_one_history(rep, size, gen)
        measurements = self._measure_one_sample(rep, sample, params, param, adm_frac, ibd_seg)
        return measurements
    
    def label_data(self, param=True, adm_frac=False, ibd_seg=False):
        """
        This method produces labels for the long format dataset produced by
        `_measure_one_history`.
        """
        labels = []
        indiv_scale = adm_frac or ibd_seg

        labels.append("replicate")

        if param:
            labels.append("endogamy")
            labels.append("fitness")
            labels.append("selected_locus")

        if adm_frac:
            labels.append("individual")
            labels.append("source")
            labels.append("adm_frac")
    
        if ibd_seg:
            labels.append("individual")
            labels.append("pair")
            labels.append("chromosome")
            labels.append("ibd_segment")
            labels.append("source")
            labels.append("length")

        return labels

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Observation methods
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def measure(self, rep_nb=1, size=1, gen=-1, param=True, adm_frac=False, ibd_seg=False):
        """
            This method simulates `rep_nb` different replicates, samples
            individuals within them, compute the requested statistics and
            returns them in long format along with corresponding labels.
        """
        if adm_frac and ibd_seg:
            raise ValueError("Measuring both admixture fractions and IBD segments is ambiguous.")
        
        measure_one_history = partial(self._measure_one_history, size=size, gen=gen, param=param, adm_frac=adm_frac, ibd_seg=ibd_seg)
        with Pool(processes=64) as pool:
            measurements = pool.map(measure_one_history, range(rep_nb))
        
        data = []
        for measurements_rep in measurements:
            for measurement in measurements_rep:
                data.append(measurement)

        labels = self.label_data(param, adm_frac, ibd_seg)
        
        return data, labels

class Simulation(SimulationFromPrior):
    """
        This class is for managing simulation of multiple independent and
        identically parametrized admixture histories.
    """

    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    #    Initialization
    # −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

    def __init__(self, contrib, pop_sizes, chr_lengths, endogamy = 0.0, intensity = 1.0, fitness = 0.0, selected_locus = 0.0):
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
        self.source_nb = len(contrib)
        self.contrib = contrib

        self.pop_sizes = pop_sizes
        self.endogamy = endogamy

        self.fitness = fitness
        self.selected_locus = selected_locus
        
        prior = partial(dirac, contrib, pop_sizes, endogamy, fitness, selected_locus)
        super().__init__(prior, chr_lengths, intensity=intensity)
