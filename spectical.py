#!/usr/bin/env python3

import os
import argparse
import os.path
import timeit
import matplotlib.pyplot as plt
import numpy as np
import sys
import csv
import gzip
import warnings
import math
import json
import re
import multiprocessing

from scipy.optimize import OptimizeWarning
from scipy import interpolate
from pyteomics import mzml, auxiliary, mass
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from numpy import exp

warnings.simplefilter("ignore", OptimizeWarning)

class MSRunPeakFinder:

    def __init__(self, file, tolerance, rows, columns, make_pdf, find_snippets):
        # Declares all of the variables from parameters
        self.tolerance = tolerance
        self.tolerance_150 = 10
        self.plot_output_rows = rows
        self.plot_output_columns = columns
        self.file_name = file
        self.make_pdf = make_pdf
        self.find_snippets = find_snippets and make_pdf

        # Declares all of the variables
        self.by_count = np.zeros((4000001,), dtype=int)
        self.by_intensity = np.zeros((4000001,), dtype=float)
        self.by_strength = np.zeros((4000001,), dtype=int)
        self.all_peaks_intensities = [] # Keeps track of all intensities
        self.triggered_peaks = []
        self.observed_peaks = [] # Keeps track of all peaks over a certain intensities, identified and unidentified
        self.known_ions = {}
        self.scan_snippets = []
        self.minimum_triggers = 50
        self.top_peaks = 0.6
        self.ppm_delta = 15
        self.t0 = timeit.default_timer()
        self.stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }
        self.initial_tolerance = 20
        self.proton_mass = 1.007276
        self.electron_mass = 0.000548
        self.isotope_mass = 1.00335 # 1.00335
        self.crude_correction = 0
        self.after_400_calibration = 0
        self.has_correction_spline = False
        self.has_second_spline = False

        if self.file_name.endswith('.gz'):
            self.infile = gzip.open(self.file_name)
            self.peak_file = self.file_name[0:len(self.file_name) - 8]
        else:
            self.infile = open(self.file_name, 'rb')
            self.peak_file = self.file_name[0:len(self.file_name) - 5]

    def process_file(self):
        last_time = self.t0
        # find the intensity of each m/z value
        self.aggregate_spectra() # 
        print(str(timeit.default_timer() - last_time) + " seconds to aggregate all spectra")
        last_time = timeit.default_timer()
        self.determine_crude_correction() # 
        if self.crude_correction >= 15:
            self.increase_initial_tolerance()
            self.determine_crude_correction()
        print(str(timeit.default_timer() - last_time) + " seconds to determine crude calibration")
        last_time = timeit.default_timer()
        if self.make_pdf:
            self.plot_crude_calibrations() # 
        # create an array of all identifiable theoretical and their masses
        self.get_theoretical_ions() #
        print(str(timeit.default_timer() - last_time) + " seconds to find all theoretical ions")
        last_time = timeit.default_timer()
        # identify all the spectra with a certain intensity
        self.find_initial_triggers() #
        self.refine_triggered_peaks() # 
        if len(self.observed_peaks) < 200: # if there aren't enough peaks, lower the minimum trigger and redo everything
            self.lower_minimum_triggers()
            self.plot_crude_calibrations() # automatically clears the plots
            self.find_initial_triggers()
            self.refine_triggered_peaks()
        print(str(timeit.default_timer() - last_time) + " seconds to find and refine peaks")
        last_time = timeit.default_timer()
        self.identify_peaks() #
        print(str(timeit.default_timer() - last_time) + " seconds to identify peaks")
        last_time = timeit.default_timer()
        # save the identified peaks to a file
        self.delta_scatterplots() #
        self.find_first_spline()
        last_time = timeit.default_timer()
        self.identify_peaks()
        print(str(timeit.default_timer() - last_time) + " seconds to identify peaks")
        last_time = timeit.default_timer()
        self.find_peak_percentage()
        try:
            self.find_second_spline()    
        except:
            print("Error with creating the second spline")
        if self.make_pdf:
            self.plot_corrected_scatterplot()
            self.plot_all_corrections()
            print(str(timeit.default_timer() - last_time) + " seconds to plot all summary plots")
            last_time = timeit.default_timer()

        # following is code to write outputs - consider putting this in its own method
        if self.find_snippets:
            self.analyze_snippets() # skipping analyzing snippets for now, to add back, uncomment 158
            print(str(timeit.default_timer() - last_time) + " seconds to plot snippet plots")
            last_time = timeit.default_timer()
        self.write_output() #
        if self.make_pdf:
            self.plot_peaks_strength()
            print(str(timeit.default_timer() - last_time) + " seconds to plot individual peaks")
            last_time = timeit.default_timer()
        self.write_json()
        # print out data, including run time, number of peaks found, etc
        return self.show_stats()

    def aggregate_spectra(self):
        print("starting to aggregate spectra")

        with mzml.read(self.infile) as reader:
            for spectrum in reader:

                #### Update counters and print progress
                self.stats['counter'] += 1

                # Populates array with each mz based on count (number of times it appears), intensity
                # (summation of all intensities of the respective mz), and strength (summation of the strength
                # of the intensities, from a scale of 1-4)
                if spectrum['ms level'] == 1:
                    self.stats['ms1spectra'] += 1
                elif spectrum['ms level'] == 2 and 'm/z array' in spectrum:
                    self.stats['ms2spectra'] += 1
                    self.smallest_peak_intensity = sys.maxsize
                    if self.find_snippets:
                        lower_index = sys.maxsize
                        upper_index = 0
                    for index in range(len(spectrum['m/z array'])):
                        peak_mz = spectrum['m/z array'][index]
                        if self.find_snippets and peak_mz > 110:
                            lower_index = min(lower_index, index)
                            if peak_mz < 130:
                                upper_index = max(upper_index, index)
                        if peak_mz > 400:
                            break
                        else:
                            intensity = spectrum['intensity array'][index]
                            self.all_peaks_intensities.append(intensity)
                            self.by_count[int(10000 * peak_mz + 0.5)] += 1
                            self.by_intensity[int(10000 * peak_mz + 0.5)] += intensity
                            self.smallest_peak_intensity = min(self.smallest_peak_intensity, intensity)

                    if self.find_snippets and upper_index != 0 and upper_index > lower_index:
                        self.scan_snippets.append([self.stats['counter'], spectrum['m/z array'][lower_index:upper_index], spectrum['intensity array'][lower_index:upper_index]])

                    # Compare the intensity to the smallest intensity, returning the strength of the intensity
                    # on a scale from 1-4
                    for index in range(len(spectrum['m/z array'])):
                        peak_mz = spectrum['m/z array'][index]
                        if peak_mz > 400:
                            break
                        else:
                            intensity = spectrum['intensity array'][index]
                            self.by_strength[int(10000 * peak_mz + 0.5)] += get_strength(intensity, self.smallest_peak_intensity)

                # Updates terminal with the progress of reading peaks
                if self.stats['counter']/1000 == int(self.stats['counter']/1000):
                    print(f"  {self.stats['counter']}")

        # print(self.scan_snippets)

    def get_theoretical_ions(self):
        print("starting to get explanations")

        # Populates an array of arrays with a list of theoretical ions with the name, mass, and a boolean
        # value. The boolean value tracks if the theoretical ion has already been identified, which prevents
        # duplicate identifications
        ion_types = ['a', 'b', 'y']
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        amino_acid_modifications = {
            'C[Carbamidomethyl]': {'mz': 57.021464, 'amino acid': 'C'},
            'M[Oxidation]': {'mz': 15.994915, 'amino acid': 'M'},
            'M[Acetyl]': {'mz': 42.010565, 'amino acid': 'M'},
            'L[Acetyl]': {'mz': 42.010565, 'amino acid': 'L'}
        }

        # for amino_acid in amino_acids: # considers Acetylation
            # amino_acid_modifications[f'{amino_acid}[Acetyl]'] = {'mz': 42.010565, 'amino acid': amino_acid}

        for key in amino_acid_modifications:
            amino_acids.append(key)

        amino_acids.sort()

        # Considers losses for certain ion types
        aa_immonium_losses = {
                'G': [],
                'A': [],
                'S': [ '+CO+NH3' ], # added +CO+NH3
                'P': [],
                'V': [ '-CH2-NH3', '-NH3', '+CO-NH3-CH2' ],
                'T': [ '+CO-NH3', '+CO+NH3' ],
                'C': [],
                'L': [ '-C3H6', '-CH2' ],
                'I': [ '-C3H6', '-CH2' ],
                'N': [ '-NH3' ],
                'D': [ '-H2O' ],
                'Q': [ '-CO-NH3', '-NH3', '+CO'],
                'K': [ '+CO-NH3', '-NH3', '+CO', '-C2H4-NH3', '+CO+H2ON2', '-NH3', '-C4H7N', '+CO+CO-C2H3N3', '+CO+H2O', '+CO+H2O-NH3' ],
                'E': [],
                'M': [ '-C2H2-NH3' ],
                'H': [ '-CH2N', '+CO-NH2', '+CO-NH3', '+CO-NH', '+CO+H2O' ],
                'F': [ '-CH3N' ],
                'R': [ '-C3H6N2', '-CH5N3', '-CH6N2', '-C2H4N2', '-CH2N2', '-CH3N', '-NH3', '-C4H7N', '+H2O+H2O-N3H7', '+CO+H2O' ],
                'Y': [ '-CO-NH3', '-CH3N', '+CO+NH3' ], # added +CO+NH3
                'W': [ '+CO', '-C4H6N2', '-C2H4N', '-CH3N', '-CHN', '+CO-NH3', '-NH3' ],
            }
        modification_deltas = {}

        # Calculatse the additional deltas from possible explanations in aa_immonium_losses and
        # turns formula into numbers
        for modification in aa_immonium_losses:
            delta_mzs = {}
            for variant in aa_immonium_losses[modification]:
                delta_mz = 0
                variant_add = variant.split('-')[0].split('+')[1:]
                variant_subtract = variant.split('-')[1:]
                for variant_formula in variant_add:
                    delta_mz += mass.calculate_mass(formula=variant_formula)
                for variant_formula in variant_subtract:
                    delta_mz -= mass.calculate_mass(formula=variant_formula)
                delta_mzs[variant] = delta_mz
            modification_deltas[modification] = delta_mzs
        
        # Populates amino_acid_mass with theoretical ions and their mass
        for acid in amino_acids:
            for ion in ion_types:
                if len(acid) == 1:
                    acid_base = acid
                    base_acid_mz = mass.calculate_mass(sequence=acid, ion_type=ion, charge=1)
                else:
                    if acid[0] == '[':
                        acid_base = acid[8]
                    else:
                        acid_base = acid[0]
                    base_acid_mz = mass.calculate_mass(sequence=acid_base, ion_type=ion, charge=1) + amino_acid_modifications[acid]['mz']

                # Adds possible modifications, but only for a ions
                if ion == 'a':
                    # Add base a ion
                    self.known_ions[f"I{acid}"] = [base_acid_mz, False]
                    # Add base a ion's isotope
                    self.known_ions[f"I{acid}+i"] = [base_acid_mz + self.isotope_mass, False]
                    # Check if there are any possible modifications, if there are, add them too
                    for modification in modification_deltas[acid_base]:
                        acid_mz = base_acid_mz
                        acid_mz += modification_deltas[acid_base][modification]
                        self.known_ions[f"I{acid}{modification}"] = [acid_mz, False]
                        self.known_ions[f"I{acid}{modification}+i"] = [acid_mz + self.isotope_mass, False]
                # add b and y ions
                else:
                    self.known_ions[f"{ion}{{{acid}}}"] = [base_acid_mz, False]
                    self.known_ions[f"{ion}{{{acid}}}+i"] = [base_acid_mz + self.isotope_mass, False]

        # Double nested for loops to identify pairs of amino acids (i.e. A and A)
        # First amino acid
        for identifier_index_1 in range(len(amino_acids)):
            for ion in ion_types:
                pair_acid_1 = amino_acids[identifier_index_1]
                # Checks if it is Carbamidomethyl or Oxidation
                if len(pair_acid_1) == 1:
                    pair_1_acid_base = pair_acid_1
                    pair_mz_1 = mass.calculate_mass(sequence=pair_acid_1, ion_type=ion, charge=1)
                else:
                    if pair_acid_1[0] == '[':
                        pair_1_acid_base = pair_acid_1[8]
                    else:
                        pair_1_acid_base = pair_acid_1[0]
                    pair_mz_1 = amino_acid_modifications[pair_acid_1]['mz'] + mass.calculate_mass(sequence=pair_1_acid_base, ion_type=ion, charge=1)
                
                # Second amino acid
                for identifier_index_2 in range(len(amino_acids) - identifier_index_1):
                    identifier_index_2 += identifier_index_1
                    pair_acid_2 = amino_acids[identifier_index_2]
                    if len(pair_acid_2) == 1:
                        pair_2_acid_base = pair_acid_2
                        pair_mz_2 = mass.calculate_mass(sequence=pair_acid_2, ion_type='b', charge=0)
                    else:
                        if pair_acid_2[0] == '[':
                            pair_2_acid_base = pair_acid_2[8]
                        else:
                            pair_2_acid_base = pair_acid_2[0]
                        pair_mz_2 = amino_acid_modifications[pair_acid_2]['mz'] + mass.calculate_mass(sequence=pair_2_acid_base, ion_type='b', charge=0)
                    
                    # Add two amino acids together
                    pair_mz = pair_mz_1 + pair_mz_2
                    self.known_ions [f"{ion}{{{pair_acid_1}{pair_acid_2}}}"]  = [pair_mz, False]
                    # Considers isotope
                    self.known_ions [f"{ion}{{{pair_acid_1}{pair_acid_2}}}+i"] = [pair_mz + self.isotope_mass, False]

                    # considers doubly charged ions
                    # pair_mz = (pair_mz + self.proton_mass) / 2
                    # self.known_ions [f"{ion}({pair_acid_1}{pair_acid_2})^2"]  = [pair_mz, False]
                    # considers isotope + doubly charged ions
                    # self.known_ions [f"{ion}({pair_acid_1}{pair_acid_2})+i^2"]  = [pair_mz + (self.isotope_mass / 2), False]

        # Triple nested for loops to identify trios of amino acids (i.e. A and A)
        # First amino acid
        for identifier_index_1 in range(len(amino_acids)):
            for ion in ion_types:
                trio_acid_1 = amino_acids[identifier_index_1]
                # Checks if it is Carbamidomethyl or Oxidation
                if len(trio_acid_1) == 1:
                    trio_1_acid_base = trio_acid_1
                    trio_mz_1 = mass.calculate_mass(sequence=trio_acid_1, ion_type=ion, charge=1)
                else:
                    if trio_acid_1[0] == '[':
                        trio_1_acid_base = trio_acid_1[8]
                    else:
                        trio_1_acid_base = trio_acid_1[0]
                    trio_mz_1 = amino_acid_modifications[trio_acid_1]['mz'] + mass.calculate_mass(sequence=trio_1_acid_base, ion_type=ion, charge=1)
                
                # Second amino acid
                for identifier_index_2 in range(len(amino_acids) - identifier_index_1):
                    identifier_index_2 += identifier_index_1
                    trio_acid_2 = amino_acids[identifier_index_2]
                    if len(trio_acid_2) == 1:
                        trio_2_acid_base = trio_acid_2
                        trio_mz_2 = mass.calculate_mass(sequence=trio_acid_2, ion_type='b', charge=0)
                    else:
                        if trio_acid_2[0] == '[':
                            trio_2_acid_base = trio_acid_2[8]
                        else:
                            trio_2_acid_base = trio_acid_2[0]
                        trio_mz_2 = amino_acid_modifications[trio_acid_2]['mz'] + mass.calculate_mass(sequence=trio_2_acid_base, ion_type='b', charge=0)
                    
                    # Third amino acid
                    for identifier_index_3 in range(len(amino_acids) - identifier_index_2):
                        identifier_index_3 += identifier_index_2
                        trio_acid_3 = amino_acids[identifier_index_3]
                        if len(trio_acid_3) == 1:
                            trio_3_acid_base = trio_acid_3
                            trio_mz_3 = mass.calculate_mass(sequence=trio_acid_3, ion_type='b', charge=0)
                        else:
                            if trio_acid_3[0] == '[':
                                trio_3_acid_base = trio_acid_3[8]
                            else:
                                trio_3_acid_base = trio_acid_3[0]
                            trio_mz_3 = amino_acid_modifications[trio_acid_3]['mz'] + mass.calculate_mass(sequence=trio_3_acid_base, ion_type='b', charge=0)

                        # Adds all 3 together
                        trio_mz = trio_mz_1 + trio_mz_2 + trio_mz_3
                        if trio_mz <= 400:
                            self.known_ions[f"{ion}{{{trio_acid_1}{trio_acid_2}{trio_acid_3}}}"] = [trio_mz, False]
                            # considers isotope
                            self.known_ions[f"{ion}{{{trio_acid_1}{trio_acid_2}{trio_acid_3}}}+i"] = [trio_mz + self.isotope_mass, False]

                        # considers doubly charged ions
                        # trio_mz = (trio_mz + self.proton_mass) / 2
                        # if trio_mz <= 400:
                            # self.known_ions[f"{ion}({trio_acid_1}{trio_acid_2}{trio_acid_3})^2"] = [trio_mz, False]
                            # considers isotope
                            # self.known_ions[f"{ion}({trio_acid_1}{trio_acid_2}{trio_acid_3})+i^2"] = [trio_mz + (self.isotope_mass / 2), False]

        # Consider water and ammonia losses. This does not account for duplicates yet!!
        water = mass.calculate_mass(formula='H2O')
        ammonia = mass.calculate_mass(formula='NH3')
        water_ammonia = water + ammonia
        additional_explanations = {}
        for key in self.known_ions:
            add_on = ""
            if key.endswith('^2'):
                add_on = "^2"
                truncated_key = key[0:-2]
                if truncated_key.endswith('+i'):
                    add_on = "+i^2"
                    truncated_key = truncated_key[0:-2]
                additional_explanations[self.simplify_explanation(truncated_key + '-H2O') + add_on] = [self.known_ions[key][0] - (water / 2), False]
                additional_explanations[self.simplify_explanation(truncated_key + '-NH3') + add_on] = [self.known_ions[key][0] - (ammonia / 2), False]
                additional_explanations[self.simplify_explanation(truncated_key + '-H2O-NH3') + add_on] = [self.known_ions[key][0] - (water_ammonia / 2), False]
            else:
                add_on = ""
                if key.endswith('+i'):
                    add_on = "+i"
                    truncated_key = key[0:-2]
                additional_explanations[self.simplify_explanation(key + '-H2O') + add_on] = [self.known_ions[key][0] - water, False]
                additional_explanations[self.simplify_explanation(key + '-NH3') + add_on] = [self.known_ions[key][0] - ammonia, False]
                additional_explanations[self.simplify_explanation(key + '-H2O-NH3') + add_on] = [self.known_ions[key][0] - water_ammonia, False]
                if key[0] == 'a' or key[0] == 'b' or 'K' in key:
                    additional_explanations[key + '+TMT'] = [self.known_ions[key][0] + 229.162932, False]

        for key in additional_explanations:
            self.known_ions[key] = additional_explanations[key].copy()

        # Adds other known theoretical ions with the formula
        # Change to a list of the formulas and change the + self.proton_mass to - self.electron_mass (0.000548), chemical formula is the
        # same. Convention is proton is around the same as a hydrogen, loop through the list and add what I need to add
        additional_explanations = [
            "C3H9N2O2",
            "CH6N4O2",
            "C5H13N2O",
            "C4H11N2O2", # after this is plasma hardcoded
            "C2H9N3O4",
            "C2H5N10",
            "C5H12N3O2",
            "C5H12N3O",
            "C6H11N4O",
            "C10H10N",
            "C4H10N3O2",
            "C11H7N2O",
            "C8H16NO",
            "C6H9O4",
            "C11H15",
            "C5H5N7",
            "C7H14N3O",
            "C9H15",
            "C6H16N5O",
            "C5H10N3O",
            "C6H15N2O",
            "C5H13N3O2",
            "C10H8NO2",
            "C6H13O3",
            "C4H8N3O",
            "C8H7O4",
            "C11H9N2O",
            "C11H14NO2",
            "C9H9O5",
            "C6H7O3",
            "C9H7O4",
            "C10H16N3",
            "C9H14N3",
            "C8H11O",
            "C7H11O4",
            "C7H5O2",
            "C3HO7",
            "C5H12N3",
            "C5H12NO",
            "C9H12O4",
            "C8H13",
            "C6H13N2",
            "C7H12N",
            "C3H3N6",
            "C7H13O2",
            "C4H5O3",
            "C6H14N",
            "C4H10NO2",
            "C12H6N2", # after this is Q2018 hardcoded
            'C5H11', # (cyclopentane)
            'CH5N4',
            'C4H9O',
            'C3H7O2',
            'C2H7N2O',
            'C4H7O2',
            'C4H9O2',
            'H6N6',
            'C5H2NO',
            'C3HNO3',
            'C6H7O2',
            'C6H9O3',
            'C6H11O3',
            'C8H9O2',
            'C8H5O3',
            'C8H11O3',
            'C5H11O6',
            'C6H2N5O2',
            'C8H17O4',
            'C6H4N5O3',
            'C6H4N5O4',
            'C6H6N5O4',
            'C4N11O2',
            "H4PO4"
        ]

        for formula in additional_explanations:
            self.known_ions[formula] = [mass.calculate_mass(formula=f"{formula}") - self.electron_mass, False]
            # print(self.known_ions[formula])

        # Monosaccharide oxonium ions
        self.known_ions["HexNAc-frag1"] = [126.05495, False]
        self.known_ions["Pentose"] = [133.0496, False]
        self.known_ions["HexNAc-frag2"] = [138.05495, False]
        self.known_ions["HexNAc-frag3"] = [144.0661, False]
        self.known_ions["Hexose-H2O"] = [145.04953, False]
        self.known_ions["Deoxyhexose"] = [147.0652, False]
        self.known_ions["Hexose"] = [163.0601, False]
        self.known_ions["HexNAc-2H2O"] = [168.0655, False]
        self.known_ions["HexNAc-H2O"] = [186.0761, False]
        self.known_ions["HexNAc"] = [204.08666, False]
        self.known_ions["Neu5Ac-2H2O"] = [256.08156, False]
        self.known_ions["Neu5Ac-H2O"] = [274.0921, False]
        self.known_ions["Neu5Ac"] = [292.1026915, False]   # Also NeuNAc
        self.known_ions["NeuGc"] = [308.0976, False]
        self.known_ions["HexHex"] = [325.1129, False]
        self.known_ions["Oxonium-NH"] = [366.1397, False]

        # TMT related ions
        self.known_ions["TMT126"] = [126.127726, False]
        self.known_ions["TMT127N"] = [127.124761, False]
        self.known_ions["TMT127C"] = [127.131081, False]
        self.known_ions["TMT128N"] = [128.128116, False]
        self.known_ions["TMT128C"] = [128.134436, False]
        self.known_ions["TMT129N"] = [129.131471, False]
        self.known_ions["TMT129C"] = [129.137790, False]
        self.known_ions["TMT130N"] = [130.134825, False]
        self.known_ions["TMT130C"] = [130.141145, False]
        self.known_ions["TMT131N"] = [131.138180, False]
        self.known_ions["TMT131C"] = [131.1445, False]
        self.known_ions["TMT132N"] = [132.141535, False]
        self.known_ions["TMT132C"] = [132.147855, False]
        self.known_ions["TMT133N"] = [133.14489, False]
        self.known_ions["TMT133C"] = [133.15121, False]
        self.known_ions["TMT134N"] = [134.148245, False]
        self.known_ions["TMT134C"] = [134.154565, False]
        self.known_ions["TMT135N"] = [135.151600, False]
        self.known_ions["TMT6plex"] = [229.162932 + self.proton_mass, False]
        self.known_ions["TMT6plex+H2O"] = [229.162932 + 18.010565 + self.proton_mass, False]
        self.known_ions["TMT6pro"] = [304.207146 + self.proton_mass, False]
        self.known_ions["TMT6pro+H2O"] = [304.207146 + 18.010565 + self.proton_mass, False]
        self.known_ions["TMT6pro+NH3"] = [304.207146 + 17.026549 + self.proton_mass, False]

        for key in self.known_ions:
            self.known_ions[key][0] = round(self.known_ions[key][0], 5)

    def simplify_explanation(self, explanation):
        # This takes a string and splits it into the + and -, reducing any redundancies
        # i.e. IE+H2O-H2O -> IE
        condensed_acid = ""
        match = re.match(r'([abyI\[\]][^\+\-]+)(.+)$', explanation)
        if match:
            ions = match.group(1)
            losses = match.group(2)
                
        match = re.findall(r'([+-][A-Z\d]+)', losses)

        if match:
            additions = []
            subtractions = []

            for modification in match:
                if modification[0] == "+":
                    additions.append(modification[1:])
                else:
                    subtractions.append(modification[1:])

            additions.sort()
            subtractions.sort()

            removed = 0
            for add_index in range(len(additions)):
                for subtract_index in range(len(subtractions)):
                    if additions[add_index - removed] == subtractions[subtract_index - removed]:
                        additions.pop(add_index - removed)
                        subtractions.pop(subtract_index - removed)
                        removed += 1
                        break
        
            condensed_acid = ions

            for add_modification in additions:
                condensed_acid += "+" + add_modification

            for subtract_modification in subtractions:
                condensed_acid += "-" + subtract_modification
        else:
            condensed_acid = explanation

        return condensed_acid

    def determine_crude_correction(self):
        # Uses popular possible explanations so in the first run through, it looks for the most 
        # intense peak within 20 PPM of each popular possible explanation. Then, it calculates the
        # delta between the theoretical mass with the measured mass, to calculate a crude correction
        first_pass_peaks = {
            'IP': 70.06512,
            'IV': 72.08078,
            'IQ-NH3': 84.04439,
            'IK-NH3': 84.08077,
            'IL': 86.09643,
            'IQ': 101.07094,
            'IE': 102.05495,
            'IM': 104.05285,
            'IH': 110.07127,
            'IR-NH3': 112.08692,
            'a{{AA}}': 115.08659,
            'IR+H2O+H2O-N3H7': 116.0706,
            'IF': 120.08078,
            'b{{AA}}-NH3': 126.05495,
            'IQ+CO': 129.06585,
            'IK+CO+H2O-NH3': 130.08625,
            'IC[Carbamidomethyl]': 133.04301,
            'IY': 136.07569,
            'a{{AP}}': 141.10224,
            'b{{AA}}': 143.0815,
            'a{{AV}}': 143.11789,
            'b{{GP}}': 155.0815,
            'IK+CO+H2ON2-NH3': 158.0924,
            'IW': 159.09167,
            'a{{DP}}-H2O': 167.0815,
            'b{{AP}}': 169.09715,
            'a{{PV}}': 169.13354,
            'a{{PT}}': 171.1128,
            'a{{TV}}': 173.12845,
            'y{{R}}': 175.11895,
            'a{{LP}}': 183.14919,
            'b{{PS}}': 185.09207,
            'b{{PV}}': 197.12845,
            'a{{DL}}': 201.12337,
            'a{{AY}}': 207.11281,
            'y{{AH}}-H2O': 209.10331,
            'a{{EL}}': 215.13902,
            'b{{EL}}-H2O': 225.12337,
            'a{{APS}}': 228.13427,
            'b{{DL}}': 229.11828,
            'y{{HV}}-H2O': 237.1346,
            'a{{APT}}': 242.14992,
            'b{{DQ}}': 244.0928,
            'y{{KP}}': 244.16557,
            'y{{HV}}': 255.14517,
            'y{{PR}}': 272.17172
        }

        self.crude_xy_scatterplot = []
        first_pass_ppms = []
        for key in first_pass_peaks:
            # start_time = timeit.default_timer()
            peak_mz = first_pass_peaks[key]
            delta = self.initial_tolerance * first_pass_peaks[key] / 1e6
            lower_bound = peak_mz - delta
            upper_bound = peak_mz + delta
            best_match = [0, 0, 0]
            i_lower_bound = int(lower_bound * 10000)
            i_upper_bound = int(upper_bound * 10000)
            for i in range(i_lower_bound, i_upper_bound): # double check this is the right range
                if self.by_count[i] > best_match[1]:
                    best_match = [i / 10000.0, self.by_count[i], (i/10000.0 - peak_mz) * 1e6 / peak_mz]
            
            if best_match[0] != 0:
                self.crude_xy_scatterplot.append(best_match)
                first_pass_ppms.append(best_match[2])
            # print(f"time to identify {key}: {start_time - timeit.default_timer()}")

        # Crude correction is the median of all the delta values between the theoretical mz and
        # measured mz value
        self.crude_correction = np.median(first_pass_ppms)

    def increase_initial_tolerance(self):
        self.initial_tolerance = 30

    def lower_minimum_triggers(self):
        # If not enough peaks were over the minimum trigger limit, it lowers the minimum triggers by 15
        self.minimum_triggers = self.minimum_triggers - 15

    def find_initial_triggers(self):
        print("starting search for initial triggers")

        # Scans through peaks to find all peaks with over the minimum triggers value
        # If it is over the minimum triggers value (which compares it to the total strength of the peak),
        # it will add it to triggered peaks
        previous_peak = [0, 0]
        counted = False
        self.triggered_peaks = []
        self.observed_peaks = []
        for i in range(len(self.by_strength)):
            if self.by_strength[i] >= self.minimum_triggers:
                if self.by_strength[i] < previous_peak[1] and not counted:
                    if self.by_strength[i + 1] < self.by_strength[i]:
                        self.triggered_peaks.append(previous_peak)
                        counted = True
                previous_peak = [i/10000, self.by_strength[i]]
            else:
                if not counted and previous_peak[0] != 0 and self.by_strength[i + 1] < self.by_strength[i]:
                    self.triggered_peaks.append(previous_peak)
                counted = False
                previous_peak = [0, 0]

    def refine_triggered_peaks(self):
        print("starting to refine triggered peaks")

        # Refines the peaks by using a gaussian fit on each triggered peak. If the peaks fit centers
        # are within 5 PPM, then the two peaks are considered the same peak. This removes duplicates
        for peak in self.triggered_peaks:
            peak_and_intensity = self.get_peak_fit_center(peak)
            if peak_and_intensity[0] == 0 and peak_and_intensity[1] == 0:
                continue
            else:
                peak_and_intensity = self.get_peak_fit_center(peak_and_intensity)
                if peak_and_intensity[0] == 0 and peak_and_intensity[1] == 0:
                    continue
                else:
                    if len(self.observed_peaks) > 0:
                        mz_delta = self.observed_peaks[-1][0] * 5 / 1e6
                        upper_bound = self.observed_peaks[-1][0] + mz_delta
                        lower_bound = self.observed_peaks[-1][0] - mz_delta
                        if peak_and_intensity[0] >= lower_bound and peak_and_intensity[0] <= upper_bound:
                            continue
                        else:
                            self.observed_peaks.append(peak_and_intensity)
                    else:
                        self.observed_peaks.append(peak_and_intensity)
            
    def get_peak_fit_center(self, peak):
        # This is a function that can calculate the peak fit center (center of the applied gaussian fit)
        mz_values = []
        intensity_values = []
        ppm_values = []
        mz_delta = int((self.ppm_delta / 1e6) * peak[0] * 10000)
        # Since the data above m/z of 400 was not stored, cannot find the peak fit center of peaks if the higher end of the gaussian fit is over 400
        if int(peak[0] * 10000 + mz_delta) > 4000000:
            return [0, 0]
    
        for index in range(mz_delta * 2 + 1):
            # change it to be 0
            add_index = int(peak[0] * 10000 + index - mz_delta - 1)
            mz_values.append(add_index / 10000 - peak[0])
            ppm_values.append((add_index / 10000 - peak[0]) * 1e6 / peak[0])
            intensity_values.append(self.by_strength[add_index])

        n = len(ppm_values)
        center = int(n/2)
        binsize = ppm_values[center]-ppm_values[center-1]

        try:
        #if 1:
            popt,pcov = curve_fit(gaussian_function,ppm_values,intensity_values,p0=[intensity_values[center],ppm_values[center],binsize])
        except:
            return [0, 0]

        peak_mz = round((peak[0] + (popt[1] * peak[0] / 1e6)), 5)
        return [peak_mz, round(popt[0], 1)]

    def binned(self):
        binned_known_ions = [ [] for i in range(4001) ]
        for key in self.known_ions:
            bin = int(self.known_ions[key][0] * 10)
            if bin < 4001:
                binned_known_ions[bin].append([key, self.known_ions[key][0], False])
        
        return binned_known_ions

    def identify_peaks(self):
        print("starting to identify peaks")
        # This uses the observed peaks and sees how many of them can be identified. It accounts for any
        # corrections applicable

        binned_known_ions = self.binned()

        # Creates an array with a spline correction for each observed peak
        if self.has_correction_spline:
            mz_values = []
            for peak in self.observed_peaks:
                mz_values.append(peak[0])
            x_values = np.array(mz_values)
            y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)
            if self.has_second_spline:
                y_values_2 = interpolate.BSpline(self.t2, self.c2, self.k2)(x_values)
                for index in range(len(y_values)):
                    y_values[index] += y_values_2[index]
        
        for index in range(len(self.observed_peaks)):
            identifications = []
            # applies the corrections
            peak = self.observed_peaks[index].copy()
            crude_correction_mz = self.crude_correction * peak[0] / 1e6

            if self.has_correction_spline:
                spline_correction_mz = y_values[index] * peak[0] / 1e6
            else:
                spline_correction_mz = 0

            peak[0] -= (crude_correction_mz + spline_correction_mz)
            self.observed_peaks[index] = self.observed_peaks[index][0:2]

            if peak[0] <= 150:
                peak_tolerance = self.tolerance_150 * peak[0] / 1e6
            else:
                peak_tolerance = self.tolerance * peak[0] / 1e6

            if int(peak[0] * 10) <= 4000:
                for known_ion_index in range(len(binned_known_ions[int(peak[0] * 10)])):
                    possible_explanation = binned_known_ions[int(peak[0] * 10)][known_ion_index]
                    amino_acid_mz = possible_explanation[1]
                    amino_acid = possible_explanation[0]
                    if not possible_explanation[2] and peak[0] > amino_acid_mz - peak_tolerance and peak[0] < amino_acid_mz + peak_tolerance:
                        identified_peak = [round(amino_acid_mz, 5), round(peak[0] - amino_acid_mz, 8), amino_acid]
                        identifications.append(identified_peak)
                        binned_known_ions[int(peak[0] * 10)][known_ion_index][2] = True
    
                # Sorts the identifications for each peak so the closest identifications are at the front
                identifications.sort(key = lambda x: len(x[2]))
                identifications.sort(key = lambda x: abs(x[1]))
                for identification in identifications:
                    self.observed_peaks[index].append(identification)

    def update_refined(self, crude_correction, t, c, k):
        mz_values = []
        for peak in self.refined_mz_values:
            mz_values.append(peak)
        x_values = np.array(mz_values)
        y_values = interpolate.BSpline(t, c, k)(x_values)

        removed = 0
        for index in range(len(self.refined_delta_ppm_values)):
            if abs(self.refined_delta_ppm_values[index - removed] - y_values[index] - crude_correction) > 10:
                del self.refined_mz_values[index - removed]
                del self.refined_delta_ppm_values[index - removed]
                removed += 1
            else:
                mz = self.refined_mz_values[index - removed]
                self.refined_mz_values[index - removed] -= (crude_correction + y_values[index]) * mz / 1e6
                self.refined_delta_ppm_values[index - removed] -= (y_values[index] + crude_correction)

    def keep_intense_remove_outliers(self, mz_values, delta_ppm_values):
        print("starting to remove outliers")
        # Creates a new variable to store a set of refined peaks, only keeping the most intense peaks in
        # a bin and removes the outliers (top/bottom 10%) in the bin
        bin_size = 25 # Uses a bin size of 25 to take outliers and top intense peaks
        xy_range = []
        xy = {}
        removable_index = []
        upper_bound = mz_values[0][0] + bin_size
        for index in range(len(mz_values)):
            xy_range.append([mz_values[index][0], mz_values[index][1], delta_ppm_values[index]])
            xy[mz_values[index][0]] = delta_ppm_values[index]
            if mz_values[index][0] > upper_bound or index == len(mz_values) - 1:
                xy_range.sort(key = lambda x: -1 * x[1])
                if len(xy_range) != 1:
                    top_peak = int(self.top_peaks * len(mz_values) + 0.5) + 1
                    xy_range = xy_range[0:top_peak]

                xy_range.sort(key = lambda x: x[2])
                length = len(xy_range)
                i = 1 / length
                if length >= 10:
                    for tuple in xy_range:
                        if i <= 0.1 or i >= 0.9:
                            removable_index.append(tuple[0])
                        i += 1 / length
                elif length >= 5:
                    for tuple in xy_range:
                        if i <= 0.2 or i >= 0.8:
                            removable_index.append(tuple[0])
                        i += 1 / length
                else:
                    for tuple in xy_range:
                        removable_index.append(tuple[0])
                upper_bound += bin_size
                xy_range = []
        
        for key in removable_index:
            xy.pop(key)

        # Creates the new array with mz values and delta ppm values
        self.refined_mz_values = []
        self.refined_delta_ppm_values = []
        for key in xy:
            self.refined_mz_values.append(key)
            self.refined_delta_ppm_values.append(xy[key])

    def plot_crude_calibrations(self):
        print("starting to plot scatterplots")
        # Plots the top left graph, which is the PPM (difference between theoretical ion masses and the
        # most intense peak within a 20 PPM range  vs the mz of the ion
        self.pdf = PdfPages(f'{self.peak_file}.calibration.pdf')
        self.fig, self.ax = plt.subplots(nrows=3, ncols=2, figsize=(8.5, 11))
        peak_mz = []
        calibration_ppm = []
        for closest_match in self.crude_xy_scatterplot:
            peak_mz.append(closest_match[0])
            calibration_ppm.append(closest_match[2])
        self.ax[0][0].plot(peak_mz, calibration_ppm, '.',c="g", markersize=2)
        self.ax[0][0].axhline(y = self.crude_correction, color = 'k', linestyle = 'dashed')
        self.ax[0][0].axhline(y = 0, color = 'k', linewidth = 1, linestyle = '-')
        self.ax[0][0].text(min(peak_mz), 0, f'tolerance: {self.initial_tolerance} ppm', fontsize='xx-small', ha='left', va='top')
        self.ax[0][0].set(xlabel='m/z', ylabel='PPM')
        plt.tight_layout()

    def delta_scatterplots(self):
        # Plots the identified peaks PPM vs the m/z without any corrections applied. The crude calibration
        # is plotted as a dashed line
        mz_values = []
        delta_values_ppm = []
        for peak in self.observed_peaks:
            if len(peak) >= 3:
                mz_values.append([peak[0], peak[1]])
                delta_values_ppm.append(((peak[2][1]) * 1e6 / peak[2][0]) + self.crude_correction)
            else:
                continue
        self.keep_intense_remove_outliers(mz_values, delta_values_ppm)
        if self.make_pdf:
            # self.ax[0,1].set_title(f"crude_correction: {self.crude_correction}", fontsize="xx-small")
            self.ax[0][1].scatter(self.refined_mz_values, self.refined_delta_ppm_values, 0.5)
            self.ax[0][1].axhline(y = self.crude_correction, color = 'k', linestyle = 'dashed')
            mz_values = [peak[0] for peak in mz_values]
            self.ax[0][1].text(min(mz_values), 0, f'tolerance: {self.tolerance} ppm', fontsize='xx-small', ha='left', va='top')
            self.ax[0][1].text(315, self.crude_correction + 0.1, f'crude_correction: {round(self.crude_correction, 2)} ppm', fontsize='xx-small', ha='right', va='bottom')
            self.ax[0][1].axhline(y = 0, color = 'k', linewidth = 1, linestyle = '-')
            self.ax[0][1].set(xlabel='m/z', ylabel='PPM')

            self.ax[2][1].scatter(mz_values, delta_values_ppm, 0.5)
            self.ax[2][1].axhline(y = 0, color = 'k', linewidth = 1, linestyle = '-')
            self.ax[2][1].set(xlabel='m/z', ylabel='PPM')
            plt.tight_layout()

    def find_first_spline(self):
        # Calculates the artificial peak point for the after 400 calibration constant
        mz_values = []
        delta_values_ppm = []
        artificial_peak_range = []
        for index in range(len(self.refined_mz_values)):
            peak = [self.refined_mz_values[index], self.refined_delta_ppm_values[index]]
            crude_correction_mz = self.crude_correction * peak[0] / 1e6
            mz_values.append(peak[0] + crude_correction_mz)
            delta_values_ppm.append(peak[1] - self.crude_correction)
            if peak[0] >= 300 and peak[0] <= 400:
                artificial_peak_range.append(peak[1] - self.crude_correction) 
        
        mz_values.sort()
        index = mz_values[-1]
        num_of_artificial_peaks = max(100, len(artificial_peak_range))
        index_increment = (500 - index) / num_of_artificial_peaks
        artificial_peak_intensity = np.median(artificial_peak_range)
        self.after_400_calibration = np.median(artificial_peak_range)
        if self.make_pdf:
            self.ax[1][0].axhline(y = 0, color = 'k', linewidth = 1, linestyle = '-')
            self.ax[1][0].set(xlabel='m/z', ylabel='PPM')

        # Adds artificial peaks above 500 with the median between 300 and 400 so the spline 
        # fit will flatten out, making the correction applicable for points above 400 if necessary

        if len(mz_values) <= 5:
            print("There are too little data points to create a spline fit.")
            self.ax[1,0].set_title("There are too little data points to create a spline fit", fontsize="xx-small")
        else:
            try: # Tries to apply a spline fit without the extension line
                x_new = np.linspace(0, 1, 5)[1:-1] # 5 = 3 knots + 2
                x_new.sort()
                q_knots = np.quantile(mz_values, x_new)
                self.t, self.c, self.k = interpolate.splrep(mz_values, delta_values_ppm, t=q_knots, s=3)
                self.has_correction_spline = True
                x_regular = np.arange(mz_values[0], mz_values[-1])
                y_values = interpolate.BSpline(self.t, self.c, self.k)(x_regular)
            except: # If the spline was not fit, add an extension line to fit the spline
                print("spline fit 1 without extension line failed, extension line applied")
                while index < 500: # extending the line so the spline fit is created with the artificial
                # peak plots
                    index += index_increment
                    mz_values.append(index)
                    delta_values_ppm.append(artificial_peak_intensity)

                x_new = np.linspace(0, 1, 5)[1:-1] # 5 = 3 knots + 2
                x_new.sort()
                q_knots = np.quantile(mz_values, x_new)
                self.t, self.c, self.k = interpolate.splrep(mz_values, delta_values_ppm, t=q_knots, s=3)
                self.has_correction_spline = True
                x_regular = np.arange(mz_values[0], mz_values[-1])
                y_values = interpolate.BSpline(self.t, self.c, self.k)(x_regular)
            if self.make_pdf:
                self.ax[1][0].plot(x_regular, y_values, 'b')
                # self.ax[1,0].set_title(f"t: {self.t}, c: {self.c}, k: {self.k}", fontsize="xx-small")
            self.update_refined(self.crude_correction, self.t, self.c, self.k)
        if self.make_pdf:
            self.ax[1][0].scatter(mz_values, delta_values_ppm, 0.5)
            plt.tight_layout()
            plt.close()

    def find_second_spline(self):
        # This calculates the coefficients for the second spline fit and plots the delta scatterplot
        # of just the first spline fit applied, with the second spline fit overlayed
        mz_values = self.refined_mz_values.copy()
        delta_values_ppm = self.refined_delta_ppm_values.copy()

        # extend out the line for the spline
        artificial_peak_range = []
        for index in range(len(mz_values)):
            peak = [mz_values[index], delta_values_ppm[index]]
            if peak[0] >= 300 and peak[0] <= 400:
                artificial_peak_range.append(peak[1]) 

        self.after_400_calibration += np.median(artificial_peak_range)
        
        # create the spline
        try:
            x_new = np.linspace(0, 1, 7)[1:-1] # 7 = 5 knots + 2
            x_new.sort()
            q_knots = np.quantile(mz_values, x_new)
            self.t2, self.c2, self.k2 = interpolate.splrep(mz_values, delta_values_ppm, t=q_knots, s=3)
            self.has_second_spline = True
            x_regular = np.arange(mz_values[0], mz_values[-1])
            y_values = interpolate.BSpline(self.t2, self.c2, self.k2)(x_regular)
            if self.make_pdf:
                self.ax[1][1].plot(x_regular, y_values, 'g')
        except:
            try:
                print("adding artificial line on second spline fit")
                # Adding the artificial points to create an extension line
                index = max(mz_values)
                num_of_artificial_peaks = max(100, len(artificial_peak_range))
                index_increment = (500 - index) / num_of_artificial_peaks
                artificial_peak_intensity = np.median(artificial_peak_range)
                while index < 500:
                    index += index_increment
                    mz_values.append(index)
                    delta_values_ppm.append(artificial_peak_intensity)

                # Plotting spline fit
                x_new = np.linspace(0, 1, 7)[1:-1] # 7 = 5 knots + 2
                x_new.sort()
                q_knots = np.quantile(mz_values, x_new)
                self.t2, self.c2, self.k2 = interpolate.splrep(mz_values, delta_values_ppm, t=q_knots, s=3)
                self.has_second_spline = True
                x_regular = np.arange(mz_values[0], mz_values[-1])
                y_values = interpolate.BSpline(self.t2, self.c2, self.k2)(x_regular)
                if self.make_pdf:
                    self.ax[1][1].plot(x_regular, y_values, 'g')
            except:
                if self.make_pdf:
                    self.ax[1][1].text(max(mz_values), max(delta_values_ppm), 'second spline fit could not be applied', fontsize='xx-small', ha='right', va='top')
                print("the second spline fit could not be applied")
        finally:
            if self.make_pdf:
                self.ax[1][1].scatter(mz_values, delta_values_ppm, 0.5)
                self.ax[1][1].axhline(y = 0, color = 'k', linewidth = 1, linestyle = '-')
                self.ax[1][1].set(xlabel='m/z', ylabel='PPM')
            self.update_refined(0, self.t2, self.c2, self.k2)

    def plot_corrected_scatterplot(self):
        # This plots the fifth graph, with both spline corrections and the crude correction applied
        # Most of the points should be around the 0 horizontal line, which indicates the corrections
        # were successful        
        self.ax[2][0].scatter(self.refined_mz_values, self.refined_delta_ppm_values, 0.5)
        self.ax[2][0].axhline(y = 0, color = 'k', linewidth = 1, linestyle = '-')
        self.ax[2][0].set(xlabel='m/z', ylabel='PPM')

    def plot_all_corrections(self):
        x_values = np.arange(self.refined_mz_values[0], self.refined_mz_values[-1])
        
        y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)
        for index in range(len(y_values)):
            y_values[index] += self.crude_correction
        if self.has_second_spline:
            y_values_2 = interpolate.BSpline(self.t2, self.c2, self.k2)(x_values)
            for index in range(len(y_values)):
                y_values[index] += y_values_2[index]
        
        self.ax[2][1].plot(x_values, y_values, 'g')
        total_after_400_correction = self.after_400_calibration + self.crude_correction
        self.ax[2][1].plot([400, 450], [total_after_400_correction, total_after_400_correction], linewidth = 1, color = 'g')
        self.fig.subplots_adjust(top=0.96, bottom=0.10, left=0.10, right=0.96, hspace=0.2, wspace=0.2)
        plt.tight_layout()
        self.pdf.savefig(self.fig)

    def analyze_snippets(self):
        print("starting to analyze snippets")
        # Send a message if the plots cannot be created
        if len(self.scan_snippets) == 0:
            print("not enough data points to create snippet plots")
            return
        
        acid_mz = {'IH': 110.07127, 'IF': 120.08078, 'IK-CO': 129.10223}
        close_matches = []

        # Creates enough elements to store all the scan numbers, delta PPM, and intensities in one
        # array per acid. The first element is the acid, the rest are the stored values
        for key in acid_mz:
            close_matches.append([key])

        # Goes through each snippet
        for snippet in self.scan_snippets:
            closest_PPM_and_logged_intensity = []
            # It looks for the most intense match within the snippet and adds the PPM and intensity to the
            # array. This will go through all the acids and add it to the result arrays so they can
            # be plotted against each other
            for key in acid_mz:
                closest_index = 0
                largest_intensity = 0
                mz = acid_mz[key]
                tolerance = 5 * mz / 1e6
                upper_bound = mz + tolerance
                lower_bound = mz - tolerance
                # snippet[0] is the scan number, snippet[1] is the range of mz and snippet[2] 
                # is the range of respective intensities
                for index_snippet_mz in range(len(snippet[1])):
                    snippet_mz = snippet[1][index_snippet_mz]

                    # Applies calibration before the tolerance is checked
                    correction_mz = self.crude_correction * snippet_mz / 1e6
            
                    if self.has_correction_spline:
                        y_values = interpolate.BSpline(self.t, self.c, self.k)([snippet_mz])
                        correction_mz += y_values[0] * snippet_mz / 1e6
                    snippet_mz -= correction_mz

                    if snippet_mz > lower_bound and snippet_mz < upper_bound:
                        if snippet[2][index_snippet_mz] > largest_intensity:
                            largest_intensity = snippet[2][index_snippet_mz]
                            closest_index = index_snippet_mz
                    elif snippet_mz > upper_bound:
                        break
                
                # Either add the PPM and intensity or if it could not find the acid, add [0, 0]
                # Adding [0, 0] ensures the indices of the closest_PPM_and_intensity and the results array
                # will line up, so the values are added to the right element in the array
                # Ex: It will make sure IK+CO is added to IK+CO element, even if IH and IF were not found
                if largest_intensity != 0:
                    observed_mz = snippet[1][closest_index]
                    # Also accounts for the correction
                    correction_mz = self.crude_correction * observed_mz / 1e6
            
                    if self.has_correction_spline:
                        y_values = interpolate.BSpline(self.t, self.c, self.k)([observed_mz])
                        correction_mz += y_values[0] * observed_mz / 1e6
                    closest_PPM_and_logged_intensity.append([(observed_mz - mz - correction_mz) * 1e6 / mz, math.log(largest_intensity)])
                else:
                    closest_PPM_and_logged_intensity.append([0, 0])

            # rework naming convention!
            # Adds the scan number, delta PPM, and intensity to the respective array so it is ready to be
            # plotted later.
            for index in range(len(closest_PPM_and_logged_intensity)):
                possible_acid = closest_PPM_and_logged_intensity[index]
                if possible_acid[1] != 0:
                    close_matches[index].append([snippet[0], possible_acid[0], possible_acid[1]])

        # Creates the subplots for the output
        # It will not have a good output if there are too many acids!! Will need to rework so there are
        # not too many plots on one page
        self.fig, self.ax = plt.subplots(nrows=len(acid_mz), ncols=2, figsize=(8.5, 11))
        self.fig.subplots_adjust(top=0.96, bottom=0.05, left=0.12, right=0.96, hspace=0.2, wspace=0.2)

        # Each acid will have its own row, one scatter plot being the scan number vs PPM, the other being
        # the intensity vs PPM to analyze the trends
        for row in range(len(acid_mz)):
            # Sorts the matches by intensity and removes bottom 10%
            sorted_matches = close_matches[row][1:].copy()
            sorted_matches.sort(key = lambda x: x[0]) # sorts sorted_matches by scan number

            upper_scan_num = 1000
            intense_sorted_matches = [] # only keeps the intense matches
            sorted_matches_bin = [] # takes all matches between 1-1000, 1001-2000, etc
            for index in range(len(sorted_matches)):
                match = sorted_matches[index]
                if match[0] > upper_scan_num or index == len(sorted_matches):
                    upper_scan_num += 1000
                    sorted_matches_bin.sort(key = lambda x:-x[2]) # sorts by negative intensity
                    sorted_matches_bin = sorted_matches_bin[0:20] # takes top 20 most intensee
                    for sorted_match in sorted_matches_bin:
                        intense_sorted_matches.append(sorted_match)
                    sorted_matches_bin = []
                
                sorted_matches_bin.append(match)

            # Extracts arrays of scan numbers, delta PPMs, and logged intensities for plotting
            scan_nums = [item[0] for item in intense_sorted_matches]
            delta_PPMs = [item[1] for item in intense_sorted_matches]
            logged_intensities = [item[2] for item in intense_sorted_matches]
            print(f"{close_matches[row][0]}: {len(scan_nums)}")
            for col in range(2):
                if col == 0:
                    # This adds additional information like the axis titles, a line at 0, and a title
                    self.ax[row][col].plot(scan_nums, delta_PPMs, '.',c="g", markersize=2)
                    self.ax[row][col].axhline(y = 0, color = 'k', linewidth = 1, linestyle = '-')
                    self.ax[row][col].set(xlabel='scan_num', ylabel='delta PPM')
                    # self.ax[row,col].set_title(f"{close_matches[row][0]} scan num vs delta PPM", fontsize="xx-small")
                if col == 1:
                    # This adds additional information like the axis titles, a line at 0, and a title
                    self.ax[row][col].plot(logged_intensities, delta_PPMs, '.',c="g", markersize=2)
                    self.ax[row][col].axhline(y = 0, color = 'k', linewidth = 1, linestyle = '-')
                    self.ax[row][col].set(xlabel='logged intensity', ylabel='delta PPM')
                    # self.ax[row,col].set_title(f"{close_matches[row][0]} intensity vs delta PPM", fontsize="xx-small")

        # This saves the plots and adds it to the output file
        self.pdf.savefig(self.fig)

    def find_peak_percentage(self):
        print("starting to find percentages")
        for index in range(len(self.observed_peaks)):
            peak = self.observed_peaks[index]
            lower_mz = int(peak[0] * 10000 - 50000 * peak[0] / 1e6 + 0.5)
            peak_search_range = int((5 * peak[0] / 1e6) * 20000) + 1

            spectra_with_peak = 0
            for mz in range(peak_search_range):
                if mz + lower_mz >= 4000000:
                    break
                spectra_with_peak += self.by_count[lower_mz + mz]

            self.observed_peaks[index].insert(2, f"{round(spectra_with_peak * 100 / self.stats['ms2spectra'], 1)}%")

    def write_json(self):
        print("starting to write json file")
        # This writes out a json file with the crude correction value, the t, c, and k values for the
        # spline fit
        # https://www.geeksforgeeks.org/reading-and-writing-json-to-a-file-in-python/
        correction_values = {
            "crude correction": self.crude_correction,
            "after_400_calibration": self.after_400_calibration
        }

        if self.has_correction_spline:
            correction_values["t1"] = self.t.tolist()
            correction_values["c1"] = self.c.tolist()
            correction_values["k1"] = self.k

        if self.has_second_spline:
            correction_values["t2"] = self.t2.tolist()
            correction_values["c2"] = self.c2.tolist()
            correction_values["k2"] = self.k2

        # Serializing json
        json_object = json.dumps(correction_values, indent=4)
        
        # Writing to sample.json
        with open(f"{self.peak_file}.calibration.json", "w") as outfile:
            outfile.write(json_object)

    def plot_peaks_strength(self):
        print("starting to plot individual peaks")
        # Prints out individual peak plots, with 15 PPM above and below each observed peak to see the
        # intensity of the mz values around each observed peak
        # There are 15 plots per page

        x = 0
        y = 0
        fig, ax = plt.subplots(self.plot_output_rows, self.plot_output_columns, figsize=(8.5,11))

        if self.has_correction_spline:
            mz_values = [peak[0] for peak in self.observed_peaks]
            x_values = np.array(mz_values)
            y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)
            if self.has_second_spline:
                y_values_2 = interpolate.BSpline(self.t2, self.c2, self.k2)(x_values)
                for index in range(len(y_values)):
                    y_values[index] += y_values_2[index]

        # Takes the most 300 most intense peaks to print out
        if len(self.observed_peaks) > 300:
            self.observed_peaks.sort(key = lambda x: -1 * x[2])
            self.observed_peaks = self.observed_peaks[0:300]
            self.observed_peaks.sort(key = lambda x: x[0])

        for index in range(len(self.observed_peaks)):
            peak = self.observed_peaks[index]
            peak.pop(1)
            fig.subplots_adjust(top=0.98, bottom=0.05, left=0.12, right=0.98, hspace=0.2, wspace=0.38)
            # mz_values = []
            intensity_values = []
            ppm_values = []
            crude_correction_mz = self.crude_correction * peak[0] / 1e6
            mz_delta = int((self.ppm_delta / 1e6) * (peak[0] + crude_correction_mz) * 10000)
            for range_index in range(mz_delta * 2 + 1):
                # change it to be 0
                add_index = int((peak[0]) * 10000 + range_index - mz_delta)
                # mz_values.append(add_index / 10000 - peak[0])
                ppm = (add_index / 10000 - peak[0]) * 1e6 / peak[0]
                intensity_values.append(self.by_strength[add_index])
                ppm_values.append(ppm)


            ax[x,y].step(ppm_values, intensity_values, where='mid')
            ax[x,y].tick_params(axis='x', labelsize='xx-small')
            ax[x,y].locator_params(axis='x', nbins=5)
            # ax[x,y].set_xticks([-0.0015, -0.0007, 0, 0.0007, 0.0015])
            ax[x,y].tick_params(axis='y', labelsize='xx-small')
            ax[x,y].set_xlim(-15, 15)

            # add gaussian fitting
            n = len(ppm_values)
            center = int(n/2)
            binsize = ppm_values[center]-ppm_values[center-1]

            try:
            #if 1:
                popt,pcov = curve_fit(gaussian_function,ppm_values,intensity_values,p0=[intensity_values[center],ppm_values[center],binsize])
                ax[x,y].plot(ppm_values,gaussian_function(ppm_values, *popt),'r:')
            except:
                # either skip it or show it, for now just skip it
                y += 1
                y %= self.plot_output_columns
                if y == 0:
                    x += 1
                    x %= self.plot_output_rows
                    if x == 0:
                        self.pdf.savefig(fig)
                        plt.close('all')
                        fig, ax = plt.subplots(self.plot_output_rows, self.plot_output_columns,figsize=(8.5,11))
                continue
            
            if len(peak) >= 4:
                # if self.has_correction_spline:
                    # spline_correction_mz = y_values[index] * peak[0] / 1e6
                # else:
                spline_correction_mz = 0
                identified_ion_name = ''
                ion_mz = peak[3][0]
                count = 3
                plotted_identification = [0, 0, 0]
                while count <= len(peak) - 1:
                    if peak[count][1] != plotted_identification[1]:
                        line_x = (peak[count][1]) # * 1e6 / peak[0]
                        plotted_identification = [peak[count][2], peak[count][1], plotted_identification[2] + 1]
                        ax[x,y].axvline(x=-line_x, color='black', lw=1, linestyle='--')
                        ax[x][y].text(-line_x, max(intensity_values), f'{plotted_identification[2]}', fontsize='xx-small', ha='left', va='top')
                        identified_ion_name += f"{plotted_identification[2]}) "
                        identified_ion_name += peak[count][2] + '\n    '
                    elif count <= 10:
                        identified_ion_name += peak[count][2] + '\n    '
                    count += 1
                if len(peak) > 11:
                    identified_ion_name += "..."

                peak_fit_center = int(100000*(peak[0] + (popt[1] * peak[0] / 1e6) - crude_correction_mz - spline_correction_mz)) / 100000
                ax[x,y].text(-14, max(intensity_values) * 1.03, f'peak fit center: \n    {peak_fit_center}\nion m/z: \n    {ion_mz}\nion id: \n    {identified_ion_name}', fontsize='xx-small', ha='left', va='top')
                ax[x,y].text(14, max(intensity_values) * 1.03, f'{peak[2]}\nof spectra', fontsize='xx-small', ha='right', va='top')
                # ax[x,y].set_title(f"Peak {peak[0]}, Amino Acid {peak[4]}", fontsize="xx-small")
            else:
                peak_fit_center = int(100000*(peak[0] + (popt[1] * peak[0] / 1e6) - crude_correction_mz)) / 100000
                ax[x,y].text(-14, max(intensity_values) * 1.03, f'peak fit center: \n    {peak_fit_center}', fontsize='xx-small', ha='left', va='top',)
                ax[x,y].text(14, max(intensity_values) * 1.03, f'{peak[2]}\nof spectra', fontsize='xx-small', ha='right', va='top')
                # ax[x,y].set_title(f"Peak {peak[0]}", fontsize="xx-small")
            # creates a new figure if full
            y += 1
            y %= self.plot_output_columns
            if y == 0:
                x += 1
                x %= self.plot_output_rows
                if x == 0:
                    self.pdf.savefig(fig)
                    plt.close('all')
                    fig, ax = plt.subplots(self.plot_output_rows, self.plot_output_columns,figsize=(8.5,11))
        
        if x != 0:
            self.pdf.savefig(fig) # maybe get rid of fig parameter in here
        self.pdf.close()
        plt.close()

    def write_output(self):
        print("starting to write tsv file")
        # Writes out all observed peaks in a TSV file, including the measured mz value from the mass 
        # spectrometer, the intensity, and all the identifications
        # Each identification has the theoretical mass, the delta in mz, and the name of the theoretical ion

        output_peaks = self.observed_peaks.copy()

        if self.has_correction_spline:
            mz_values = [peak[0] for peak in output_peaks]
            x_values = np.array(mz_values)
            y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)
            if self.has_second_spline:
                y_values_2 = interpolate.BSpline(self.t2, self.c2, self.k2)(x_values)
                for index in range(len(y_values)):
                    y_values[index] += y_values_2[index]

        for index in range(len(output_peaks)):
            peak_mz = output_peaks[index][0]
            correction = (y_values[index] + self.crude_correction) * peak_mz / 1e6
            output_peaks[index].insert(1, round(peak_mz - correction, 5))
            for identification_index in range(len(output_peaks[index]) - 4):
                identification_index += 4
                output_peaks[index][identification_index][1] = round(output_peaks[index][identification_index][1] * 1e6 / (peak_mz - correction), 2)

        output_peaks.insert(0, ['uncorrected m/z', 'corrected m/z', 'intensity', 'percentage of spectra', 'primary identification', 'other identifications'])
        with open(f'{self.peak_file}.tsv', 'w') as file:
        # with open('common_peaks.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n')
            writer.writerows(output_peaks)

    def show_stats(self):
        # Prints out the stats and how long it took to run the file
        t1 = timeit.default_timer()
        print(f"INFO: Read {self.stats['counter']} spectra from {self.file_name}")
        print(f"The number of ms1spectra is {self.stats['ms1spectra']}")
        print(f"The number of ms2spectra is {self.stats['ms2spectra']}")

        print(f"INFO: Elapsed time: {t1-self.t0}")
        print(f"INFO: Processed {self.stats['counter']/(t1-self.t0)} spectra per second")
        return self.stats['ms2spectra']

    # Gaussian function used for curve fitting

def gaussian_function(x,a,x0,sigma):
    # a = amplitude
    # x0 = center point of the gaussian
    # sigma = measure of the width of the gaussian
    return a*exp(-(x-x0)**2/(2*sigma**2))

def get_strength(intensity, smallest_peak_intensity):
    # Returns the strength based on the smallest intensity and a given intensity value
    if intensity <= 3 * smallest_peak_intensity:
        return 1
    elif intensity <= 10 * smallest_peak_intensity:
        return 2
    elif intensity <= 30 * smallest_peak_intensity:
        return 3
    else:
        return 4

def process_job(job_parameters):
    peak_finder = MSRunPeakFinder(job_parameters["file"], job_parameters["tolerance"], job_parameters["rows"], job_parameters["columns"], job_parameters["make_pdf"], job_parameters["find_snippets"])
    print(f"working on {job_parameters['file']}")
    return peak_finder.process_file()

# put main at the end of the program, define identify peaks method first
def main():
    argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
    argparser.add_argument('--rows', action='store', default=5, type=int, help='Number of rows for graphs')
    argparser.add_argument('--columns', action='store', default=3, type=int, help='Number of columns for output')
    argparser.add_argument('--tolerance', action='store', default=5, type=float, help='Tolerance for identifying peaks in ppm')
    argparser.add_argument('--n_threads', action='store', type=int, help='Set the number of files to process in parallel (defaults to number of cores)')
    argparser.add_argument('--make_pdf', action='store', default=False, type=bool, help='Determines if the pdf should be created. Auto assumes False')
    argparser.add_argument('--find_snippets', action='store', default=False, type=bool, help='Determines if the program should search for IH snippets. Auto assumes False')
    argparser.add_argument('files', type=str, nargs='+', help='Filenames of one or more mzML files to read')

    params = argparser.parse_args()

    for file in params.files:
        file_stats = os.stat(file)
        if file_stats.st_size == 0:
            print(f"ERROR: File '{file}' has a size of 0")
            return
        if not os.path.isfile(file): # add something that can glob the asterisks
            print(f"ERROR: File '{file}' not found or not a file")
            return
        if file[-4:] != "mzML" and file[-7:] != "mzML.gz":
            print(f"ERROR: File '{file}' is not an mzML nor a mzML.gz")
            return
        
    start_time = timeit.default_timer()
    total_ms2_spectra = 0
    
    jobs = []
    for file in params.files:
        jobs.append({"file": file, "tolerance": params.tolerance, "rows": params.rows, "columns": params.columns, "make_pdf": params.make_pdf, "find_snippets": params.find_snippets})

    if len(jobs) == 1:
        process_job(jobs[0])
    else:
        #### Process the jobs in parallel
        n_threads = params.n_threads or multiprocessing.cpu_count()
        print(f"Processing files with n_threads={n_threads} (one mzML per thread)", end='', flush=True)
        pool = multiprocessing.Pool(processes=n_threads)
        # results = pool.map_async(process_job, jobs)
        results = pool.map(process_job, jobs)
        print(results)
        #results = pool.map(process_job, jobs)
        pool.close()
        pool.join()
        print("")

    print(f"Total time to run {len(jobs)} files is {timeit.default_timer() - start_time}")

if __name__ == "__main__": main()
