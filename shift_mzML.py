#!/usr/bin/env python3

import sys
import os
import argparse
import os.path
import json
import gzip
import numpy as np

from pyteomics import mzml, auxiliary
from scipy import interpolate
from psims.transform.mzml import MzMLTransformer

def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


class MyMzMLTransformer(MzMLTransformer):
    def __init__(self, input_stream, output_stream, transform=None, transform_description=None,
                 sort_by_scan_time=False, ppm_shift=None, alter_ms_levels=None, correction_values=None):
        super().__init__(input_stream, output_stream, transform=transform, transform_description=transform_description,
                 sort_by_scan_time=sort_by_scan_time)

        self.has_spline_correction = False
        self.has_second_spline = False
        self.correction_cache = {}

        self.alter_ms_levels = None
        if alter_ms_levels is not None:
            self.alter_ms_levels = alter_ms_levels.split(',')
            for i in range(len(self.alter_ms_levels)):
                self.alter_ms_levels[i] = int(self.alter_ms_levels[i])

        if ppm_shift is not None:
            self.ppm_shift = float(ppm_shift)
        else:
            self.ppm_shift = float(correction_values['crude correction'])
            self.after_400_calibration = correction_values['after_400_calibration']
            if len(correction_values) >= 5:
                self.has_spline_correction = True
                self.t = correction_values['t1']
                self.c = correction_values['c1']
                self.k = int(correction_values['k1'])
                if len(correction_values) >= 8:
                    self.has_second_spline = True
                    self.t2 = correction_values['t2']
                    self.c2 = correction_values['c2']
                    self.k2 = int(correction_values['k2'])

    def format_spectrum(self, spectrum):
        new_spectrum = super().format_spectrum(spectrum)

        #### Shift all the m/z values
        vfunc = np.vectorize(self.correct_mz)
        new_spectrum['mz_array'] = vfunc(new_spectrum['mz_array'])

        #### Get the MS level
        ms_level = None
        for param in new_spectrum['params']:
            if param['name'] == 'MS:1000511':
                ms_level = param['value']

        #### Correct the precursor m/z values by the requested shift
        if ms_level is not None and ms_level > 1:
            precursor_mz = new_spectrum['precursor_information'][0]['mz']
            new_spectrum['precursor_information'][0]['mz'] = self.correct_mz(precursor_mz)

        return new_spectrum

    def correct_mz(self, input_mz):
        int_input_mz = int(input_mz)
        if int_input_mz in self.correction_cache:
            return input_mz + self.correction_cache[int_input_mz]
        else:
            total_correction_offset = -1 * self.ppm_shift * input_mz / 1e6
            if self.has_spline_correction:
                if input_mz <= 400:
                    x_values = np.array([input_mz])
                    y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)
                    if self.has_second_spline:
                        y_values_2 = interpolate.BSpline(self.t2, self.c2, self.k2)(x_values)
                        for index in range(len(y_values)):
                            y_values[index] += y_values_2[index]
                    total_correction_offset -= y_values[0] * input_mz / 1e6
                else:
                    total_correction_offset -= self.after_400_calibration * input_mz / 1e6
            self.correction_cache[int_input_mz] = total_correction_offset
            return input_mz + total_correction_offset
                    

def main():

    argparser = argparse.ArgumentParser(description='Read an input mzML file and write out a new one with the m/zs shifted')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--input_filename', type=str, action='store', required=True, help='Name of the input mzML file')
    argparser.add_argument('--output_filename', type=str, action='store', required=True, help='Name of the output mzML file shifted m/zs')
    argparser.add_argument('--ppm_shift', type=str, action='store', help='Offset to shift all spectra in units of PPM')
    argparser.add_argument('--mslevels', type=str, action='store', help='MS Levels to shift, comma separated (e.g. 1 or 2 or 2,3)')
    argparser.add_argument('--json_filename', type=str, action='store', help='Name of the json file with calibration values')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    #### Check that the filenames are defined and exist
    if not os.path.isfile(params.input_filename):
        print(f"ERROR: --input_filename '{params.input_filename}' not found or not a file")
        return

    #### Ensure the output is not the same as the input
    if params.output_filename == params.input_filename:
        print(f"ERROR: --output_filename '{params.input_filename}' may not be the same as the --input_filename")
        return

    #### Open in the input file
    print(f"INFO: Opening {params.input_filename} for reading")
    if params.input_filename.endswith('.gz'):
        infile = gzip.open(params.input_filename)
    else:
        infile = open(params.input_filename, 'rb')

    #### Open the output file
    print(f"INFO: Opening {params.output_filename} for writing")
    try:
        outfile = open(params.output_filename, 'wb')
    except:
        print(f"ERROR: --output_filename '{params.output_filename}' is not writable. Permissions problem?")
        return
    
    selected_spectra = 'all'
    if params.mslevels is not None:
        selected_spectra = params.mslevels

    # Determines if json file is used or ppm shift value
    if params.json_filename is None or params.json_filename == "":
        if params.ppm_shift is None:
            print("ERROR: --json_filename and --ppm_shift were not found")
            return
        MyMzMLTransformer(infile, outfile, ppm_shift=params.ppm_shift, alter_ms_levels=params.mslevels,
            transform_description=f"Shifted {selected_spectra} spectra by {params.ppm_shift}").write()
    else:
        if not os.path.isfile(params.json_filename):
            print(f"ERROR: --json_filename '{params.json_filename}' not found or not a file")
            return
            
        correction_file = open(f'{params.json_filename}')
        correction_values = json.load(correction_file)
        MyMzMLTransformer(infile, outfile, correction_values=correction_values, alter_ms_levels=params.mslevels,
            transform_description=f"Shifted {selected_spectra} spectra by values {correction_values['crude correction']} and spline fit").write()

#### For command line usage
if __name__ == "__main__": main()
