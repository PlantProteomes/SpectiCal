# SpectiCal

In order to run the spectical.py script, it is required to provide an input mzML or gz file.
The command to run the script is as follows:
python spectical.py filename.mzML
Multiple files can be inputted at once as well. Example:
python spectical.py *.mzML
Additional commands can be added, such as:
--rows n -> the number of rows in the pdf output. By default, n = 5
--columns n -> the number of columns in the pdf output. By default, n = 3
--tolerance n -> the tolerance for identifying the peak's ppm. By defualt, n = 5
--n_threads n -> the number of threads you want to run the program with. By default, n = number of cores
--make_pdf n -> toggles whether the pdf is generated. By default, n = false
--find_snippets n -> toggles whether the snippets of IH, IF, and IK+CO are collected. By default, n = false

In order to run the combine_list.py script, it is required to include all the tsv files to be combined.
The command to run the script is as follows:
python combine_list.py file_one.tsv, file_two.tsv
You may also choose to filter the TSV file so only the most commonly found ions are included. By default, it will filter the file
The command is as follows:
python combine_list.py --filter_tsv True file_one.tsv, file_two.tsv

To run the shift_mzML.py, the input file to be corrected, the name of the new file, and the correction constants are required.
There are two ways to do this:
1. provide a constant PPM shift
The script is as follows:
python shift_mzML.py --input_filename input.mzML --output_filename input_calibrated.mzML --ppm_shift 5
2. provide the JSON file generated from spectical.py
An example script is as follows:
python shift_mzML.py --input_filename input.mzML --output_filename input_calibrated.mzML --json_filename input.calibration.json
