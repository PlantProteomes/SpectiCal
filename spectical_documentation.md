# SpectiCal Capabilities

# Running SpectiCal and Commands
In order to run the spectical.py script, it is required to provide an input mzML or mzML.gz file.<br>
The command to run the script is as follows:<br>
```python spectical.py filename.mzML```<br>
Multiple files can be inputted at once as well. Example:<br>
```python spectical.py *.mzML```<br>
Additional commands can be added, such as:<br>
```--rows n -> the number of rows in the pdf output. By default, n = 5<br>
--columns n -> the number of columns in the pdf output. By default, n = 3<br>
--tolerance n -> the tolerance for identifying the peak's ppm. By defualt, n = 5<br>
--n_threads n -> the number of threads you want to run the program with. By default, n = number of cores<br>
--make_pdf n -> toggles whether the pdf is generated. By default, n = false<br>
--find_snippets n -> toggles whether the snippets of IH, IF, and IK+CO are collected. By default, n = false<br>```
<be>

# Combining TSV outputs from SpectiCal
In order to run the combine_list.py script, it is required to include all the tsv files that you want to be combined.<br>
The command to run the script is as follows:<br>
python combine_list.py file_one.tsv, file_two.tsv<br>
All tsv files in a folder can also be inputted at once. The command is as follows:
python combine_list.py *.tsv
You may also choose to filter the TSV file so only the most commonly found ions are included. By default, it will not filter the file<br>
The command is as follows:<br>
python combine_list.py --filter_tsv True file_one.tsv, file_two.tsv<br>
<be>

# Correcting input files using calibration generated from SpectiCal
To run the shift_mzML.py, the input file to be corrected, the name of the new file, and the correction constants are required.<br>
There are two ways to do this:<br>
1. provide a constant PPM shift<br>
The script is as follows:<br>
python shift_mzML.py --input_filename input.mzML --output_filename input_calibrated.mzML --ppm_shift 5<br>
2. provide the JSON file generated from spectical.py<br>
An example script is as follows:<br>
python shift_mzML.py --input_filename input.mzML --output_filename input_calibrated.mzML --json_filename input.calibration.json<br>
