# SpectiCal

SpectiCal is a script that computes an m/z calibration based on known low-mass ions, then prints out a PDF displaying the low-mass ions and common peaks in the run.

To run SpectiCal, follow the steps:
1. Clone the repository - "git clone https://github.com/PlantProteomes/SpectiCal/tree/main"
2. Install required packages - "pip install -r requirements.txt"
3. Download the data files to run. Here is a sample 10,000 spectra subset data file of a full run (28 MB) - "curl -O https://peptideatlas.org/refdata/HFX_9850_GVA_DLD1_2_180719_subset.mzML.gz"
4. Run the data file - "python spectical.py HFX_9850_GVA_DLD1_2_180719_subset.mzML.gz"

Once the data file has been run, you may open the TSV, PDF, and JSON files to view the results of the script.
For more information about additional metrics and capabilities, please refer to spectical_documentation.md
