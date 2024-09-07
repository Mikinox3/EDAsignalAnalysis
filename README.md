# EDAsignalAnalysis
Python script for EDA signal preprocessing and analysis (data acquisition with Biopac device and ACQKNOWLEDGE software).
Includes Neurokit2 functions for filtering, artifact detection and parameters calculation. 
Neurokit function eda_preprocess extracts from the raw signal the SCL and SCR, detects the peaks in the signal and returns their amplitude, latency, peak onset, peak rise time
The spectral power of the signal and AUC are also calculated.

See https://neuropsychology.github.io/NeuroKit/ for more info.

Please note : This script has been created for a very specific experimental paradigm with precise event triggers during EDA recording. It cannot be used as it is currently written and has to be adapted for your own experimental paradigm (especially for the Cut_the_file.py).
