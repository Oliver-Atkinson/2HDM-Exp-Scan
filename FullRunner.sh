#!/bin/bash
echo -e "Generating the data, plots and scans for the Type 1 and Type 2 2HDMs with $1 points per scan\n"
echo -e "Warning: Ensure the relevant functions in the HB_FileWriter.py and ScanResults.py are switched on if BR plots are desired\n"
bash Scanner.sh $1 1 Degen 0
#bash Scanner.sh $1 1 H+ 0
#bash Scanner.sh $1 1 H0 0
#bash Scanner.sh $1 1 A0 0
bash Scanner.sh $1 2 Degen 0
#bash Scanner.sh $1 2 H+ 0
#bash Scanner.sh $1 2 H0 0
#bash Scanner.sh $1 2 A0 0