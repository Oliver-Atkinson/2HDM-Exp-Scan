#!/bin/bash

echo "Check for EWPT? (y/n)"
read EWPT

if [ "$EWPT" = "y" ]
then
    echo -e "\e[1AUsing 2HDecay and HiggsBounds to check if the first $2 points in file $1 are allowed by LHC data, extrapolated to a luminosity of $3fb^-1, and BSMPT to check if these points give a SFOEWPT"
else        
    echo -e "\e[1AUsing 2HDecay and HiggsBounds to check if the first $2 points in file $1 are allowed by LHC data, extrapolated to a luminosity of $3fb^-1"
fi

echo "Cleaning files..."
rm 2HDecay_Inputs/*
rm 2HDecay_Outs/*
rm /home/oliver/HEPTools/2HDECAY-master/Input/*
rm /home/oliver/HEPTools/2HDECAY-master/Results/*
rm HB_Inputs/*
rm HB_Outs/*
rm /home/oliver/HEPTools/higgsbounds-5.10.1/build/Inputs/*
if [ "$EWPT" = "y" ]
then
    rm BSMPT_Output.dat
fi

#Extrapolation of the data
cd /home/oliver/HEPTools/
if [ -d  "/home/oliver/HEPTools/higgsbounds-5.10.1_$3" ]; then
     echo "A HiggsBounds file with this luminosity already exists, this will be used"
     mv higgsbounds-5.10.1_$3 higgsbounds-5.10.1
else        
    echo "Unpacking the fresh HiggsBounds code..."
    tar -xzf higgsbounds-5.10.1.tar.gz
    cd higgsbounds-5.10.1

    echo "Performing the extrapolation of the data..."
    cd /home/oliver/2HDMII_NewP/
    python Extrapolator.py $3
    cp -a ExtrapolationData/ATLAS_$3/. /home/oliver/HEPTools/higgsbounds-5.10.1/data/Expt_tables/ATLtables
    cp -a ExtrapolationData/CMS_$3/. /home/oliver/HEPTools/higgsbounds-5.10.1/data/Expt_tables/CMStables
    cd /home/oliver/HEPTools/higgsbounds-5.10.1/

    echo "Building the new HiggsBounds with extrapolated data..."
    mkdir build && cd build
    cmake .. > /dev/null
    make 2> /dev/null
fi



echo "Calculating cross sections..."
cd /home/oliver/2HDMII_NewP/
python Interpolater.py $2 0 $1 

echo "Writing new input files for 2HDecay..."
python 2HDecay_FileWriter.py $1 $2  
cp -a 2HDecay_Inputs/. /home/oliver/HEPTools/2HDECAY-master/Input/
cd /home/oliver/HEPTools/2HDECAY-master/
echo "Running 2HDecay..."
python 2HDECAY.py > /dev/null 2> /dev/null
echo "Moving 2HDecay results..."
cp -a Results/. /home/oliver/2HDMII_NewP/2HDecay_Outs/
cd /home/oliver/2HDMII_NewP/

echo "Extracting BRs, calculating couplings and creating HiggsBounds input files..."
python HB_FileWriter.py $1 $2
cp -a HB_Inputs/. /home/oliver/HEPTools/higgsbounds-5.10.1/build/Inputs/
cd /home/oliver/HEPTools/higgsbounds-5.10.1/build/
echo "Running HiggsBounds..."
for (( i=1; i<=$2; i++ ))
do  
   ./HiggsBounds LandH SLHA 3 1 Inputs/Point$i.slha  > /dev/null 2> /dev/null
done
echo "Moving HiggsBounds results..."
cp -a Inputs/. /home/oliver/2HDMII_NewP/HB_Outs/
cd /home/oliver/2HDMII_NewP/

if [ "$EWPT" = "y" ]
then
    echo "Beginning the BSMPT analysis"
    python BSMPT_Converter.py $1 $2
    mv BSMPT_Input.dat //home/oliver/HEPTools/BSMPT/build/2HDM_Input.dat
    echo "Conversion to lambda basis completed"
    cd  //home/oliver/HEPTools/BSMPT/build
    echo "Running BSMPT (may take a little time)..."
    let n=$2+1
    ./bin/BSMPT r2hdm 2HDM_Input.dat 2HDM_Output.dat 2 $n
    mv 2HDM_Output.dat //home/oliver/2HDMII_NewP/BSMPT_Output.dat
    echo "Calculations complete - see BSMPT_Output.dat"
    cd //home/oliver/2HDMII_NewP
fi

echo "Extracting the top-level results..."
python ResultExtractor.py $1 $2 $EWPT

if [ "$EWPT" = "y" ]
then
    cd  //home/oliver/HEPTools/BSMPT/build
    echo "Cleaning BSMPT - this will take a few minutes (minor warnings expected)..."
    make clean
    make #> /dev/null
    cd //home/oliver/2HDMII_NewP
fi


mv /home/oliver/HEPTools/higgsbounds-5.10.1 /home/oliver/HEPTools/higgsbounds-5.10.1_$3
echo "All computations completed"