#!/bin/bash
echo -e "Performing a scan of the tanb-m$3 space in the Type $2 2HDM with $1 random points at an integrated luminosity of $4fb^-1, using 2HDecay and HiggsBounds to check if they are allowed by LHC data"
echo "Warning: must be in an eviornment with cmake installed, such as hep_env, for a new extrapolation to be run"


echo -e "\nCleaning files (likely that some don't yet exist)..."
rm -r ScanPointSets
mkdir ScanPointSets
rm InterpolationData/ScanData/*

echo "Generating the random points and storing in files of 20 points"
python ScanGenerator.py $1 $2 $3
let n=($1/20)

#Extrapolation of the data
cd /home/oliver/HEPTools/
if [ -d  "/home/oliver/HEPTools/higgsbounds-5.10.1_$4" ]; then
     echo "A HiggsBounds file with this luminosity already exists, this will be used"
     mv higgsbounds-5.10.1_$4 higgsbounds-5.10.1
else        
    echo "Unpacking the fresh HiggsBounds code..."
    tar -xzf higgsbounds-5.10.1.tar.gz
    cd higgsbounds-5.10.1

    echo "Performing the extrapolation of the data..."
    cd /home/oliver/2HDMII_NewP/
    python Extrapolator.py $4
    cp -a ExtrapolationData/ATLAS_$4/. /home/oliver/HEPTools/higgsbounds-5.10.1/data/Expt_tables/ATLtables
    cp -a ExtrapolationData/CMS_$4/. /home/oliver/HEPTools/higgsbounds-5.10.1/data/Expt_tables/CMStables
    cd /home/oliver/HEPTools/higgsbounds-5.10.1/

    echo "Building the new HiggsBounds with extrapolated data..."
    mkdir build && cd build
    cmake .. > /dev/null
    make 2> /dev/null
fi


echo "Calculating cross sections..."
cd /home/oliver/2HDMII_NewP/
python Interpolater.py $1 1

echo -e "Beginning scan...\n"
SECONDS=0
INTERVAL=0
for (( i=1; i<=$n; i++ ))
do
    STARTTIME=$(date +%s)
    python 2HDecay_FileWriter.py ScanPointSets/Set$i.dat 20
    cp -a 2HDecay_Inputs/. /home/oliver/HEPTools/2HDECAY-master/Input/
    cd /home/oliver/HEPTools/2HDECAY-master/
    python 2HDECAY.py > /dev/null 2> /dev/null
    cp -a Results/. /home/oliver/2HDMII_NewP/2HDecay_Outs/
    cd /home/oliver/2HDMII_NewP/

    python HB_FileWriter.py ScanPointSets/Set$i.dat 20
    cp -a HB_Inputs/. /home/oliver/HEPTools/higgsbounds-5.10.1/build/Inputs/
    cd /home/oliver/HEPTools/higgsbounds-5.10.1/build/
    for (( j=1; j<=20; j++ ))
    do  
       ./HiggsBounds LandH SLHA 3 1 Inputs/Point$j.slha > /dev/null
    done
    cp -a Inputs/. /home/oliver/2HDMII_NewP/HB_Outs/
    cd /home/oliver/2HDMII_NewP/
    
    python ResultExtractor.py ScanPointSets/Set$i.dat 20 n
    ENDTIME=$(date +%s)
    INTERVAL=$(($ENDTIME - $STARTTIME))
    REMAIN=$((($n - $i)*$SECONDS/$i))
    echo -e "\e[1ASet $i of $n done in $INTERVAL s, with a total elapsed time of $SECONDS s and an estimated finish at in $REMAIN seconds"
done

echo -e "\nPlotting the results.."
python ScanResults.py $1 $2 $3
echo "All computations completed, cleaning files..."
mv /home/oliver/HEPTools/higgsbounds-5.10.1 /home/oliver/HEPTools/higgsbounds-5.10.1_$4
rm -r ScanPointSets
rm InterpolationData/ScanData/*