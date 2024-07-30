"""Extrapolates the LHC data used in HiggsBounds to a given luminosity"""

from Functions import *
import os

#The luminosity - in fb^-1
Lum = float(sys.argv[1])

#File names
HB_dir = "/home/oliver/HEPTools/higgsbounds-5.10.1/"
Atlas_dir = HB_dir+"data/Expt_tables/ATLtables/"
CMS_dir = HB_dir+"data/Expt_tables/CMStables/"


#Getting the filenames with the data in
ATLfiles = [f for f in os.listdir(Atlas_dir) if os.path.isfile(os.path.join(Atlas_dir, f))]
CMSfiles = [f for f in os.listdir(CMS_dir) if os.path.isfile(os.path.join(CMS_dir, f))]
    
#Extrapolates the values in each file to a new luminosity
def Extrapolator(filename, expdir, newLum):
    #The extrpolation scale factor
    ogLum = LumGetter(filename)
    scale = np.sqrt(ogLum/newLum)
    
    #Getting the data
    readname = expdir+filename
    with open(readname, 'r') as file:
        content = file.readlines()
    dat = content[5:] #Line 6 onwards for the data values
    
    #New data to write over the file, keeping the header (first 5 lines) the same
    newDat = []
    for i in range(5):
        newDat.append(content[i])
    
    #Going line by line to do the extrapolation
    for line in dat:
        if len(line.strip()) != 0: #Only doing work for non-empty lines
            comps = line.strip().split()
            m, obs, exp = comps[0], float(comps[1]), float(comps[2])
            obs_ext, exp_ext = exp*scale, exp*scale #obs and exp should be the same in extrapolations
            newLine = "         "+m+"   "+'{:.7g}'.format(obs_ext)+"       "+'{:.7g}'.format(exp_ext)+'    \n'
            newDat.append(newLine)
    
    #Writing the new data to file
    with open(readname, 'w') as file:
        for line in newDat:
            file.write(line)
        
    return

#Running the extrapolation
print("Beginning extrapolation...")
for filename in ATLfiles:
    Extrapolator(filename, Atlas_dir, Lum)
for filename in CMSfiles:
    Extrapolator(filename, CMS_dir, Lum)