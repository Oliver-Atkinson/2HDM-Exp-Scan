"""Uses an interpolation function for H+ production cross sections from LHC data stored in Higgs Bounds to generate cross sections for all the scan points"""

from Functions import *

npoints = int(sys.argv[1])
scan = int(sys.argv[2]) #Switch for scan vs analysis - value of 1 means scan
nfiles = int(npoints/20) 

#Directories for the various file sets
DatDir = "InterpolationData/"
ScanDatDir = "InterpolationData/ScanData/"
ReadDir = "ScanPointSets/"

#Creating the interpolation functions from a specified for each channel
t_b_obs_calc = IntMaker(DatDir+"2020039_ATLAS_Hp_tb_139fb-1.txt", 5, 0, 2) 
t_b_mg_calc = IntMaker(DatDir+"Mg_tbhxsScan.dat", 1, 0, 1)

#Function to write out the interpolated result for a given input file
def IntWriter(readfile, nfile, mg_func, obs_func):
    #Reading in
    m_basisDat = np.loadtxt(readfile, skiprows=1, delimiter=' ')
    
    #The outfile
    outfile_name = ScanDatDir+"Set"+str(nfile+1)+"_xs.dat"
    outfile = open(outfile_name, 'w') 
    outfile.write("mHp tanb(beta) xs_tbHpm_mg (xs_tbHpm)(BR_H+_tb)_obs\n")
    
    mt, mb, v = 173, 4.1, 246
    #Getting and writing the cross section
    for line in m_basisDat:
        typ, tanb, mhp = line[0], line[1], line[2]
        if typ == 1:
            ep_d = 1/tanb
        elif typ == 2:
            ep_d = tanb
        xsBR = 10**obs_func(mhp)
        xs_tb0 = 10**mg_func(mhp) #Conversion from log for observational xs
        Yuk = 2*((mt/tanb)**2+(mb*ep_d)**2)/(v**2)
        xs_tb_Yuk = xs_tb0*Yuk #Accounting for the Yukawa factor
        outfile.write('{:.6g}'.format(mhp)+" "+'{:.3g}'.format(tanb)+" "+'{:.8e}'.format(xs_tb_Yuk)+" "+'{:.8e}'.format(xsBR)+"\n")
    
    outfile.close()
    
    return

#Calculating the xs for each file and writing out to a file, depending on if a scan is to be performed
if scan == 1:
    for i in range(nfiles):
        IntWriter("ScanPointSets/Set"+str(i+1)+".dat", i, t_b_mg_calc, t_b_obs_calc)
else:
    readfile_name = sys.argv[3]
    IntWriter(readfile_name, 0, t_b_mg_calc, t_b_obs_calc)

