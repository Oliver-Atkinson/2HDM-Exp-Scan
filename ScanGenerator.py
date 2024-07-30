"""Generating a given number of random points in the tanb-mhp space to be checked against LHC data"""
from Functions import *

#The number of points and corresponding files to generate
npoints = int(sys.argv[1])
typ = int(sys.argv[2])
m_config = sys.argv[3]
nfiles = int(npoints/20)+1
if npoints % 20 == 0:
    nfiles -= 1

#Directory to store the possibly large number of files
scandir = "ScanPointSets"    
    
#The tanb-mass space to investigate (in log space)
#tb_min, tb_max = -1.5, 2.5
#m_min, m_max = 1, 4
tb_min, tb_max = np.log10(1.8), np.log10(2.1)
m_min, m_max = np.log10(88), np.log10(104)
tbran_l, tbran_h = int(1000*npoints*tb_min), int(1000*npoints*tb_max)
mran_l, mran_h = int(1000*npoints*m_min), int(1000*npoints*m_max)

tbrans = random.sample(range(tbran_l, tbran_h), npoints)
tbs = [ran/(1000*npoints) for ran in tbrans]

mrans = random.sample(range(mran_l, mran_h), npoints)
ms = [ran/(1000*npoints) for ran in mrans]

#For checking purposes and BR generation
#tbs = list(np.zeros(npoints))
#ms = np.linspace(1, 2.5, npoints)
#tbs = np.linspace(-1.5, 2.5, npoints)
#ms = list(np.zeros(npoints)+3)

#d1points = int(np.sqrt(npoints))
#tbs0 = np.linspace(np.log10(1.8), np.log10(2.1), d1points)
#ms0 = np.linspace(np.log10(88), np.log10(104), d1points)
#tbs=np.repeat(tbs0,d1points)
#ms = [ms0[i%d1points] for i in range(npoints)]

#Fixed mass values (GeV)
m1 = 5000
m2 = 500


#Generating a file of 20 points
def FileGen(typ, nfile, nps, config=m_config):
    #Output file
    
    outfile_name = scandir+"/Set"+str(nfile+1)+".dat"
    outfile = open(outfile_name, 'w') 
    outfile.write("type tanb mHp mH0 mA0 cos(b-a) \n")
    
    #Setting up for the configuration as a variable input
    ms_in = [m1,m1,m1]
    m_dict = {
        "H+": [0],
        "H0": [1],
        "A0": [2],
        "Degen": [0,1,2]}
    
    #Validity check, maybe not required
    #if config not in m_dict:
    #    print("Invalid mass configuration input (argument 3)")
    #    exit()
        
    #Getting tanb and mhp
    for i in range(nps):
        tanb, m = 10**tbs[(nfile*20)+i], 10**ms[(nfile*20)+i]
        for mID in m_dict[config]:
            ms_in[mID] = m
               
        #Assuming exact alignment limit
        outfile.write(str(typ)+" "+"{0:.4f}".format(tanb)+" "+"{0:.2f}".format(ms_in[0])+" "+"{0:.2f}".format(ms_in[1])+" "+"{0:.2f}".format(ms_in[2])+" "+"0.0\n")
    
    outfile.close()
    
    return

#Writing the files
for i in range(nfiles):
    if npoints >= 20:
        FileGen(typ, i, 20)
    else:
        FileGen(typ, i, npoints)
    npoints -= 20