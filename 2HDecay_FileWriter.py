"""Takes the basic input file, does simple calculations (m12^2 and alpha), then writes out an input file for use in 2HDecay calculations, using Example2HDecay.in as a template"""

from Functions import *

#Reading in masses and tanbs and cos(b-a)s
readfile_name = sys.argv[1]
npoints = int(sys.argv[2])
m_basisDat = np.loadtxt(readfile_name, skiprows=1, delimiter=' ')

#Creating a new file from Example2HDecay.in and the paramters of a point
def HDecay_FileWriter(npoint, pointdat):
    typ, tanb, m12sq, alpha, mH0, mA, mHp = pointdat
    
    #The output file
    outfile = "2HDecay_Inputs/In"+str(npoint)+".in"
    
    #Copying the template file
    shutil.copy2("Example2HDecay.in", outfile)

    #Opening the template file, then changing the relevant data and writing out to the new file 
    with open('Example2HDecay.in', 'r') as file:
        data = file.readlines()
    
    data[56] = "TYPE     = "+"{0:.1g}".format(typ)+"\n"
    data[60] = "TGBET2HDM= "+"{0:.2f}".format(tanb)+"\n"
    data[61] = "M_12^2   = "+"{0:.3f}".format(m12sq)+"D0"+"\n"
    data[65] = "ALPHA_H  = "+"{0:.12g}".format(alpha)+"D0"+"\n"
    data[67] = "MHH      = "+"{0:.2f}".format(mH0)+"D0"+"\n"
    data[68] = "MHA      = "+"{0:.2f}".format(mA)+"D0"+"\n"
    data[69] = "MH+-     = "+"{0:.2f}".format(mHp)+"D0"+"\n"

    with open(outfile, 'w') as file:
        file.writelines(data)
        
    return

#For each point in the input, calculate the parameters then write out to a new file    
for i in range(npoints):
    params = HDecayCalc(m_basisDat[i])
    HDecay_FileWriter(i+1, params)

