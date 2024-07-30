"""File to convert a given input of points (tanb, mH+, mH0, mA0, cos(b-a)) to the format required to input to BSMPT"""
import numpy as np
import sys

#Constants
mh = 125.25
v = 246.21965079413735

#Reading in masses and tanbs
readfile_name = sys.argv[1]
npoints = int(sys.argv[2])
m_basisDat = np.loadtxt(readfile_name, skiprows=1, delimiter=' ')

#The output file here is to be used as input in BSMPT
outfile_name = "BSMPT_Input.dat"
outfile = open(outfile_name,"w") 
outfile.write("type	L1	L2	L3	L4	L5	m12sq	tbeta	mH	mh	mA	mHc	alpha")

#Calculating the additional parameters for BSMPT
def params_calc(point):
    typ, tanb, mhp, mH0, mA0, cosba = point[0], point[1], point[2], point[3], point[4], point[5]
    
    b = np.arctan(tanb)
    ba = np.arccos(cosba)
    a = b-ba
    
    m_12_2 = (mH0**2)*np.sin(b)*np.cos(b) #Not always valid
    
    #Eq. 13 in 2003.06170 - lambda_1,2 = 0 in the approx that m_12^2 = mH0^2*sinb*cosb - Requires updating
    lam1n2 = (mh/v)**2
    lam3 = (mh**2 + 2*(mhp**2) - 2*(mH0**2))/ (v**2) 
    lam4 = (mA0**2 - 2*(mhp**2) + mH0**2) / (v**2)
    lam5 = (-mA0**2 + mH0**2) / (v**2)
    
    return lam1n2, lam1n2, lam3, lam4, lam5, m_12_2, a

#Calculating the parameters and writing them to the BSMPT input file (note the ordering, particularly for the masses)
for i in range(npoints):
    params = params_calc(m_basisDat[i])
    outfile.write("\n"+str(m_basisDat[i,0])+" "+str(params[0])+" "+str(params[1])+" "+str(params[2])+" "+str(params[3])+" "+str(params[4])+" "+str(params[5])+" "+str(m_basisDat[i,1])+" "+str(m_basisDat[i,3])+" "+str(mh)+" "+str(m_basisDat[i,4])+" "+str(m_basisDat[i,2])+" "+str(params[6]))

