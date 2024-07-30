"""Takes the 2HDecay output and uses it, along with some further calculations, to write out HiggsBounds input files, based on the template ExampleHBInput.slha"""

from Functions import *

#Reading in masses and tanbs and cos(b-a)s, as well as the cross sections
readfile_name = sys.argv[1]
npoints = int(sys.argv[2])
m_basisDat = np.loadtxt(readfile_name, skiprows=1, delimiter=' ')

#Reading in in a different way for the scan and single file analysis
if "ScanPointSets" in readfile_name:
    n_end = readfile_name.find(".")
    n = readfile_name[17:n_end]
    xs_file = "InterpolationData/ScanData/Set"+str(n)+"_xs.dat"
else:
    xs_file = "InterpolationData/ScanData/Set1_xs.dat"
xsDat = np.loadtxt(xs_file, skiprows=1, delimiter=' ')

#xs_checker(xsDat, readfile_name) #Can remove this once the check is complete
#for p in ["h0", "H0", "H+", "A0"]: #Can remove this once BRs are plotted
#    BR_DatGetter(readfile_name, npoints, p)

#Creating a new file from ExampleHBInput.slha and the paramters of a point
def HB_FileWriter(npoint):
    #Getting the data
    typ, tanb, mHp, mH0, mA0, cosba = m_basisDat[npoint]
    BR_Blocks, GamList = HDecay_InfoGetter(npoint)
    m, tb, xs_tb, xs_BR = xsDat[npoint]
    
    #Slight redundancy as we calculate this elsewhere but it's not too bad really
    b = np.arctan(tanb)
    ba = np.arccos(cosba)
    alpha = b-ba
    
    #Calculating the couplings
    Kappas = CouplingCalc(typ, tanb, cosba)

    #The output file
    outfile = "HB_Inputs/Point"+str(npoint+1)+".slha"
    
    #Copying the template file
    shutil.copy2("ExampleHBInput.slha", outfile)

    #Opening the template file, then changing the relevant data and writing out to the new file 
    with open('ExampleHBInput.slha', 'r') as file:
        NewDat = file.readlines()
        

    #Basic parameters - lines are fixed for these
    NewDat[26] = "         3     "+'{:.8E}'.format(tanb)+"   # TB\n"
    NewDat[30] = "        35     "+'{:.8E}'.format(mH0)+"   # MHH\n"
    NewDat[31] = "        36     "+'{:.8E}'.format(mA0)+"   # MA0\n"
    NewDat[32] = "        37     "+'{:.8E}'.format(mHp)+"   # MHp\n"
    NewDat[34] = "              "+'{:.8E}'.format(alpha)+"   # Alpha\n"
    
    #Cross sections - so far (and possibly ever) only t-b-Hpm
    NewDat[66] = "5 6 37 "+'{:.8E}'.format(xs_tb)+" # t-b-Hpm production\n"
    
    #Maybe get a function/itertor for these...
    #Decay Widths and BRs - always starts at line 70 but then things get a little nasty...
    #Setting the 125 GeV scalar to be exactly that of the SM, using numbers from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
    
    #If all the new Higgses are heavy (>2 TeV), set h^0 to match SM behaviour (2HDecay struggles in extreme tanb regions otherwise)
    if mHp < 2e3 and mH0 < 2e3 and mA0 < 2e3:
        s = 69 #Important dyanmic line counter 
        NewDat[s] = "DECAY        25     "+GamList[0]+"   # Gamma(h0)\n"
        s += 1
        f = s+len(BR_Blocks[0])
        NewDat[s:f] = BR_Blocks[0]
        s = f
    else:
        s = 80
    
    
    NewDat[s] = "DECAY        35     "+GamList[1]+"   # Gamma(H0)\n"
    s += 1
    f = s+len(BR_Blocks[1])
    NewDat[s:f] = BR_Blocks[1]
    
    s = f
    NewDat[s] = "DECAY        36     "+GamList[2]+"   # Gamma(A0)\n"
    s += 1
    f = s+len(BR_Blocks[2])
    NewDat[s:f] = BR_Blocks[2]
    
    s = f
    NewDat[s] = "DECAY        37     "+GamList[3]+"   # Gamma(Hp)\n"
    s += 1
    f = s+len(BR_Blocks[3])
    NewDat[s:f] = BR_Blocks[3]
   
    s = f
    NewDat[s] = "DECAY        6     "+GamList[4]+"   # Gamma(top)\n"
    s += 1
    f = s+len(BR_Blocks[4])
    NewDat[s:f] = BR_Blocks[4]
        
    #Getting the back end of the original data to restore to the file
    with open('ExampleHBInput.slha', 'r') as file:
        datBack = file.readlines()[118:149]
    
    NewDat[f:f+31] = datBack

    
    #Now to write out the Boson couplings, assuming that the Hgg effective coupling is purely from the top quark triangle loop (not necessarily valid)
    s=f+3

    for i in range(3):
        NewDat[s+i] = KappaBWriter(NewDat, Kappas, s+i, i, 2)     #Higgs-W-W
        NewDat[s+3+i] = KappaBWriter(NewDat, Kappas, s+3+i, i, 2) #Higgs-Z-Z
        NewDat[s+6+i] = KappaBWriter(NewDat, Kappas, s+6+i, i, 0) #Higgs-gluon-gluon
    #Higgs-higgs-Z couplings from Higgs Hunters guide pg. 360, which has a factor 2 vs the HiggsBounds normalisation, which is then squared
    NewDat[s+12] = "    "+'{:.5f}'.format(cosba/2)+"         3    36    25    23 # higgs-higgs-Z effective coupling, normalised\n"
    NewDat[s+13] = "    "+'{:.5f}'.format(np.sin(ba)/2)+"         3    36    35    23 # higgs-higgs-Z effective coupling, normalised\n"
    
    #Finally, for the fermion couplings (here the coupling values are straightforward)
    s += 19
    
    for i in range(3):
        NewDat[s+i] = KappaFWriter(NewDat, Kappas, s+i, i, 1)     #Higgs-b-b 
        NewDat[s+3+i] = KappaFWriter(NewDat, Kappas, s+3+i, i, 0) #Higgs-top-top
        NewDat[s+6+i] = KappaFWriter(NewDat, Kappas, s+6+i, i, 1) #Higgs-tau-tau
   
    #Clearing the spare lines at the end of the file
    del NewDat[s+9:]
    
    #Writing the new data to the ouput file to use as input for HiggsBounds
    with open(outfile, 'w') as file:
        file.writelines(NewDat)
        
    return

#Iterating over all the points and creating the new files
for i in range(npoints):
    HB_FileWriter(i)
