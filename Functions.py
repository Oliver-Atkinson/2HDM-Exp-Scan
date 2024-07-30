"""Gathering together all the various imports and functions for use in this chain of calculations etc"""

import numpy as np
import shutil
import sys
import random
from scipy.interpolate import interp1d
import matplotlib as mpl
import matplotlib.pyplot as plt
from itertools import cycle
from pathlib import Path

#Doing the straightforward calculations of the parameters for 2HDecay
def HDecayCalc(point):
    typ, tanb, mHp, mH0, mA0, cosba = point[0], point[1], point[2], point[3], point[4], point[5]
    
    b = np.arctan(tanb)
    ba = np.arccos(cosba)
    alpha = b-ba
    
    m12sq = (mH0**2)*np.sin(b)*np.cos(b) #Valid in degenerate + alignment limit, not in all other cases
    
    return typ, tanb, m12sq, alpha, mH0, mA0, mHp


#Coupling calculations for 2HDM of a given type compared to the SM 
def CouplingCalc(typ, tanb, cosba):
    sinba = np.sin(np.arccos(cosba))
    
    #Getting the correct couplings for the type of 2HDM
    if typ == 1:
        #h^0 kappas
        k_h0_u = sinba+cosba/tanb
        k_h0_dl = k_h0_u
        #H^0 kappas
        k_H_u = cosba-sinba/tanb
        k_H_dl = k_H_u
        #A^0 kappas
        k_A_u = 1/tanb
        k_A_dl = -k_A_u
    
    elif typ == 2:
        k_h0_u = sinba+cosba/tanb
        k_h0_dl = sinba-cosba*tanb
        
        k_H_u = cosba-sinba/tanb
        k_H_dl = cosba+sinba*tanb
    
        k_A_u = 1/tanb
        k_A_dl = tanb
            
    #Vector couplings are type independent 
    k_h0_V = sinba
    k_H_V = cosba
    k_A_V = 0
    
    #Hpm couplings not required
    ks_h0 = [k_h0_u, k_h0_dl, k_h0_V] 
    ks_H = [k_H_u, k_H_dl, k_H_V]
    ks_A = [k_A_u, k_A_dl, k_A_V]
    
    Kappas = np.array((ks_h0, ks_H, ks_A))
    
    return Kappas


#Gets the BRs and width for a given particle ID in a 2HDecay outfile
def BRGamGetter(dat, ID):
    #The start line of the relevant block
    start = "DECAY QCD    "+str(ID)
    
    #Scanning over the file to find the start of the block
    n = 0
    for line in dat:
        if start in line:
            #Getting the width from the start line
            stripped = line.strip()
            split = stripped.split()
            Gam = split[3]
            
            BR_start_ID = n + 6 #6 lines between start of block and start of BRs
            
            break 
        #The full file line counter 
        else:
            n += 1
    
    #Now finding all the BR block, which is of variable length, but always ends on a line beginning with "#"
    m = 0
    for line in dat[BR_start_ID:]:
        if line[0] == "#":
            BR_end_ID = BR_start_ID + m
            break
        #BR block line counter
        else:
            m += 1
    
    BRs = dat[BR_start_ID:BR_end_ID]
        
    return BRs, Gam

#Gets the BRs and width for the top quark only (it appears quite differently in the 2HDecay output)
def TopBRGamGetter(dat):
    #The start line of the relevant block
    start = "DECAY         6"
    
    #Scanning over the file to find the start of the block
    n = 0
    for line in dat:
        if start in line:
            #Getting the width from the start line
            stripped = line.strip()
            split = stripped.split()
            Gam = split[2]
            
            BR_start_ID = n + 2 #1 line between start of block and start of BRs
            
            break 
        #The full file line counter 
        else:
            n += 1
    
    #Now finding all the BR block, which is of variable length, but always ends at the end of the file
    m = 0
    for line in dat[BR_start_ID:]:
        if line[5] != " ":
            #BR block line counter
            m += 1
    BR_end_ID = BR_start_ID + m
            
    BRs = dat[BR_start_ID:BR_end_ID]
        
    return BRs, Gam

#Function to write out the bosonic Higgs couplings
def KappaBWriter(dat, sqKappas, line, KapH, KapP):
    newList =list(dat[line])
    #Boson couplings to 5 d.p always have some fixed elements
    newList[4:11] = '{:.5f}'.format(sqKappas[KapH, KapP])
    return ''.join(newList)

#Function to write out the fermionic Higgs couplings
def KappaFWriter(dat, Kappas, line, KapH, KapP):
    newList =list(dat[line])
    #To differentiate between scalar and pseudoscalar couplings
    if KapH == 2:
        newList[29:47] = '{:.16f}'.format(Kappas[KapH, KapP])
    else:
        newList[3:21] = '{:.16f}'.format(Kappas[KapH, KapP])
    return ''.join(newList)

    
#Extracting the useful data from the 2HDecay output file
def HDecay_InfoGetter(npoint):
    #Opening up the 2HDecay output file and extracting the full data
    with open('2HDecay_Outs/In'+str(npoint+1)+"_BR.out", 'r') as file:
        data = file.readlines()
 
    h_BRs, Gamh = BRGamGetter(data, 25)
    H0_BRs, GamH0 = BRGamGetter(data, 35)
    A0_BRs, GamA0 = BRGamGetter(data, 36)
    Hp_BRs, GamHp = BRGamGetter(data, 37)
    t_BRs, Gamt = TopBRGamGetter(data)
    
    BRs = [h_BRs, H0_BRs, A0_BRs, Hp_BRs, t_BRs]
    Gams = [Gamh, GamH0, GamA0, GamHp, Gamt]
        
    return BRs, Gams


#Gets the pertinent information from the output files
def HBExtractor(npoint):
    #Opening up the HiggsBounds output file and extracting the full data
    with open('HB_Outs/Point'+str(npoint)+".slha", 'r') as file:
        data = file.readlines()
    
    #Getting the results
    ratline, chanline = data[-19], data[-17] #Always these line for the ratio and channel details
    ratio = float(ratline.strip().split()[2]) 
    chan = " ".join(chanline.strip().split()[2:-5])
  
    #Including the hard limit from 1303.6065 that HB misses below 43 GeV
    mhp = float(data[32].strip().split()[1])
    if ratio < 1 and mhp <= 43:
        ratio = 10

    return ratio, chan
            
#New data from id cuts
def new_data(dat_list, ids):
    return [dat_list[i] for i in ids]

#Function to take a readfile in the HiggsBounds style and return an interpolation function
def IntMaker(readfile, lskip, mcol, xscol):
    #Opening the given file and reading in the data
    with open(readfile) as f:
        content = f.readlines()
        
    #Extracting the relevant information (mhp and expected xs) 
    mhps, expxss = [], []
    inp = [x.strip() for x in content] 
    for line in inp[lskip:]: #Skipping heading lines
        split = line.split()
        mhps.append(float(split[mcol]))
        expxss.append(float(split[xscol])) 

    #Do interpolation to get an interpolation function - Note that this is in log space, so intfunc returns log(xs)
    intfunc = interp1d(mhps, np.log10(expxss), kind='linear', bounds_error=False, fill_value='extrapolate') 
    
    #mhps_t = np.linspace(10, 5000, 25)
    #xs_t = (intfunc(mhps_t))
    #plt.plot(mhps_t, xs_t)
    #plt.scatter(mhps, np.log10(expxss))
    #plt.savefig("InterObs.pdf")

    return intfunc 

#Function for consistency check between calculated xs and obs xs
def xs_checker(xs_Dat, readfile):
    #Get set number
    nstart, nend = readfile.find("/Set"), readfile.find(".")
    nset = int(readfile[nstart+4:nend])
    outfile_name = "InterpolationData/ScanData/Set"+str(nset)+"_Check.dat"
    
    #Write out headings to file
    outfile = open(outfile_name, 'w')
    outfile.write("mHp tanb xs_mg BR_tb xs*BR_calc xs*BR_inter ratio\n")
    
    #Getting all the info and writing out
    for i in range(20):
        BR_Hpm_tb = -1 #To allow for possible zero branching ratios
        #Extracting BR(H+ -> t b)
        HpmBR_Block = HDecay_InfoGetter(i)[0][3]
        for line in HpmBR_Block:
            if "BR(H+ -> t       bb    )" in line:
                BR_Hpm_tb = float(line.strip().split()[0])
        m, tb, xs_tb, xs_BR = xs_Dat[i]
        prod = xs_tb*BR_Hpm_tb
        rat = prod/xs_BR
        
        #Write out to file mhp, tanb, xs_tb * BR(H+ -> tb), xs_BR
        outfile.write('{:.6g}'.format(m)+" "+'{:.3g}'.format(tb)+" "+'{:.8e}'.format(xs_tb)+" "+'{:.8e}'.format(BR_Hpm_tb)+" "+'{:.8e}'.format(prod)+" "+'{:.8e}'.format(xs_BR)+" "+'{:.8e}'.format(rat)+"\n")
        
    outfile.close()
    
    return

#Get the xs check results
def xscheck_reader(save_name, nfs):
    #The individual files only exist immediately after the scan is done, but this code may be run after this
    rats_HB, rats_Calc = [], []
    
    InterDatDir = "InterpolationData/ScanData"
    scandir = "ScanPointSets"
    
    #Checking if the relevant file exists
    if Path(InterDatDir+"/Set1_Check.dat").is_file():
        savefile = open(save_name,"w") 
        savefile.write("type tanb mHp mH0 mA0 cos(b-a) | ratio_HB Allowed_HB | ratio_calc Allowed_calc | Agree\n")
        
        nBad = 0 #Counters
        nErr = 0
        #For every set, getting the relevant data from the relevant file and then writing it all out to the new file
        for i in range(nfs):
            xsfile = InterDatDir+"/Set"+str(i+1)+"_Check.dat"
            HBfile = scandir+"/Set"+str(i+1)+"_Out.dat"
            xsDat = np.loadtxt(xsfile, skiprows=1, delimiter=' ')
            
            #Opening the files and getting the data
            with open(HBfile) as f:
                content1 = f.readlines()
                
            for j in range(1,21):
                perm_calc, agr = "No", "No"
                nBad +=1
                
                line_p1 = content1[j][:-1] 
                rat_calc = xsDat[j-1][-1]
                
                xsBR_calc = xsDat[j-1][-3]
                rat_obs = float(line_p1.strip().split()[-2])
                xsBR_obs = xsBR_calc / rat_obs
                
                #Checking permittance and agreement
                if rat_calc < 1:
                    perm_calc = "Yes"
                if rat_calc < 0:
                    perm_calc = "Err"
                    rat_calc = -1
                    nErr += 1
                    nBad -= 1
                if line_p1.strip().split()[-1] == perm_calc:
                    agr = "Yes"
                    nBad -= 1
                    
                #print(line_p1.strip().split()[2], '{:.6f}'.format(xsBR_obs))                           
                
                savefile.write(line_p1+"  |  "+'{:.6f}'.format(rat_calc)+" "+perm_calc+"  |  "+agr+"\n")
        savefile.close()
        
        print("There are",nBad,"/",nfs*20,"points with differing results, and",nErr,"for which the BR is not calculated")

    return

#Getting the BRs as a function of mass and writing out to a file 
def BR_DatGetter(readfile_name, npoints, part):
    #Outfile setup
    outdir = "InterpolationData/ScanData/"
    n_end = readfile_name.find(".")
    n = readfile_name[17:n_end]
    outfile_name = outdir+"Set"+str(n)+"_"+part+"BRs.dat" 
    outfile = open(outfile_name, 'w')
    
    #Writing out header to file 
    outfile.write("tanb mHp mH0 mA0 ") 
    
    #Dictionary of numbers for each mass config and the dictionaries of BRs as they appear in the code - a bit nasty but seems to work well
    part_dict = {
        "h0": [0, {"BR(h -> b bb )": 1,"BR(h -> tau+ tau- )": 3,"BR(h -> mu+ mu- )": 4,"BR(h -> s sb )": 2,"BR(h -> c cb )": 0,"BR(h -> g g )": 5,"BR(h -> gam gam )": 6,"BR(h -> Z gam )": 7,"BR(h -> W+ W- )": 9,"BR(h -> Z Z )": 8,"BR(h -> A A )": 10,"BR(h -> Z A )": 11,"BR(h -> H+ H- )": 12, "BR(h -> W- H+ )": 13,"BR(h -> W+ H- )": 14}],
        "H+": [3, {"BR(H+ -> c bb )": 3, "BR(H+ -> tau+ nu_t au )": 8, "BR(H+ -> mu+ nu_m u )": 9, "BR(H+ -> u bb )": 6, "BR(H+ -> u sb )": 7, "BR(H+ -> c db )": 5, "BR(H+ -> c sb )": 4,"BR(H+ -> t bb )": 0,"BR(H+ -> t sb )": 1, "BR(H+ -> t db )": 2, "BR(H+ -> W+ h )": 10, 'BR(H+ -> W+ H )': 11, 'BR(H+ -> W+ A )': 12}], 
        "H0": [1, {"BR(H -> b bb )": 2,"BR(H -> tau+ tau- )": 4,"BR(H -> mu+ mu- )": 5,"BR(H -> s sb )": 3,"BR(H -> c cb )": 1, "BR(H -> t tb )": 0,"BR(H -> g g )": 6,"BR(H -> gam gam )": 7,"BR(H -> Z gam )": 8,"BR(H -> W+ W- )": 10, "BR(H -> Z Z )": 9,"BR(H -> h h )": 11,"BR(H -> A A )": 12,"BR(H -> Z A )": 13, "BR(H -> H+ H- )": 14, "BR(H -> W- H+ )": 15,"BR(H -> W+ H- )": 16}],
        "A0": [2, {"BR(A -> b bb )": 2, "BR(A -> tau+ tau- )": 4, "BR(A -> mu+ mu- )": 5,"BR(A -> s sb )": 3,"BR(A -> c cb )": 1,"BR(A -> t tb )": 0,"BR(A -> g g )": 6, "BR(A -> gam gam )": 7,"BR(A -> Z gam )": 8,"BR(A -> Z h )": 9,"BR(A -> Z H )": 10,"BR(A -> W+ H- )": 11,"BR(A -> W- H+ )": 12}]
    }
    
    
    #Getting the BR dictionary 
    ID_dict = part_dict[part][1]
    
    #Writing out the column headings to file
    for br in range(len(ID_dict)):
        name = list(ID_dict.keys())[list(ID_dict.values()).index(br)]
        outfile.write(name.replace(" ", "")+" ")
    outfile.write("\n")
       
    #For each point, getting the info and writing out to file
    m_Dat = np.loadtxt(readfile_name, skiprows=1, delimiter=' ')
    for i in range(npoints):    
        #Default to zeros for the BRs
        BRs = list(np.zeros(len(ID_dict)))
    
        #Getting all the data
        BR_Blocks = HDecay_InfoGetter(i)[0]
        BR_Block = BR_Blocks[part_dict[part][0]] 
        t, tanb, mHp, mH0, mA0 = m_Dat[i][0], m_Dat[i][1], m_Dat[i][2], m_Dat[i][3], m_Dat[i][4]
    
        #Getting the information from the BR block and inserting it into the right slot in the list
        for line in BR_Block:
            comps = line.strip().split()
            BR_name = " ".join(comps[5:])
            BR_val = float(comps[0])
               
            BR_ID = ID_dict[BR_name]
            BRs[BR_ID] = BR_val
            
        #Writing out to file
        outfile.write('{:.3g}'.format(tanb)+" "+'{:.6g}'.format(mHp)+" "+'{:.6g}'.format(mH0)+" "+'{:.6g}'.format(mA0)+" ")
        for j in range(len(BRs)-1):
            outfile.write('{:.8E}'.format(BRs[j])+" ")
        outfile.write('{:.8E}'.format(BRs[-1])+"\n") #Needed to avoid trailing " " on lines that prevents nice read in later
    
    outfile.close()
    
    return 


#Plots the BRs as a function of mass
def BR_plotter(part, nfs, mconf, t):
    #The individual files only exist immediately after the scan is done, but this code may be run after this    
    InterDatDir = "InterpolationData/ScanData"
    DatsaveDir = "Results/BRs/Data/"
    PlotsaveDir = "Results/BRs/"
    save_name = "T"+t+mconf+"_"+part+"BRs.dat"
    
    #Dictionary of things for each particle; BR names, plot labels, etc.
    part_dict = {"h0": [" BR(h->ccb) BR(h->bbb) BR(h->ssb) BR(h->tau+tau-) BR(h->mu+mu-) BR(h->gg) BR(h->gamgam) BR(h->Zgam) BR(h->ZZ) BR(h->W+W-) BR(h->AA) BR(h->ZA) BR(h->H+H-) BR(h->W-H+) BR(h->W+H-)", r"$\mathcal{B}(h^0 \to xy)$", [r"$c\bar{c}$", r"$b\bar{b}$", r"$s\bar{s}$", r"$\tau^+\tau^-$", r"$\mu^+\mu^-$", r"$gg$", r"$\gamma\gamma$",r"$Z\gamma$", r"$ZZ$", r"$W^+W^-$", r"$A^0A^0$", r"$ZA^0$", r"$H^+H^-$", r"$W^-H^+$", r"$W^+H^-$"]],
                 "H+": [" BR(H+->tbb) BR(H+->tsb) BR(H+->tdb) BR(H+->cbb) BR(H+->csb) BR(H+->cdb) BR(H+->ubb) BR(H+->usb) BR(H+->tau+nu_tau) BR(H+->mu+nu_mu) BR(H+->W+h) BR(H+->W+H) BR(H+->W+A)", r"$\mathcal{B}(H^+ \to xy)$", [r"$t\bar{b}$", r"$t\bar{s}$", r"$t\bar{d}$", r"$c\bar{b}$", r"$c\bar{s}$", r"$c\bar{d}$", r"$u\bar{b}$", r"$u\bar{s}$", r"$\tau^+\nu_\tau$", r"$\mu^+\nu_\mu$", r"$W^+h^0$", r"$W^+H^0$", r"$W^+A^0$"]],
                 "H0": [" BR(H->ttb) BR(H->ccb) BR(H->bbb) BR(H->ssb) BR(H->tau+tau-) BR(H->mu+mu-) BR(H->gg) BR(H->gamgam) BR(H->Zgam) BR(H->ZZ) BR(H->W+W-) BR(H->hh) BR(H->AA) BR(H->ZA) BR(H->H+H-) BR(H->W-H+) BR(H->W+H-)", r"$\mathcal{B}(H^0 \to xy)$", [r"$t\bar{t}$",r"$c\bar{c}$", r"$b\bar{b}$", r"$s\bar{s}$", r"$\tau^+\tau^-$", r"$\mu^+\mu^-$", r"$gg$", r"$\gamma\gamma$",r"$Z\gamma$", r"$ZZ$", r"$W^+W^-$", r"$hh$", r"$A^0A^0$", r"$ZA^0$", r"$H^+H^-$", r"$W^-H^+$", r"$W^+H^-$"]],
                 "A0": [" BR(A->ttb) BR(A->ccb) BR(A->bbb) BR(A->ssb) BR(A->tau+tau-) BR(A->mu+mu-) BR(A->gg) BR(A->gamgam) BR(A->Zgam) BR(A->Zh) BR(A->ZH) BR(A->W+H-) BR(A->W-H+)", r"$\mathcal{B}(A^0 \to xy)$", [r"$t\bar{t}$",r"$c\bar{c}$", r"$b\bar{b}$", r"$s\bar{s}$", r"$\tau^+\tau^-$", r"$\mu^+\mu^-$", r"$gg$", r"$\gamma\gamma$",r"$Z\gamma$", r"$Zh$", r"$ZH^0$", r"$W^+H^-$", r"$W^-H^+$"]]}

    #Dictionary for the configuration setup
    config_dict = {"H+": [1, r"$\log_{10}[m_{H^+}$ / GeV]", r"$m_{H^+} =1$ TeV, $m_{H^0}=m_{A^0}=5$ TeV, $\cos(\beta-\alpha) = 0$"],
                   "H0": [2, r"$\log_{10}[m_{H^0}$ / GeV]", r"$m_{H^0} =1$ TeV, $m_{H^+}=m_{A^0}=5$ TeV, $\cos(\beta-\alpha) = 0$"],
                   "A0": [3, r"$\log_{10}[m_{A^0}$ / GeV]", r"$m_{A^0} =1$ TeV, $m_{H^+}=m_{H^0}=5$ TeV, $\cos(\beta-\alpha) = 0$"],
                   "Degen": [1, r"$\log_{10}[m_{H^+,H^0,A^0}$ / GeV]", r"$m_{H^+}=m_{H^0}=m_{A^0}=1$ TeV, $\cos(\beta-\alpha) = 0$"]} 
    
    #Checking if the relevant file exists
    if Path(InterDatDir+"/Set1_"+part+"BRs.dat").is_file():
        savefile = open(DatsaveDir+save_name,"w") 
        savefile.write("tanb mHp mH0 mA0"+part_dict[part][0]+"\n")
        
        #For every set, getting the relevant data from the relevant file and then writing it all out to the new file
        for i in range(nfs):
            BRfile = InterDatDir+"/Set"+str(i+1)+"_"+part+"BRs.dat"
            with open(BRfile, 'r') as file:
                content = file.readlines()[1:]
            for j in range(20):
                savefile.write(content[j])
        
        savefile.close()
        
    #Reading in from the save file with all the data in and creating a plot
    BR_Dat = np.loadtxt(DatsaveDir+save_name, skiprows=1, delimiter=' ')
    logms = np.log10(BR_Dat[:,config_dict[mconf][0]]) #Keep masses in log space
    logtbs = np.log10(BR_Dat[:,0])
    
    #Plotting
    lines, nl = ["-","--","-.",":"], 0
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['black', '#0F08E0', '#44EBDB', '#FF0001','#F5D561', '#FF0AE0'])#, 'black', 'purple', 'brown','pink', 'orange', 'coral', 'darkgreen','lightblue', 'lime', 'darkgreen', 'turquoise']) 
    plt.figure()
    for i in range(4,len(part_dict[part][2])+4):
        logBRs = np.log10(BR_Dat[:,i])
        logBRs[logBRs == -np.inf] = -500 #Dealing with the log(0) issue
        plt.plot(logtbs, logBRs, lines[int(nl)])
        nl += 1/len(lines)
        if nl >= 4: 
            nl = 0
    #plt.xlabel(config_dict[mconf][1])
    plt.xlabel(r"$\log_{10}[\tan\beta]$")
    plt.ylabel(part_dict[part][1]) 
    plt.xlim([min(logtbs), max(logtbs)])
    plt.ylim([-5, 0.1])
    plt.legend(part_dict[part][2],bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
    plt.yticks(ticks=np.array([-5., -4., -3., -2., -1.,  0.]), labels = [r"$10^{-5}$",r"$10^{-4}$",r"$10^{-3}$",r"$10^{-2}$",r"$10^{-1}$","1"])
    plt.xticks(ticks=np.array([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2, 2.5]))
    plt.title(config_dict[mconf][2])
    plt.savefig(PlotsaveDir+"T"+t+mconf+"_"+part+"BRs.pdf",bbox_inches='tight')
    
    return

#Checks if string is float
def is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

#Function to extract the luminosity (as a float in fb^-1) from a given filename
def LumGetter(filename):
    #pb or fb
    if "pb-1" in filename:
        Lend = filename.find("pb-")
        unit=1e-3
    elif "fb-1" in filename:
        Lend = filename.find("fb-")
        unit=1
    else:
        print("Neither pb nor fb found in",filename)
        return
    
    #The few (~5) characters preceding the unit should include the numerical value of the luminosity
    num = filename[Lend-5:Lend]
    #Iterating over these characters until only a pure number is left - sometimes picks up a "-" as negative, hence abs()
    while is_number(num) == False:
        num = num[1:]
    
    #Scaling for units and making a float
    lum = abs(float(num))*unit

    return lum