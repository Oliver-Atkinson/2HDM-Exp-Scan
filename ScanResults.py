"""Extracts the numerous results from the scan and plots them"""

from Functions import *

npoints = int(sys.argv[1])
nfiles = int(npoints/20)
m_config = sys.argv[3]
typ = sys.argv[2]
scandir = "ScanPointSets"

#Mass configuration
mdict = {
  "H+": 2,
  "H0": 3,
  "A0": 4,
  "Degen": 2} #Aribitary really 

#Bringing all the results together, if not already done so
PlotSaveDir = "Results/"
DatSaveDir = "Results/Data/"
savefile_name = DatSaveDir+"T"+typ+"_"+m_config+"_Points.out"

#xscheck_reader("xsCheck.out", nfiles) #Can remove once xs checks are complete
#for p in ["h0", "H0", "H+", "A0"]: #Can remove this once BRs are plotted
#    BR_plotter(p, nfiles, m_config, typ)  

#The individual files only exist immediately after the scan is done, but this code may be run after this
if Path(scandir+"/Set1_Out.dat").is_file():
    savefile = open(savefile_name,"w") 
    savefile.write("type tanb mHp mH0 mA0 cos(b-a) | Obsratio Allowed Channel\n")
    for i in range(nfiles):
        file = scandir+"/Set"+str(i+1)+"_Out.dat"
        with open(file) as f:
            content = f.readlines()
        for j in range(1,21):
            savefile.write(content[j])
    savefile.close()

#Getting the pertinent data
#Include sorting by mass of particle to make scan output more useful?
tanbs, ms, rats = [], [], []
with open(savefile_name) as f:
    content = f.readlines()[1:]
    for j in range(npoints):
        stripped = content[j].strip()
        split = stripped.split()
        tanb, rat = float(split[1]), float(split[7])
        m = float(split[mdict[m_config]])
        tanbs.append(tanb)
        ms.append(m)
        rats.append(rat)
      
        
#Finding the good and bad datasets, based on if the point is allowed or not
good_ids, bad_ids = [], []
for i in range(len(rats)):
    if rats[i] < 1:
        good_ids.append(i)
    else:
        bad_ids.append(i)
        
tbs_good, ms_good = new_data(tanbs, good_ids), new_data(ms, good_ids)  
tbs_bad, ms_bad = new_data(tanbs, bad_ids), new_data(ms, bad_ids)
#print("Minimum allowed mass:", min(ms_good))

#Plotting
#Label and title dictionaries
capdict = {
  "H+": [r"$m_{H^0}=m_{A^0}=5$ TeV, $\cos(\beta-\alpha) = 0$", r"$\log_{10}[m_{H^+} / GeV]$"],
  "H0": [r"$m_{H^+}=m_{A^0}=5$ TeV, $\cos(\beta-\alpha) = 0$", r"$\log_{10}[m_{H^0} / GeV]$"],
  "A0": [r"$m_{H^+}=m_{H^0}=5$ TeV, $\cos(\beta-\alpha) = 0$", r"$\log_{10}[m_{A^0} / GeV]$"],
  "Degen": [r"$m_{H^+}=m_{H^0}=m_{A^0}$, $\cos(\beta-\alpha) = 0$", r"$\log_{10}[m_{H^+,H^0,A^0} / GeV]$"] }

plt.figure(figsize=(6,5))
#plt.scatter(np.log10(tbs_bad), np.log10(ms_bad),c=[[0.86,0.2,0.12]],s=1)
#plt.scatter(np.log10(tbs_good), np.log10(ms_good),c=[[0,0.35,0.71]],s=1)
plt.scatter(tbs_bad, ms_bad,c=[[0.86,0.2,0.12]],s=3)
plt.scatter(tbs_good, ms_good,c=[[0,0.35,0.71]],s=3)
#plt.xlim([-1.5, 2.5])
#plt.ylim([1, 4])
plt.xlim([1.8, 2.1])
plt.ylim([88, 104])
plt.title(capdict[m_config][0])
#plt.xlabel(r"$\log_{10}[\tan\beta]$")
#plt.ylabel(capdict[m_config][1])
plt.xlabel(r"$\tan\beta$")
plt.ylabel(r"$m_{H^+,H^0,A^0}$ (GeV)")
plt.savefig(PlotSaveDir+"T"+typ+"_"+m_config+".pdf",bbox_inches='tight')
