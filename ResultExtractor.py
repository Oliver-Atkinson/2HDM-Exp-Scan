"""Takes the HiggsBounds output files and extracts the top-level result - i.e. is the point allowed and how close is it"""

from Functions import *


#Reading in masses and tanbs and cos(b-a)s
readfile_name = sys.argv[1]
npoints = int(sys.argv[2])
EWPT = sys.argv[3]

#Output file
outfile = readfile_name[:-4]+"_Out.dat"

#Reading in BSMPT output if required
if EWPT == "y":
    BSMPT_file = "BSMPT_Output.dat"
    with open(BSMPT_file) as f:
        content = f.readlines()
    BSMPT_out = [x.strip() for x in content] 

#Copying the original output file and adding in the new data
with open(readfile_name, 'r') as file:
    NewDat = file.readlines()[:npoints+1]

        
#Including the final results
NewDat[0] = NewDat[0][:-1]+" | Obsratio Allowed Channel\n" #Added colummns
if EWPT == "y":
    NewDat[0] = NewDat[0][:-1]+" | omega_c/T_c SFOEWPT\n"
    
#Each data point
for i in range(1,npoints+1):
    #Extracting the ratio and checking if the point is allowed
    rat, chan = HBExtractor(i)
    perm = "No"
    if rat < 1:
        perm ="Yes"
    
    #Extracting the ratio of critical temps and checking if a SFOEWPT is allowed (if required)
    if EWPT == "y":
        outdat = BSMPT_out[i].split()
        xi_c = float(outdat[-5])
        SFO = "No"
        if xi_c > 1:
            SFO ="Yes"
    
        NewDat[i] = NewDat[i][:-1]+"    | "+'{:.6g}'.format(rat)+"   "+perm+"   "+chan+"    | "+'{:.6g}'.format(xi_c)+"     "+SFO+"\n"
    
    else:
        NewDat[i] = NewDat[i][:-1]+"    | "+'{:.6g}'.format(rat)+"   "+perm+"   "+chan+"\n"
    
#Writing the final results to the ouput file 
with open(outfile, 'w') as file:
    file.writelines(NewDat)

#print("Completed! See "+outfile+" for top-level results")