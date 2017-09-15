#!/usr/bin/python

from __future__ import print_function
import time
import sys
import os

def div():
    consoleprint("-------------------------")

def consoleprint(*args, **kwargs):
    print(*args, file=sys.stderr, end="\n", **kwargs)

def progress(*args, **kwargs):
    print(*args, file=sys.stderr, end="\r", **kwargs)
    
def quit():
    div()
    exit(1)

def title():
    consoleprint("\n           \033[92m> Multi Conformational Analysis Tool <\033[0m")

def help():
    title()
    consoleprint("           Last updated : 30 May 2017 by Ze Manel\n")
    consoleprint("  \033[93mUSAGE\033[0m:\n  amc.py -i input.pdb [-n] [-bcc] [-v]\n\n")
    consoleprint("  \033[93mPARAMETERS\033[0m:\n  -n   [int] : Save the [int] frames with the lowest potential energies\n  -bcc [int] : Use AM1-BCC charges in topology file, with [int] total charge\n  -v         : Be loud and noisy! Run in verbose mode\n")
    quit()

def confirm(): #Request user confirmation to proceed with the command
    while True:
        confirmation = raw_input(" WARNING : This action will delete all files and folders in %s. Continue ? (y/n) " % (os.getcwd()))
        if confirmation == 'n': 
            div()
            exit(1)
        elif confirmation == 'y': 
            div()
            break
        else:
            div()
            consoleprint("Please enter y or n only.")

def getCurrentDir(): #Get current dir path
    return os.path.split(os.getcwd())[len(os.path.split(os.getcwd())) - 1]

def deleteContentsInCurrentDirExcept(inputPdb): #Delete all contents of current dir except inputPdb
    thisDir = getCurrentDir()
    os.system("rm -r ../%s_safety >/dev/null 2>&1" % (thisDir))
    os.system("mkdir ../%s_safety && mv %s ../%s_safety" % (thisDir, inputPdb, thisDir))
    os.system("rm -rf * >/dev/null 2>&1")
    os.system("mv ../%s_safety/* . && rm -r ../%s_safety/" % (thisDir, thisDir))

def getPdbFromInput(): #Parses input for pdb file, checking if it is valid
    inputPdb = ''
    for argIndex, arg in enumerate(sys.argv):
        if arg == '-h': help()
        elif arg == '-i': 
            inputPdb = sys.argv[argIndex + 1]
            if not inputPdb[-3:] == 'pdb':
                consoleprint("amc.py can only read pdb files. Please check input files.")
                quit()
    if inputPdb == '':
        consoleprint("No pdb file detected for input. Please use '-i' flag to select input file.")
        quit()
    return inputPdb

def setGlobalsFromInput(): # Sets global variables from defaults/input
    global v
    global bcc
    global toSave
    #Set defaults
    v = '>/dev/null 2>&1'
    bcc = ''
    toSave = 10
    #Parse input for parameters
    for argIndex, arg in enumerate(sys.argv):
        if arg == '-v': v = ''
        elif arg == '-bcc': bcc = '-c bcc -nc %d' % (int(sys.argv[argIndex + 1]))
        elif arg == '-n': toSave = int(sys.argv[argIndex + 1])

def generateMol2(inputPdb): #Returns the name of the generated mol2 file (with Antechamber)
    os.system("antechamber -i %s -fi pdb -o %s.mol2 -fo mol2 %s %s" % (inputPdb, inputPdb[:-4], bcc, v))
    os.system("rm -f A*")
    os.system("rm -f sqm.* >/dev/null 2>&1")
    return inputPdb[:-4]+'.mol2'

def generateFrcmod(inputMol2): #Returns the name of the generated frcmod file (with Parmchk)
    os.system("parmchk -i %s -f mol2 -o %s.frcmod" % (inputMol2, inputMol2[:-5]))
    return inputMol2[:-5]+'.frcmod'

def writeSource(inputMol2, inputFrcmod): #Create a default source file for tLeap
    with open('source', 'ab') as filout:
        filout.write("source leaprc.gaff\n")
        filout.write("frcmod = loadamberparams %s\n" % (inputFrcmod))
        filout.write("molecule = loadmol2 %s\n" % (inputMol2))
        filout.write("saveamberparm molecule %s.top %s.crd\n" %(inputMol2[:-5], inputMol2[:-5]))
        filout.write("quit")

def generateTop(inputMol2): #Returns the name of the generated topology file (with tLeap+ACPype)
    os.system("tleap -f %s %s" % ('source', v))
    os.system("acpype -p %s.top -x %s.crd %s" % (inputMol2[:-5], inputMol2[:-5], v))
    os.system("rm -f *.log *.mdp *.crd *.gro *.frcmod *.mol2")
    os.system("rm -f %s.top && rm -f source" % (inputMol2[:-5]))
    outputTop = inputMol2[:-5]+'.top'
    os.system("mv *.top %s" % (outputTop))
    fix(outputTop)
    return outputTop

def fix(top): #Removes numbers from element symbols in topology file (i.e. H11 -> H)
    data = []
    #Parse data
    with open(top, 'r') as filin:
        for line in filin:
            data.append(line)
    #Find [ atoms ] sub-directory
    for n, line in enumerate(data):
        elem = line.split()
        if len(elem) <= 1: continue
        if elem[1] == 'atoms': 
            nstart = n + 2
        elif elem[1] == 'bonds': nend = n -1
    #Remove numbers from element symbols
    for n, line in enumerate(data[nstart:nend]):
        elem = line.split()
        if len(elem[4]) > 1:
            elem[4] = elem[4][:1]
            data[n+nstart] = "%6s %4s %5s %5s %5s %4s %12s %12s %s %s %s\n" % (elem[0], elem[1], elem[2], elem[3], elem[4], elem[5], elem[6], elem[7], elem[8], elem[9], elem[10])
    #Re-write topology file
    with open(top, 'w') as filout:
        for line in data:
            filout.write(line)

#Correct storing of generated files during run in pre-determined directories, according to steps
def storeStep0():
    os.system("mkdir 0_generation && mv *.mdp *.cpt *.edr *.gro *.log *.log *.top *.xtc 0_generation/ >/dev/null 2>&1")

def storeStep1(inputPdb):
    os.system("mkdir 1_heating && mv * 1_heating >/dev/null 2>&1")
    os.system("mv 1_heating/0_generation .")
    os.system("mv 1_heating/*.trr 1_heating/*.tpr 0_generation")
    os.system("mv 1_heating/%s ." % (inputPdb))

def storeStep2(inputPdb):
    os.system("mkdir 2_minimization && mv * 2_minimization >/dev/null 2>&1")
    os.system("mv 2_minimization/0_generation 2_minimization/1_heating 2_minimization/minimizationEnergies.txt 2_minimization/%s ." % (inputPdb))
    
def extractMinimizationEnergy(frame): #Extract final potential energy after minimization from log
    with open(frame+'.log', 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) < 2: continue
            if elem[0] == 'Potential' and elem[1] == 'Energy' and elem[2] == '=': energy = elem[3]
    with open('minimizationEnergies.txt', 'a') as fillout:
        fillout.write("%5s %13s\n" % (frame, energy))

def cleanMdrun(): #Remove unwanted files from the working directory after MDRun
    os.system("rm *.edr *.log *.tpr *.trr")

def sortEnergies(): #Sort extracted minimization energies, returning the first [int] energies
    data = []
    #Parse data
    with open('minimizationEnergies.txt', 'r') as filin:
        for line in filin:
            elem = line.split()
            data.append([elem[0], float(elem[1])])
    #Sort data
    data_sorted = sorted(data, key = lambda x: x[1])
    #Write final txt file
    with open('minimizationEnergies.txt', 'w') as filout:
        filout.write("\n%3s %6s %10s\n" % ('', 'Frame', 'Energy'))
        filout.write(" ____________________\n\n")
        for n, elem in enumerate(data_sorted):
            filout.write("%3d %6s %10.4f\n" % (n, elem[0], elem[1]))
            if n == toSave: filout.write(" ____________________\n\n")
        filout.write(" ____________________\n\n")
    #Return lowest potential energy frames
    lpef = data_sorted[:toSave]
    return lpef

def saveFrames(lpef):
    os.system("mkdir 3_saved")
    for frame in lpef:
        os.system("cp 2_minimization/%s.gro 3_saved" % (frame[0]))
    os.system("mv minimizationEnergies.txt 3_saved")    

if __name__ == '__main__':

    div()
    
    #MDP files location
    heatingMdp = '/home/jpereira/carbons/1_models/0_mdps/hca.mdp'
    minimizationMdp = '/home/jpereira/carbons/1_models/0_mdps/mca.mdp'
    
    #Input parsing and checking
    inputPdb = getPdbFromInput()
    setGlobalsFromInput()
    
    #Pre-cleaning of current directory
    confirm()
    title()
    deleteContentsInCurrentDirExcept(inputPdb)
    
    #Generate input.mol2
    consoleprint("(Step  1/10) Generating mol2 file from input.pdb")
    inputMol2 = generateMol2(inputPdb)

    #Generate input.frcmod
    consoleprint("(Step  2/10) Generating frcmod file from input.mol2")
    inputFrcmod = generateFrcmod(inputMol2)
    
    #Generate source
    consoleprint("(Step  3/10) Writing source code for tleap")
    writeSource(inputMol2, inputFrcmod)
    
    #Generate GMX.top and GMX.crd
    consoleprint("(Step  4/10) Running tleap + acpype for top file generation")
    inputTop = generateTop(inputMol2)

    #Run grompp
    consoleprint("(Step  5/10) Running grompp and generating input.tpr file")
    os.system("gmx grompp -f %s -p %s -c %s -o %s.tpr %s" % (heatingMdp, inputTop, inputPdb, inputPdb[:-4], v))

    
    #Run mdrun
    consoleprint("(Step  6/10) Running mdrun")
    os.system("gmx mdrun -deffnm %s %s" % (inputPdb[:-4], v))
    
    #Generate frames
    consoleprint("(Step  7/10) Generating individual frames in 1_heating/")
    storeStep0()
    os.system("yes 0 2>/dev/null 2>&1 | gmx trjconv -f %s.trr -s %s.tpr -o %s -sep %s" % (inputPdb[:-4], inputPdb[:-4], inputPdb, v))
    storeStep1(inputPdb)
    
    #Minimize each frame
    consoleprint("(Step  8/10) Minimizing frames in 2_minimization/")
    frames = [frame for frame in os.listdir('1_heating')]
    with open('minimizationEnergies.txt', 'w'): pass
    remainingTime = 0.
    totalTimeStart = time.time()
    for frameIndex, frame in enumerate(frames):
        startingTime = time.time()
        progress("             Running minimization on frame %d out of %d (%3.1f%%) - Will end in %2.1fs" % (frameIndex, len(frames)-1, (float(frameIndex)/(len(frames)-1))*100, remainingTime))
        #Run grompp
        os.system("gmx grompp -f %s -p %s -c %s -o %s.tpr %s" % (minimizationMdp, '0_generation/'+inputTop, '1_heating/'+frame, frame[:-4], v))
        os.system("rm mdout.mdp")
        os.system("gmx mdrun -deffnm %s %s" % (frame[:-4], v))
        extractMinimizationEnergy(frame[:-4])
        cleanMdrun()
        remainingTime = (time.time() - startingTime)*(len(frames)-1-frameIndex)
    elapsedTime = time.time() - totalTimeStart
    progress("             Sucessefully minimized %d out of %d frames in %2.1fs                                  " % (len(frames) - 1, len(frames) - 1, elapsedTime))
    storeStep2(inputPdb)
    consoleprint("\n(Step  9/10) Picking conformations with the lowest potential energy")
    lpef = sortEnergies()
    
    #Save frames with the lowest potential energies
    consoleprint("(Step 10/10) Saving final conformations to 3_saved/")
    saveFrames(lpef)
    consoleprint("")
    div()    
        
