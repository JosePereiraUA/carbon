#!/usr/bin/python

from __future__ import print_function
import random
import time
import sys
import os

class Atom:
    
    def __init__(self, n = 0., symbol = "NaN", xyz = [0.,0.,0.], res = 0, charge = 0., atomtype = "NaN"):
        self.n = n
        self.symbol = symbol
        self.xyz = xyz
        self.res = res
        self.charge = charge
        self.atomtype = atomtype

def div():
    consoleprint("-------------------------")

def consoleprint(*args, **kwargs):
    print(*args, file=sys.stderr, end="\n", **kwargs)
    
def log(*args, **kwargs):
    logFile = ("../%s_safety/run.log" % (thisDir))
    with open(logFile, 'ab') as f:
        print(*args, file=f, end="\n", **kwargs)
        

def getPosInt(string): #Returns a positive integer
    while True:
        try: out = int(raw_input(string))
        except ValueError:
            consoleprint("Please enter numeric values.")
            continue
        else:
            if out <= 0:
                consoleprint("Please enter numeric values greater than 0.")
                continue
            else: return out

def changeTitle(title): #Change the title of the system to a user defined input
    new_data = []
    with open("topol.top", 'r') as fin:
        data = fin.readlines()
    for line in data:
        elem = line.split()
        if line == 'CRV\n':
            new_data.append(title+'\n')
        else: new_data.append(line)
    with open("topol.top", 'w') as filout:
        filout.writelines(new_data)
        
def checkAdd(current): #Checks if the amount of atoms in the tmp.pdb file increased
    this = 0
    with open("tmp.pdb", 'r') as fin:
        for line in fin:
            elem = line.split()
            if elem[0] == 'ATOM': this += 1
        if this > current: return True
    return False

def getNumberOfAtoms(fileName='tmp.pdb'): #Returns current number of atoms
    this = 0
    with open(fileName, 'r') as fin:
        for line in fin:
            elem = line.split()
            if elem[0] == 'ATOM': this += 1
    return this

def validToContinue(current, fileName, maxAtom):
    nextResidue = getNumberOfAtoms(fileName)
    if current+nextResidue > maxAtom: return False
    else: return True
    
def updateTopol(res, insertionType): #Updates topol to include the itp/rest_itp of res/pharma
    if insertionType == 'res':
        log("Updating topol.top to insert new residue")
        res1, res2, res3, n = str(res), str(res), 'R'+str(res), 1
        lineFinish = ';[ pharmaceuticals ]'
    elif insertionType == 'na':
        log("Updating topol.top to insert new Na+")
        res1 = '%s3_ions' % (pharma_lib)
        res2 = 'na'
        lineFinish = ';[ solvent ]'
    else:
        log("Updating topol.top to insert new", insertionType)
        n = res
        if insertionType == 'cbz':
            res1 = '%s1_cbz' % (pharma_lib)
            res2 = 'cbz'
            res3 = 'CBZ'
        else:
            res1 = '%s2_smx' % (pharma_lib)
            res2 = 'smx-1/smx'
            res3 = 'SMX'
        lineFinish = ';[ ions ]'
    new_data = []
    with open("topol.top", 'r') as fin:
        data = fin.readlines()
    for index, line in enumerate(data):
        if line == lineFinish+'\n':
            if insertionType == 'na': new_data.append('#include "%s/%s.itp"\n\n' % (res1, res2))
            else: new_data.append('#include "%s/%s.itp"\n\n#ifdef REST\n#include "%s/%s_rest.itp"\n#endif\n\n' % (res1, res2, res1, res2))
            new_data.append(line)
        else: new_data.append(line)
    log(" > Inserting itp and rest_itp ...")
    if not insertionType == 'na': new_data.append("%-4s %13d\n" % (res3, n))
    if not insertionType == 'na': log(" > Inserting [ molecules ] reference to", n, "molecules of", res3)
    with open("topol.top", 'w') as filout:
        filout.writelines(new_data)

def updateIonsMdp():
    ionsMdp = mdps+'ions.mdp'
    data = []
    with open(ionsMdp, 'r') as filin:
        for line in filin:
            data.append(line)
    for n, line in enumerate(data):
        elem = line.split()
        if len(elem) <= 0: continue
        if elem[0] == 'include': data[n] = data[n][:-1] + " -I%s\n" % (lib)
    with open('ions.mdp', 'w') as filout:
        for line in data:
            filout.write(line)

def countResidues(carbon):
    residueList = []
    uncountable = ['SOL']
    with open(carbon, 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) <= 1: continue
            if elem[0] == 'ATOM':
                if elem[3][1:] not in residueList and elem[3] not in uncountable:
                    residueList.append(elem[3][1:])
    consoleprint(residueList)
    return len(residueList)
    
def addPharma(pharma, carbonCharge, carbon='no_pharma/6/6.pdb'): #Outputs a .pdb file with a number of pharmaceutical from user input added to the last frame of carbon minimization
    div()
    consoleprint(" > Adding pharmaceuticals", pharma)
    div()
    if not carbon == 'no_pharma/6/6.pdb': os.system("mv %s tmp.pdb" % (carbon))
    else: os.system("cp %s . && mv 6.pdb tmp.pdb" % (carbon))
    log("Accessing %s library ..." % (pharma_lib))
    for ID, number in enumerate(pharma):
        current = getNumberOfAtoms()
        if ID == 0: pharma_type = "1_cbz/cbz"
        else: pharma_type = "2_smx/smx-1/smx"
        #Insert pharmaceutical
        log("Run > gmx insert-molecules -f tmp.pdb -ci %s%s.pdb -o tmp.pdb -try 100 -nmol %d %s" % (pharma_lib, pharma_type, number, v))
        os.system("gmx insert-molecules -f tmp.pdb -ci %s%s.pdb -o tmp.pdb -try 100 -nmol %d %s" % (pharma_lib, pharma_type, number, v))
        if checkAdd(current):
            updateTopol(number, pharma_type[-3:])
            if not ID == 0: carbonCharge += -1
    
    #Solvate
    log("Re-solvating carbon after pharmaceutical insertion ...")
    os.system("gmx solvate -cp tmp.pdb -cs spc216.gro -o 1.pdb -p topol.top %s && rm *# && rm tmp.pdb" % (v))

    #Add ions to keep charge 0
    carbonCharge = 1
    consoleprint("CHARBON CHARGEEEEEEEEEEEEEEEEEEEEE:", carbonCharge)
    updateIonsMdp()
    yesN = countResidues('1.pdb') + 3
    os.system("gmx grompp -f ions.mdp -c 1.pdb -p topol.top -o ions.tpr -maxwarn 1")
    os.system("yes %d | gmx genion -s ions.tpr -o 1.pdb -p topol.top -pname NA -np %d" % (yesN, carbonCharge * -1))
    updateTopol('NaN', 'na')
    os.system("rm -f *# >/dev/null 2>&1")
    os.system("mkdir 1 && mv ions.tpr ions.mdp 1/")
    
def process():
    #2-Res_Minimization | 3-Unres_Minimization | 4-Heating | 5-Pressure | 6-Collect
    for step in range(2, 7):
        #Skip step 5 if in vacuum
        if step == 5 and not sol: continue
        elif step == 6 and not sol:
            log(" >>Skipping step 5 because simulation is being run in vacuum")
            prevFile = str(step-2)
        else:
            log(" >>Starting step", step)
            prevFile = str(step-1)
        if step == 2: inputFile = prevFile+'.pdb'
        else: inputFile = prevFile+'.gro'
        log("Running grompp ...")
        if not sol: vacuum = '_vac'
        else: vacuum = ''
        mdpFile = str(step)+vacuum+'.mdp'
        updateMDP(mdps+mdpFile, mdpFile)
        log("Run > gmx grompp -f %s -p topol.top -c %s -o %s.tpr -maxwarn 1" % (mdpFile, inputFile, step))
        os.system("gmx grompp -f %s -p topol.top -c %s -o %s.tpr -maxwarn 1" % (mdpFile, inputFile, step))
        log("Running mdrun ...")
        rdd = ''
        #if step == 6: rdd = '-rdd 1'
        log("Run > gmx mdrun -deffnm %s -v %s" % (step, rdd))
        os.system("gmx mdrun -deffnm %s -v %s" % (step, rdd))
        log("Writing trajectory (no h20)")
        pbc = '-pbc mol'
        if sol: os.system("yes 1 2>/dev/null 2>&1 | gmx trjconv -f %s.gro -s %s.tpr -o %s_h2o.pdb %s -conect 1" % (step, step, step, pbc))
        log("Writting trajectory (with h2o)")
        os.system("yes 0 2>/dev/null 2>&1 | gmx trjconv -f %s.trr -s %s.tpr -o %s_trj.pdb %s -conect 1" % (step, step, step, pbc))
        os.system("mv mdout.mdp %s_mdout.mdp" % (step))
        log("Storing files in dir", prevFile)
        os.system("mkdir %s >/dev/null 2>&1" % (prevFile))
        os.system("mv %s* %s/ >/dev/null 2>&1 && rm *# >/dev/null 2>&1" % (prevFile, prevFile))
    log("Writting last position pdb file as 6/6.pdb")
    os.system("yes 0 2>/dev/null 2>&1 | gmx trjconv -f 6.gro -s 6.tpr -o 6.pdb %s -conect 1" % (pbc))
    log("Storing files in dir 6")
    os.system("mkdir 6 && mv 6* 6/ >/dev/null 2>&1 && rm *# >/dev/null 2>&1")

def updateMDP(fileName, outputName):
    data = []
    with open(fileName, 'r') as filin:
        for line in filin:
            data.append(line)
    for n, line in enumerate(data):
        elem = line.split()
        if len(elem) > 0 and elem[0] == 'include': 
            data[n] = line[:-1] + lib +'\n'
    with open(outputName, 'w') as filout:
        for line in data:
            filout.write(line)
            

def updatePDB(selectedRes, title, remark): #Updates final pdb to reflect which residue each atom belongs to
    new_data = []
    with open("1.pdb", 'r') as filin:
        for line in filin:
            elem = line.split()
            if elem[0] == 'CRYST1':
                size = line
            if elem[0] == 'ATOM':
                if not elem[3] == 'SOL': resIdentifier = "R"+str(selectedRes[int(elem[4])-1])
                else: resIdentifier = elem[3]
                new_data.append([elem[0], int(elem[1]), elem[2], resIdentifier, int(elem[4]), float(elem[5]), float(elem[6]), float(elem[7]), float(elem[8]), float(elem[9]), elem[2]])
    with open("1.pdb", 'w') as filout:
        filout.write("%-9s %s\n" % ('TITLE', title))
        filout.write("%-9s %s\n" % ('REMARK', remark))
        filout.write("%s" % size)
        filout.write("%-9s %4d\n" % ('MODEL', 1))
        for line in new_data:
            filout.write("%-5s %5d  %-3s %-3s %5d %11.3f %7.3f %7.3f %5.2f %5.2f %11s\n" % (line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10]))
        filout.write("TER\nENDMDL")

def generateCIF(carbon):
    consoleprint("Generating CIF")
    atoms = []
    residues = {}
    #Obtain CRYST1 values and atom N and POS
    with open(carbon, 'r') as filin:
        for line in filin:
            elem = line.split()
            if elem[0] == 'CRYST1':
                size = [float(elem[1]), float(elem[2]), float(elem[3])]
                angle = [float(elem[4]), float(elem[5]), float(elem[6])]
            elif elem[0] == 'ATOM':
                res = int(elem[3][1:])
                if res not in residues: residues[res] = 1
                else: residues[res] += 1
                atom = Atom(residues[res], elem[2], [float(elem[6]), float(elem[7]), float(elem[8])], res)
                atoms.append(atom)
    cellVolume = size[0]*size[1]*size[2]
    #Obtein atom CHARGE and ATOMTYPE
    l = [int(x) for x in os.listdir(lib)]
    atmDict = {}
    for dirIndex in l:
        atmDict[dirIndex] = []
        itp = '%s/%d/%d.itp' % (lib, dirIndex, dirIndex)
        with open(itp, 'r') as filin:
            for _ in range(6): next(filin)
            for line in filin:
                elem = line.split()
                if len(elem) < 1: break
                atmDict[dirIndex].append([int(elem[0]), elem[1], float(elem[6])])
    for atom in atoms:
        for atomIndex in atmDict[atom.res]:
            if atom.n == atomIndex[0]:
                atom.atomtype = atomIndex[1]
                atom.charge = atomIndex[2]
    #Make sure final charge is 0 or save RESIDUAL_CHARGE
    totalCharge = 0
    for atom in atoms:
        totalCharge += atom.charge
    totalCharge = int(totalCharge)
    RESIDUAL_CHARGE = 0
    if not totalCharge == 0:
        RESIDUAL_CHARGE = totalCharge
    #Write carbon.cif            
    with open("carbon.cif", 'w') as filout:
        filout.write("# CIF generated at %s on %s by carbCluster.py (Ze Manel)\n\n" % (time.strftime("%H:%M:%S"), time.strftime("%d/%m/%Y")))
        filout.write("data_00000001\n_symmetry_cell_setting           triclinic\n_symmetry_space_group_name_H-M   'P 1'\n_symmetry_Int_Tables_number      1\n_space_group_name_Hall           'P 1'\nloop_\n_symmetry_equiv_pos_site_id\n_symmetry_equiv_pos_as_xyz\n1\nx,y,z\n")
        filout.write("_cell_length_a %24.3f\n_cell_length_b %24.3f\n_cell_length_c %24.3f\n_cell_angle_alpha %17.0f\n_cell_angle_beta %18.0f\n_cell_angle_gamma %17.0f\n_cell_volume %27.1f\n" % (size[0], size[1], size[2], angle[0], angle[1], angle[2], cellVolume))
        filout.write("loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_charge\n")
        for atom in atoms:
            filout.write("%s %s %8.6f %8.6f %8.6f %+7.4f\n" % (atom.atomtype, atom.symbol, atom.xyz[0]/size[0], atom.xyz[1]/size[1], atom.xyz[2]/size[2], atom.charge))
        filout.write("\n#END")
    
    return int(round(RESIDUAL_CHARGE))

def storeRaspaFiles(step):
    os.system("mkdir %s && mv * %s >/dev/null 2>&1" % (step, step))
    for i in range(1, 7):
        os.system("mv %s/%d . >/dev/null 2>&1" % (step, i))
    if step == 'surfaceArea': os.system("mv %s/carbon.cif ." % (step))
    else: os.system("mv %s/surfaceArea ." % (step))
    os.system("mv %s/topol.top . && mv %s/carbon.dat . && mv %s/input.dat ." % (step, step, step))
        
def extract(string):
    output = []
    with open("Output/System_0/output_carbon_1.1.1_298.000000_0.data", 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) > 0 and elem[0] == string: output.append(line[1:-1])
    return output

def getSurfaceArea(carbon='6/6.pdb'):
    div()
    consoleprint(" > Determining surface area ...")
    div()
    generateCIF(carbon)
    os.system("cp %smaster_SA.input . && mv master_SA.input simulation.input && simulate simulation.input" % (masters))   
    output = extract('Surface')
    storeRaspaFiles('surfaceArea')
    return output

def getVoidFraction(carbon='6/6.pdb'):
    div()
    consoleprint(" > Determining void fraction ...")
    div()
    generateCIF(carbon)
    os.system("cp %smaster_VF.input . && mv master_VF.input simulation.input && simulate simulation.input" % (masters))
    output = extract('[helium]')
    storeRaspaFiles('voidFraction')
    elem = output[0].split()
    output = ["Void Fraction:  %-8.6f +/- %-9.6f [-]" % (float(elem[3]), float(elem[5]))]
    return output

def getAdsorption(CRVtitle, carbon='6/6.pdb'): #Writes 'output.dat' file and extracts adsorption values
    div()
    consoleprint(" > Determining relative adsorption ...")
    div()
    RESIDUAL_CHARGE = generateCIF(carbon)
    Concentrations = []
                           #name  Concentrations
    Concentrations.append(['cbz', ['0', '1', '1', '0', '10', '10', '0', '100', '100', '0', '1000', '1000']])
    Concentrations.append(['smx', ['1', '0', '1', '10', '0', '10', '100', '0', '100', '1000', '0', '1000']])
                      #Title     Concentrations  Press   T    Charge          pH  Cyc  In  Print
    generateMasterAds(CRVtitle, Concentrations, 101235, 300, RESIDUAL_CHARGE, 7, 5000, 50, 50)
    os.system("gcmc.py")

def generateMasterAds(title, concs, pressure, temp, charge, pH, nCycles, inits, printEvery):
    div = '  |  '
    with open('inputAds.dat', 'w') as filout:
        filout.write("\n[ Simulation Title ]\n%s\n\n" % (title))
        filout.write("[ Molecule Definition ]\n%s%s%s\n" % ('Molecule', div, 'Concentration [ppb]'))
        for c in concs:
            filout.write("%-8s%s%s\n" % (c[0], div, ' '.join(c[1])))
        filout.write("\n[ Simulation Parameters ]\n%-18s%s%d\n%-18s%s%d\n%-18s%s%d\n%-18s%s%d\n%-18s%s%d\n%-18s%s%d\n%-18s%s%d" % ('Pressure [Pa]', div, pressure, 'Temperature [K]', div, temp, 'Charge', div, charge, 'pH', div, pH, 'NumberOfCycles', div, nCycles, 'Initialization [%]', div, inits, 'Print Every', div, printEvery))

def addToOutput(body, header): #Adds a body of text to the final 'carbon.dat' file
    log(" > Adding", header, "tab to carbon.dat")
    path = 'carbon.dat'
    with open(path, 'ab') as filout:
        filout.write("%s\n\n" % (header))
        for line in body:
            filout.write("%s\n" % (line))
        filout.write("\n")

def getInfo(carbon='6/6.pdb'): #Get general info regarding the current carbon sample
    body = []
    mass = 0
    weight = {'C': 12.0107, 'O': 15.9994, 'H': 1.00794}
    RESIDUAL_CHARGE = generateCIF(carbon)
    consoleprint(RESIDUAL_CHARGE)
    with open(carbon, 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) > 0 and elem[0] == 'CRYST1':
                size =  [float(elem[1]), float(elem[2]), float(elem[3])]
            if len(elem) > 0 and elem[0] == 'ATOM':
                count = elem[1]
                if len(elem[5]) < 3: res = elem[5]
                else: res = elem[4]
                mass += weight[elem[2]]
    body.append("Sample was created at %s on %s" % (time.strftime("%H:%M:%S"), time.strftime("%d/%m/%Y")))
    body.append("Box volume: %6.2f ang2 (%6.3f x %6.3f x %6.3f)" % (size[0]*size[1]*size[2], size[0], size[1], size[2]))
    body.append("Atom count: %d atoms in %d different residues" % (int(count), int(res)))
    body.append("Mass: %d g mol-1" % (mass))
    body.append("Charge: %d" % (RESIDUAL_CHARGE))
    return body

def elementalAnalysis(carbon='6/6.pdb'):
    body = []
    ea = {'C': 0, 'O': 0, 'H': 0}
    weight = {'C': 12.0107, 'O': 15.9994, 'H': 1.00794}
    with open(carbon, 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) > 0 and elem[0] == 'ATOM':
                ea[elem[2]] += 1
                count = elem[1]
    totalW = 0
    for elem in ea: totalW += ea[elem]*weight[elem]
    body.append("(mol/mol)")
    for elem in ea:
        n = (float(ea[elem])/int(count))*100
        body.append("%s : %5.2f%% (%3d/%3d)" % (elem, n, ea[elem], int(count)))
    body.append("\n(mass/mass)")
    for elem in ea: 
        m = ((ea[elem]*weight[elem])/totalW)*100
        body.append("%s : %5.2f%%" % (elem, m))
    return body

def functionalGroups(recoverDict=False, carbon='6/6.pdb'):
    body, resList = [], []
    fg = {'sp2': 0, 'ch': 0, 'ror': 0, 'co': 0, 'coo': 0, 'cooh': 0}
    with open(carbon, 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) > 0 and elem[0] == 'ATOM':
                res = int(elem[3][1:])
                if res not in resList: resList.append(res)
    for res in resList:
        path = "%s/%d/%d.dat" % (lib, res, res)
        with open(path, 'r') as filin2:
            for line in filin2:
                elem = line.split()
                if len(elem) >= 1 and elem[1] == 'Functional': break
            for line in filin2:
                elem = line.split()
                if len(elem) <= 2: break
                fg[elem[0]] += int(elem[4])
    if recoverDict: return fg
    total = 0
    for group in fg: total += fg[group]
    body.append("(mol/mol)")
    for group in fg:
        if float(fg[group]) == 0: continue
        n = (float(fg[group])/total)*100
        body.append("%s : %5.2f%% (%3d/%3d)" % (group, n, fg[group], total))
    return body

def getParameter(parameter, name):
    if parameter == 'Yes': return True
    elif parameter == 'No': return False
    else: discrepancyError(name)

def error(string):
    consoleprint("\033[91m", string, "\033[0m")
    div()
    exit(1)

def getCarbonCharge(pdbFile):
    fg = functionalGroups(True, pdbFile)
    return fg['coo']*-1

def getCharge(res):
    datFile = "%s/%s/%s.dat" % (lib, res, res)
    with open(datFile, 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) <= 0 or elem[0][0] == ';': continue
            elif elem[0] == 'charge': return int(elem[2])
    error("No charge information found on %s" % (datFile))

def clearTopol():
    dataToRemove = ['NA', 'SOL', 'CBZ', 'SMX']
    data = []
    with open('topol.top', 'r') as filin:
        for n, line in enumerate(filin):
            elem = line.split()
            if len(elem) <= 1:
                data.append(line)
                continue
            if not elem[0] in dataToRemove: data.append(line)
            if elem[1] == 'pharmaceuticals': phaN = n
            elif elem[1] == 'solvent': solN = n
            elif elem[1] == 'ions': ionN = n
    with open('topol.top', 'w') as filout:
        for n, line in enumerate(data):
            if n not in range(phaN+1, ionN) and n not in range(ionN+1, solN):
                filout.write(line)
            elem = line.split()
            if len(elem) > 1 and elem[1] == 'pharmaceuticals': filout.write('\n')
            elif len(elem) > 1 and elem[1] == 'ions': filout.write('\n')          

def processTrueRun(carbon):
    clearTopol()
    carbonCharge = getCarbonCharge(carbon)
    if (pharma[0] + pharma[1]) > 0:
        log("Trying to add", pharma[0], "cbz and", pharma[1], "smx to", carbon)
        addPharma(pharma, carbonCharge, carbon)
        sol = True
        process()
        newDirName = "pharma_cbz_%d_smx_%d" % (pharma[0], pharma[1])
        os.system("mkdir %s" % (newDirName))
        os.system("mv * %s >/dev/null 2>&1" % (newDirName))
        os.system("mv %s/topol.top . && mv %s/input* . && mv %s/carbon.dat . >/dev/null 2>&1" % (newDirName, newDirName, newDirName))
        os.system("mv %s/no_pharma . >/dev/null 2>&1" % (newDirName))

def discrepancyError(parameter):
    consoleprint("\033[91mDiscrepancy detected in %s. Please verify input.dat file\033[0m" % (parameter))
    div()
    exit(1)

def parseMaster(masterCC):
    outList = ['ncrv', 'nres', 'size', 'ncbz', 'nsmx', 'maxAtoms', 'pzc', 'pH']
    out = {}
    path = {}
    lib, preExistingCarbonPDB = '', ''
    preExistingCarbon, sol, surf, void = False, False, False, False
    title, pharma = [], []
    with open(masterCC, 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) <= 0 or elem[0][0] == ';': continue
            if elem[0] == 'STRUCTURE' and len(elem) > 2:
                preExistingCarbon = True
                preExistingCarbonPDB = elem[2]
            if elem[0] in outList:
                if len(elem) <= 1: discrepancyError(elem[0])
                out[elem[0]] = int(elem[1])
            if elem[0] == 'Library':
                if len(elem) <= 1: discrepancyError('library')
                lib = elem[1]
            if elem[0] == 'Solvate': sol = getParameter(elem[1], 'solvate')
            if elem[0] == 'Surface_area': surf = getParameter(elem[1], 'surface area')
            if elem[0] == 'Void_fraction': void = getParameter(elem[1], 'void fraction')
            if elem[0] == 'Adsorption': ads = getParameter(elem[1], 'adsorption')
            if elem[0] == 'Architecture':
                path['masters'] = next(filin).split()[2]
                path['lib'] = next(filin).split()[2]
                path['pharma'] = next(filin).split()[2]
                path['mdps'] = next(filin).split()[2]
        for index in range(out['ncrv']):
            title.append('CRV'+str(index+1))
        pharma = [out['ncbz'], out['nsmx']]
        if out['pH'] > out['pzc']: pzcLib = '-'
        else: pzcLib = '+'
    return [out['ncrv'], out['nres'], out['size'], pharma, title, lib, sol, surf, void, path, preExistingCarbon, preExistingCarbonPDB, out['maxAtoms'], pzcLib, ads]
    
def help():
    consoleprint("\n  > \033[92mFlora\033[0m: Clusters residues from library to create carbon model, then runs minimization. \n    Adds pharmaceuticals and performs molecular mechanics simulation\n  > Last updated: 12 Jun 2017 by Ze Manel\n")
    consoleprint("  \033[93mUSAGE\033[0m:")
    consoleprint("  CAUTION: \033[91mRequires input.dat file.\033[0m\n  flora.py [-\033[92mf\033[0m] [-\033[92ml\033[0m] [-\033[92mo\033[0m] [-\033[92mr\033[0m] [-\033[92ma\033[0m] [-g] [-s] [-p] [-v]\n\n  \033[93mPARAMETERS\033[0m:\n  -h (--help)               : Displays this help screen\n                            : \033[91mUse following options to override values in input.dat\033[0m\n  -\033[92mf\033[0m (--frac)   [int]       : Flip boolean value for void fraction determination in input.dat\n  -\033[92ml\033[0m (--lib)    [str]       : Use [str] library instead\n  -\033[92mo\033[0m (--odd)    [int]       : Try [int] insertions of residues for each carbon sample\n  -\033[92mr\033[0m (--replic) [int]       : Build [int] replicas of the input carbon\n  -\033[92ma\033[0m (--area)               : Flip boolean value for surface area determination in input.dat\n  -g (--gcmc)               : Flip boolean value for adsorption simulation in input.dat\n  -s (--size)   [int]       : Use a cubic box with [int] nm in size for simulation\n  -p (--pharma) [int] [int] : Use [int] cbz and [int] smx molecules on the produced carbon\n  -v (--verb)               : Be loud and noisy! Run script in verbose mode\n")
    div()

if __name__ == '__main__': #-----------------------------------------------------------------------

    div()

    #Global variables
    global RESIDUAL_CHARGE
    RESIDUAL_CHARGE = 0
    
    #Input parsing and verification
    l, surf, void = False, False, False
    v = ">/dev/null 2>&1"
    title = ['NaN']
    for index, arg in enumerate(sys.argv):
        if arg == '-h' or arg == '--help':
            help()
            exit(1)
    masterCC = 'input.dat'
    #Get parameters from input.dat
    if os.path.isfile('input.dat') == False: error("No input.dat file detected in current directory.")
        
    inP = parseMaster(masterCC)
    n, s, p, l = True, True, True, True
    nCrv, nRes, size, pharma, title, path, lib, sol, surf, void, preExistingCarbon, preExistingCarbonPDB, maxAtoms, pzcLib, ads = inP[0], inP[1], inP[2], inP[3], inP[4], inP[9], inP[9]['lib'] + str(inP[5]), inP[6], inP[7], inP[8], inP[10], inP[11], inP[12], inP[13], inP[14]
    masters, mdps, pharma_lib = path['masters'], path['mdps'], path['pharma']
    #Overide parameters from optional flags
    for index, arg in enumerate(sys.argv):
        if arg == '-r' or arg == '--replic': nCrv = int(sys.argv[index + 1])
        if arg == '-o' or arg == '--odd': nRes = [int(sys.argv[index + 1])]
        if arg == '-s' or arg == '--size': size = [int(sys.argv[index + 1])]
        if arg == '-a' or arg == '--area': surf = not surf
        if arg == '-f' or arg == '--frac': void = not void
        if arg == '-g' or arg == '--gibbs': ads = not ads
        if arg == '-p' or arg == '--pharma': pharma = [[int(sys.argv[index + 1])], [int(sys.argv[index + 1])]]
        if arg == '-l' or arg == '--lib': lib = path['lib'] + str(sys.argv[index + 1])
        if arg == '-v' or arg == '--verb': v = ""
    #Check for the presence of topol.top and preExistingCarbonPDB
    if preExistingCarbon:
        if os.path.isfile('topol.top') == False: error("No topol.top file detected. Can't use pre existing carbon.")
        elif os.path.isfile(preExistingCarbonPDB) == False: error("%s not detected. Can't use pre existing carbon." % (preExistingCarbonPDB))
    #Use correct library relative to pH of the simulation
    lib = lib+pzcLib
    
    #Confirmation
    consoleprint("\n\033[92m        > Current input <\033[0m\n")
    if preExistingCarbon: consoleprint(" Using pre existing carbon: %s\n" % (preExistingCarbonPDB))
    consoleprint("%26s : %d\n%26s : %s\n%26s : %s\n%26s : %s\n%26s : %s\n%26s : %s\n%26s : %s\n%26s : %s\n%26s : %s\n%26s : %s\n%26s : %s\n%26s : %s\n" % ('Number of carbon replicas', nCrv, 'Number of res per sample', nRes, 'Cell size per sample', size, 'Number of cbz molecules', pharma[0], 'Number of smx molecules', pharma[1], 'Maximum number of atoms:', maxAtoms, 'Title', ', '.join(str(a) for a in title), 'Library', lib, 'Add water as solvent', sol, 'Compute surface area', surf, 'Compute void fraction', void, 'Simulate adsorption', ads))
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
            
    #Remove existing folders and create safety storage
    this = os.path.split(os.getcwd())
    thisDir = this[len(this)-1]
    os.system("rm -r ../%s_safety >/dev/null 2>&1")
    os.system("mkdir ../%s_safety" % (thisDir))
    if preExistingCarbon: os.system("mv %s ../%s_safety/ && mv topol.top ../%s_safety/" % (preExistingCarbonPDB, thisDir, thisDir))
    os.system("mv input.dat ../%s_safety/ && rm -rf * >/dev/null 2>&1" % (thisDir))
    os.system("mv ../%s_safety/* ." % (thisDir))
    os.system("cp %s input_backup.pdb >/dev/null 2>&1" % (preExistingCarbonPDB))
    
    if not preExistingCarbon:
        #Start log
        log("[ Run log ]\nCarbon production run started at %s on %s by Ze Manel" % (time.strftime("%H:%M:%S"), time.strftime("%d/%m/%Y")))
        log("\nWill attempt to produce", nCrv, "carbon replicas, with the following parameters:")
        log("Number of residues per carbon:", nRes, "\nSize of box per carbon:", size, "\nNumber of pharmaceuticals to add per carbon:", pharma, "\nTitle of each carbon", title, "\n")
        log("--------------")
        
        if sol == True: startingSol = True
        else: startingSol = False
        for carbIndex in range(nCrv):
            carbonCharge = 0
            sol = startingSol
            consoleprint("Starting carbon production: replica %d out of %d" % (carbIndex, nCrv))
            #Select nRes random residues from master pool
            log("Production run number", carbIndex, "\nWorking on carbon", title[carbIndex])
            l = [int(x) for x in os.listdir(lib)]
            l = random.sample(l, len(l))
            log("Acessing residue library at", lib)
            selectedRes = l[0:nRes]
            log("Selected residues:", selectedRes)
            
            #Copy master.top to current directory
            os.system("cp %smaster.top . && mv master.top topol.top" % (masters))
            log("Sucessefully copied topol.top from", masters)
            if not title[carbIndex] == "NaN":
                changeTitle(title[carbIndex])
                log("Changed topol.top title to", title[carbIndex])
            
            #Box first residue
            firstRes = str(selectedRes[0])
            log("Boxing the first residue:", selectedRes[0])
            os.system("gmx editconf -f %s/%s/%s.pdb -o tmp.pdb -c -box %d -bt cubic %s" % (lib, firstRes, firstRes, size, v))
            log("Run > gmx editconf -f %s/%s/%s.pdb -o tmp.pdb -c -box %d -bt cubic %s" % (lib, firstRes, firstRes, size, v))
            current = getNumberOfAtoms()
            log("Carbon has", current, "atoms")
            updateTopol(selectedRes[0], 'res')
            carbonCharge += getCharge(firstRes)

            #Try to insert the following residues
            curRes, success = 0, 0
            for res in selectedRes[1:]:
                log("Boxing residue", res)
                curRes += 1
                nextResidue = '%s/%s/%s.pdb' % (lib, res, res)
                if validToContinue(current, nextResidue, maxAtoms) == False: 
                    log("Adding %s would increase atom count over %d" % (res, maxAtoms))
                    continue
                else: log("Inserting", res)
                os.system("gmx insert-molecules -f tmp.pdb -ci %s -o tmp.pdb -try 100 -nmol 1 %s" % (nextResidue, v))
                log("Run > gmx insert-molecules -f tmp.pdb -ci %s -o tmp.pdb -try 100 -nmol 1 %s" % (nextResidue, v))
                os.system("rm *# %s" % (v))
                if checkAdd(current):
                    log(res, "was sucessefully introduced")
                    carbonCharge += getCharge(res)
                    updateTopol(res, 'res')
                    current = getNumberOfAtoms()
                    log("Carbon has", current, "atoms")
                    success += 1
                else: log(res, "was not introduced in the current carbon")
            log("Molecule insertion finalized. Sucesseful insertions: %d out of %d requested." % (success+1, len(selectedRes)))
            
            #Solvate
            if sol:
                log("Solvating carbon ...")
                log("Run > gmx solvate -cp tmp.pdb -cs spc216.gro -o 1.pdb -p topol.top %s && rm *#" % (v))
                os.system("gmx solvate -cp tmp.pdb -cs spc216.gro -o 1.pdb -p topol.top %s && rm *#" % (v))
            else: os.system("mv tmp.pdb 1.pdb")
            if title[carbIndex] == "NaN": titleP = "Activated Carbon Model"
            else: titleP = title[carbIndex]
            remark = "Model built at %s on %s by Ze Manel" % (time.strftime("%H:%M:%S"), time.strftime("%d/%m/%Y"))
            log("Updating 1.pdb to reflect residue numbers ...")
            updatePDB(selectedRes, titleP, remark)
            os.system("rm -f tmp.pdb")
            
    #Minimization---------------------------------------------------------------------------------
            log("\nInitializing minimization of carbon:\n")
            #2-Res_Minimization | 3-Unres_Minimization | 4-Heating | 5-Pressure | 6-Collect
            process()
            log("Finished minimization steps")       
            t = open('carbon.dat', 'w+')
            t.close()
            addToOutput(getInfo(), '[ Carbon 1 ]')
            addToOutput(elementalAnalysis(), '[ Elemental Analysis ]')
            addToOutput(functionalGroups(), '[ Functional content ]')
            if ads: getAdsorption("CRV%d" % (carbIndex + 1))
            if surf: addToOutput(getSurfaceArea(), '[ Surface Area ]')
            if void: addToOutput(getVoidFraction(), '[ Void Fraction ]')
            if (pharma[0] + pharma[1]) > 0:
                log("Trying to add", pharma[0], "cbz and", pharma[1], "smx to carbon")
                log("Storing existing files in dir 'no_pharma/'")
                os.system("mkdir no_pharma && mv * no_pharma >/dev/null 2>&1")
                os.system("mv no_pharma/topol.top . && mv no_pharma/carbon.dat .")
                addPharma(pharma, carbonCharge)
                sol = True
                process()
                log("Storing existing files in dir 'pharma'")
                os.system("mkdir pharma && mv * pharma >/dev/null 2>&1")
                os.system("mv pharma/topol.top . && mv pharma/no_pharma . && mv pharma/carbon.dat . && mv pharma/input.dat . >/dev/null 2>&1")
            else: log("No pharmaceuticals to add. Continuing process.")
            if nCrv > 1:#Safely store produced carbons
                log("Storing all files in ../%s_safety/CRV%d" % (thisDir, carbIndex+1))
                os.system("mkdir ../%s_safety/CRV%d/ && mv * ../%s_safety/CRV%d/" % (thisDir, carbIndex+1, thisDir, carbIndex+1))
        #Recover stored directories
        log("Recovering files saved in ../%s_safety/\n-------------" % (thisDir))
    
    else:
        processTrueRun(preExistingCarbonPDB)
    
    #Finish
    log("Carbon production run finished at %s on %s" % (time.strftime("%H:%M:%S"), time.strftime("%d/%m/%Y")))
    os.system("mv ../%s_safety/* ." % (thisDir))
    os.system("rm -r ../%s_safety/" % (thisDir))
    os.system("mv CRV1/input.dat . >/dev/null 2>&1")
    os.system("mv no_pharma/input.dat . >/dev/null 2>&1")
    os.system("mv input_backup.pdb input.pdb >/dev/null 2>&1")
    div()
    consoleprint("\n\nCarbon production run sucessefull. %d out of %d requested carbons built." % (nCrv, nCrv))
    div()
            
            
