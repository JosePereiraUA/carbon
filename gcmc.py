#!/usr/bin/python

from __future__ import print_function
from decimal import Decimal as dec
import math
import time
import sys
import os

#Classes
class Atom: 
    def __init__(self, n, ID, symbol, xyz, rigid=True):
        self.n = n
        self.ID = ID
        self.symbol = symbol
        self.xyz = xyz
        self.rigid = rigid

    def addTypeAndCharge(self, atomtype, charge):
        self.atomtype = atomtype
        self.charge = charge
        
    def addLJ(self, V, W):
        self.V = V
        self.W = W

    def addMassAndConnectivity(self, mass, connectivity):
        self.mass = mass
        self.connectivity = connectivity

class Molecule:
    def __init__(self, n, name, quantity, cation, concentrations, fractions):
        self.n = n
        self.name = name
        self.quantity = quantity
        self.cation = cation
        self.concentrations = concentrations
        self.fractions = fractions

    def calculateFraction(self):
        molarMass = {'cbz': 236.3, 'smx': 253.3, 'h2o': 18.02}
        waterDensity = 996.6 #g L-1 at 300K and 1 atm
        fractions = []
        for concentration in self.concentrations:
            fractions.append((((concentration*math.pow(10, -6))/molarMass[self.name])/(waterDensity/molarMass['h2o']))*H2OFRACTION)
        self.fractions = fractions

class Simulation:
    def __init__(self, pressure, temperature, carbonCharge, pH, collectCycles, initCycles, printEvery):
        self.pressure = pressure
        self.temperature = temperature
        self.carbonCharge = carbonCharge
        self.pH = pH
        self.collectCycles = collectCycles
        self.initCycles = initCycles
        self.printEvery = printEvery

class CriticalConstants:
    def __init__(self, temperature, pressure, accentricFactor):
        self.temperature = temperature
        self.pressure = pressure
        self.accentricFactor = accentricFactor

class EnergyValue:
    def __init__(self, name, total, vdw, coulomb):
        self.name = name
        self.total = total
        self.vdw = vdw
        self.coulomb = coulomb

class AdsorptionValue:
    def __init__(self, name, value, meanDeviation):
        self.name = name
        self.value = value
        self.meanDeviation = meanDeviation

class Bond:
    def __init__(self, at1, at2, potential, r, k):
        self.at1 = at1
        self.at2 = at2
        self.potential = potential
        self.r = r
        self.k = k

class Bend:
    def __init__(self, at1, at2, at3, potential, r, k):
        self.at1 = at1
        self.at2 = at2
        self.at3 = at3
        self.potential = potential
        self.r = r
        self.k = k

class ProperDihedral:
    def __init__(self, at1, at2, at3, at4, potential, c1, c2, c3, c4, c5, c6):
        self.at1 = at1
        self.at2 = at2
        self.at3 = at3
        self.at4 = at4
        self.potential = potential
        self.c1 = c1
        self.c2 = c2 
        self.c3 = c3
        self.c4 = c4
        self.c5 = c5
        self.c6 = c6

class ImproperDihedral:
    def __init__(self, at1, at2, at3, at4, potential, phase, kd, pn):
        self.at1 = at1
        self.at2 = at2
        self.at3 = at3
        self.at4 = at4
        self.potential = potential
        self.phase = phase
        self.kd = kd
        self.pn = pn

        
#ToolBox
def consoleprint(*args, **kwargs):
    print(*args, file=sys.stderr, end="\n", **kwargs)

def div():
    consoleprint("-------------------------")

def error(string):
    consoleprint("\033[91m", string, "\033[0m")
    div()
    exit(1)

def checkPath(path):
    if os.path.isfile(path) == False: error("%s file was not found." % (path))

def help():
    consoleprint("\n  > \033[92mGCMC\033[0m: Performs Monte Carlo simulations using RASPA 2.0 software on Grand Canonical Ensemble\n  > Last updated: 14 Jun 2017 by Ze Manel\n")
    consoleprint("  \033[93mUSAGE\033[0m:")
    consoleprint("  Requires carbon.cif and inputAds.dat files in the current working directory\n  gcmc.py [-\033[92mh\033[0m] [-\033[92msmx\033[0m] [-\033[92mcbz\033[0m] [-\033[92md\033[0m] [-\033[92mm\033[0m] [-\033[92mp\033[0m] [-\033[92mc\033[0m]\n\n  \033[93mPARAMETERS\033[0m:\n  -\033[92mh\033[0m         : Displays this help screen\n  -\033[92mcbz\033[0m       : Use cbz to define necessary files\n  -\033[92msmx\033[0m       : Use smx to define necessary files\n  -\033[92md\033[0m         : Write pharma.def file in pharmaceutical library\n  -\033[92mm\033[0m         : Write force_field_mixing_rules.def file in pharmaceutical library\n  -\033[92mp\033[0m         : Write pseudo_atoms.def file in pharmaceutical library\n  -\033[92mpH\033[0m        : Define molecules for this pH\n  -\033[92mt\033[0m         : Run GCMC simulation in TEST mode (Default: 10 steps, print every 1, at 10 Pa)")   
    div()
    exit(1)
    
#Functions
def parseInputFile(inputFile): #Parse inputAds.dat to get simulation settings
    molecules = []
    #Get simulation title
    with open(inputFile, 'r') as filin:
        for line in filin:
            elem = line.split()
            if elem[1:3] == ['Simulation', 'Title']:
                title = next(filin)[:-1]
            if elem[1:3] == ['Molecule', 'Definition']: break
        next(filin)
        #Add pharmaceuticals
        for n, line in enumerate(filin):
            elem = line.split()
            if len(elem) <= 0: break
            concentrations = []
            for c in elem[2:]:
                if c == '|': break
                concentrations.append(int(c))
            molecules.append(Molecule(n, elem[0], 0, 'no', concentrations, []))
        #Verify pharmaceutical concentration count
        verifyConCount(molecules)
        #Calculate fractions based on ppb concentration
        for molecule in molecules:
            molecule.calculateFraction()
        #Get simulation parameters
        for line in filin:
            elem = line.split()
            if elem[1:3] == ['Simulation', 'Parameters']: break
        press = int(next(filin).split()[3])
        temp = int(next(filin).split()[3])
        charge = abs(int(next(filin).split()[2]))
        pH = int(next(filin).split()[2])
        cycles = float(next(filin).split()[2])
        inits = float(next(filin).split()[3])
        printEvery = int(next(filin).split()[3])
        simulation = Simulation(press, temp, charge, pH, (cycles*(1-(inits/100))), (cycles*(inits/100)), printEvery)
        #Add water and cations
        h2oFraction = [H2OFRACTION] * len(molecules[0].concentrations)
        molecules.append(Molecule(2, 'h2o', 0, 'yes', 0, h2oFraction))
        naFraction = getNaFraction(molecules, simulation.pH)
        molecules.append(Molecule(3, 'Na', simulation.carbonCharge, 'yes', 0, naFraction))
        
        return [title, molecules, simulation]
        

def verifyConCount(molecules): #Verify that all mol in inputAds.dat have the same number of essays
    conCount = len(molecules[0].concentrations)
    for molecule in molecules:
        if not len(molecule.concentrations) == conCount: error("Please verify concentration count on %s" % (molecule.name))


def getNaFraction(molecules, pH): #Get Na+ fraction based on smx fractions and simulation pH
    if pH < SMX_PKA: return [0] * len(molecules[0].concentrations)
    else:
        for molecule in molecules:
            if molecule.name == 'smx': return molecule.fractions
    

def parseInput(): #Parse user input to check what files to define
    moleculesToGather = []
    pH = 0
    test, molDef, mix, pseudo = False, False, False, False
    for argIndex, arg in enumerate(sys.argv):
        if arg == '-h': help()
        elif arg == '-smx': moleculesToGather.append('smx')
        elif arg == '-cbz': moleculesToGather.append('cbz')
        elif arg == '-t': test = True
        elif arg == '-d': molDef = True
        elif arg == '-m': mix = True
        elif arg == '-p': pseudo = True
        elif arg == '-pH': pH = int(sys.argv[argIndex + 1])
    if (molDef or pseudo) and pH == 0: error("You must define pH with -pH [int] flag.")
    return [moleculesToGather, test, molDef, mix, pseudo, pH]

            
def writeMolDef(molecules, smxCharge): #Writed mol.def for each molecule requested 
    for molecule in molecules:
        appendix = getAppendix(molecule)
        if molecule == 'smx':
            consoleprint(" > Writting %s.def in %s%s%s/%s%s" % (molecule, pathToLibrary, appendix, molecule, molecule, smxCharge))
            pdbFile = "%s%s%s/%s%s/%s.pdb" % (pathToLibrary, appendix, molecule, molecule, smxCharge, molecule)
            itpFile = "%s%s%s/%s%s/%s.itp" % (pathToLibrary, appendix, molecule, molecule, smxCharge, molecule)
            datFile = "%s%s%s/%s%s/%s.dat" % (pathToLibrary, appendix, molecule, molecule, smxCharge, molecule)
        else:
            consoleprint(" > Writting %s.def in %s%s%s" % (molecule, pathToLibrary, appendix, molecule))
            pdbFile = "%s%s%s/%s.pdb" % (pathToLibrary, appendix, molecule, molecule)
            itpFile = "%s%s%s/%s.itp" % (pathToLibrary, appendix, molecule, molecule)
            datFile = "%s%s%s/%s.dat" % (pathToLibrary, appendix, molecule, molecule)
        #Get AtomList
        molAtoms = getAtoms(pdbFile)
        #Update AtomList with Rigid/Non-Rigid groups
        molAtoms = updateRigidness(molAtoms, datFile)
        groupsOfAtoms = getGroupsOfAtoms(molAtoms)
        #Get critical constants
        critCts = getCriticalConstants(datFile)
        #Get bond list
        bonds = getBondList(itpFile, molAtoms)
        #Get bend list
        bends = getBendList(itpFile, molAtoms)
        #Get propers list
        propers = getPropersList(itpFile, molAtoms)
        #Get impropers list
        impropers = getImpropersList(itpFile, molAtoms)
        #Write .def file
        constructDef(molecule, smxCharge, critCts, molAtoms, groupsOfAtoms, bonds, bends, propers, impropers)
        


def getAppendix(molecule): #welp
    if molecule == 'cbz': return '1_'
    return '2_'


def getAtoms(path): #Return an array of Atoms from a pdb file in the pharmaceutical library
    checkPath(path)
    output = []
    rank = 0
    with open(path, 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) <= 1: continue
            elif elem[0] == 'ATOM':
                output.append(Atom(rank, int(elem[1]), elem[2], [float(elem[5]), float(elem[6]), float(elem[7])]))
                rank += 1
    return output


def updateRigidness(atomArray, datFile): #Updates list of atoms with rigidness based on dat file
    checkPath(datFile)
    flexGroup = []
    with open(datFile, 'r') as filin:
        for line in filin:
            elem = line.split()
            if len(elem) > 0 and elem[0] == "FLEXIBLE_ATOMS": break
        elem = next(filin).split()
        for atomID in elem: flexGroup.append(int(atomID))
    for atom in atomArray:
        if atom.ID in flexGroup: atom.rigid = False
    return atomArray


def getCriticalConstants(datFile): #Parses the pharmaceutical dat file to get it's critical cts
    output = CriticalConstants(0, 0, 0)
    with open(datFile) as filin:
        for line in filin:
            elem = line.split()
            if len(elem) > 0 and elem[0] == "CRITICAL_CONSTANTS": break
        output.temperature = float(next(filin).split()[0])
        output.pressure = float(next(filin).split()[0])
        output.accentricFactor = float(next(filin).split()[0])
    return output


def getGroupsOfAtoms(atomArray): #Checks for the existence of flexible groups in the molecule
    current = atomArray[0].rigid
    output, thisGroup = [], []
    for atom in atomArray:
        if atom.rigid == current: thisGroup.append(atom)
        else:
            current = atom.rigid
            output.append(thisGroup)
            thisGroup = []
            thisGroup.append(atom)
    output.append(thisGroup)
    return output    
    

def getRank(ID, atomArray): #welp 2.0
    for atom in atomArray:
        if atom.ID == ID: return atom.n


def getBondList(itp, atomArray): #Returns a list of bonds with correct parameters from .itp file
    output = []
    with open(itp, 'r') as filin:
        #Find [ bonds ] header
        for line in filin:
            elem = line.split()
            if len(elem) > 1 and elem[1] == 'bonds': break
        next(filin)
        #Parse list of bonds and parameters
        for line in filin:
            elem = line.split()
            if len(elem) == 0: break
            at1, at2, rk = int(elem[0]), int(elem[1]), ((float(elem[3])*10), ((float(elem[4])*10)/(6.022*1.38)))
            flexible = checkIfRigid([at1, at2], atomArray)
            if not flexible: bondType = 'RIGID_BOND'
            else: bondType = 'HARMONIC_BOND'
            at1 = getRank(at1, atomArray)
            at2 = getRank(at2, atomArray)
            output.append(Bond(at1, at2, bondType, rk[0], rk[1]))
    return output


def getBendList(itp, atomArray): #Returns a list of bends with correct parameters from .itp file
    output = []
    with open(itp, 'r') as filin:
        #Find [ angles ] header
        for line in filin:
            elem = line.split()
            if len(elem) > 1 and elem[1] == 'angles': break
        next(filin)
        #Parse list of angles and parameters
        for line in filin:
            elem = line.split()
            if len(elem) == 0: break
            at1, at2, at3, rk = int(elem[0]), int(elem[1]), int(elem[2]), (float(elem[4]), float(elem[5])*1000/(6.022*1.38))
            #Only use bends who have at least one flexible atom
            if checkIfRigid([at1, at2, at3], atomArray):
                bendType = 'HARMONIC_BEND'
                at1, at2, at3 = getRank(at1, atomArray), getRank(at2, atomArray), getRank(at3, atomArray)
                thisBend = Bend(at1, at2, at3, bendType, rk[0], rk[1])
                output.append(thisBend)
    return output


def getPropersList(itp, atomArray): #Returns a list of propers w/ correct parameters from itp file
    output = []
    with open(itp, 'r') as filin:
        #Find [ dihedrals ] header (propers)
        for line in filin:
            elem = line.split()
            if len(elem) >= 5 and elem[4] == 'propers': break
        for _ in range(2): next(filin)
        #Parse list of propers and parameters
        for line in filin:
            elem = line.split()
            if len(elem) == 0: break
            at1, at2, at3, at4 = int(elem[0])-1, int(elem[1])-1, int(elem[2])-1, int(elem[3])-1
            if checkIfRigid([at1, at2, at3, at4], atomArray):
                torsionType = 'SIX_COSINE_DIHEDRAL'
                params = [((float(x)*1000)/(1.38*6.022)) for x in elem[5:11]]
                at1, at2, at3, at4 = getRank(at1, atomArray), getRank(at2, atomArray), getRank(at3, atomArray), getRank(at4, atomArray)
                output.append(ProperDihedral(at1, at2, at3, at4, torsionType, params[0], params[1], params[2], params[3], params[4], params[5]))
    return output


def getImpropersList(itp, atomArray): #Returns a list of propers w/ correct parameters from itp file
    output = []
    with open(itp, 'r') as filin:
        #Find [ dihedrals ] header (impropers)
        for line in filin:
            elem = line.split()
            if len(elem) >= 5 and elem[4] == 'impropers': break
        for _ in range(2): next(filin)
        #Parse list of impropers and parameters
        for line in filin:
            elem = line.split()
            if len(elem) == 0: break
            at1, at2, at3, at4 = int(elem[0])-1, int(elem[1])-1, int(elem[2])-1, int(elem[3])-1
            if checkIfRigid([at1, at2, at3, at4], atomArray, 2):
                torsionType = 'CVFF_IMPROPER_DIHEDRAL'
                phase, kd, pn = float(elem[5]), float(elem[6])*1000/(1.38*6.022), int(elem[7])
                at1, at2, at3, at4 = getRank(at1, atomArray), getRank(at2, atomArray), getRank(at3, atomArray), getRank(at4, atomArray)
                output.append(ImproperDihedral(at1, at2, at3, at4, torsionType, phase, kd, pn))
    return output


def checkIfRigid(atomsToSearch, atomArray, x=1): #Checks if at least X atoms in the array are flexible
    flex = 0
    for atom in atomArray:
        for at in atomsToSearch:
            if atom.ID == at and not atom.rigid:
                flex += 1
                if flex == x: return True
    return False
        

def getTypeOfAtom(atomArray, groupDef): #Return an array of atoms of the same type (rigid/flexible)
    if groupDef == 'Rigid': requested = True
    else: requested = False
    return [x for x in atomArray if x.rigid == requested]


def constructDef(molecule, smxCharge, critCts, atomsArray, groupsOfAtoms, bonds, bends, propers, impropers):
    appendix = getAppendix(molecule)
    numberOfAtoms = len(atomsArray)
    if molecule == 'smx': path = "%s%s%s/%s%s/%s.def" % (pathToLibrary, appendix, molecule, molecule, smxCharge, molecule)
    else: path = "%s%s%s/%s.def" % (pathToLibrary, appendix, molecule, molecule)
    with open(path, 'w') as filout:
        #Critical constants
        filout.write("# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [-]\n")
        filout.write("%.2f\n%.2f\n%.4f\n" % (critCts.temperature, critCts.pressure, critCts.accentricFactor))
        #Number of atoms
        filout.write("# Number Of Atoms\n%d\n" % (numberOfAtoms))
        #Number of groups
        filout.write("# Number Of Groups\n%d\n" % (len(groupsOfAtoms)))
        #Groups
        for groupID, group in enumerate(groupsOfAtoms):
            #Group definition
            if group[0].rigid: groupDef = 'Rigid'
            else: groupDef = 'Flexible'
            filout.write("# Alkane-group\n%s\n" % (groupDef))
            #Number of atoms in this group
            filout.write("# number of atoms\n%d\n" % (len(group)))
            #Atomic positions of atoms in this group
            filout.write("# atomic positions\n")
            for atom in group:
                if molecule == 'smx': atomName = "%s%s%s%s" % (molecule.upper(), smxCharge, atom.symbol, atom.ID)
                else: atomName = "%s%s%s" % (molecule.upper(), atom.symbol, atom.ID)
                if atom.rigid: filout.write("%-2d %-8s %6.3f %6.3f %6.3f\n" % (atom.n, atomName, atom.xyz[0], atom.xyz[1], atom.xyz[2]))
                else: filout.write("%-2d %-8s\n" % (atom.n, atomName))
        #Number of parameters
        filout.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % ('# Chiral centers', 'Bond', 'BondDipoles', 'Bend', 'UrayBradley', 'InvBend', 'Torsion', 'Imp. Torsion', 'Bond/Bond', 'Stretch/Bend', 'Bend/Bend', 'Stretch/Torsion', 'Bend/Torsion', 'IntraVDW', 'IntraCoulomb'))
        filout.write("%16d %4d %11d %4d %11d %7d %7d %12d %9d %12d %9d %15d %12d %8d %12d\n" % (0, len(bonds), 0, len(bends), 0, 0, len(propers), len(impropers), 0, 0, 0, 0, 0, 0, 0))
        #Bonds
        filout.write("# Bond stretch: atom n1-n2, type, parameters\n")
        for bond in bonds:
            if bond.potential == 'RIGID_BOND': filout.write("%-2d %-2d %s\n" % (bond.at1, bond.at2, bond.potential))
            else: filout.write("%-2d %-2d %s %10.3f %6.4f\n" % (bond.at1, bond.at2, bond.potential, bond.k, bond.r))
        #Bends
        filout.write("# Bond bending: atom n1-n2-n3, type, parameters\n")
        for bend in bends:
            filout.write("%-2d %-2d %-2d %s %10.3f %6.4f\n" % (bend.at1, bend.at2, bend.at3, bend.potential, bend.k, bend.r))
        #Propers
        filout.write("# Torsion n1-n2-n3-n4 type\n")
        for proper in propers:
            filout.write("%-2d %-2d %-2d %-2d %s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" % (proper.at1, proper.at2, proper.at3, proper.at4, proper.potential, proper.c1, proper.c2, proper.c3, proper.c4, proper.c5, proper.c6))
        #Impropers
        filout.write("# Imp. Torsion n1-n2-n3-n4 type\n")
        for improper in impropers:
            filout.write("%-2d %-2d %-2d %-2d %s %10.5f %d %5.1f\n" % (improper.at1, improper.at2, improper.at3, improper.at4, improper.potential, improper.kd, improper.pn, improper.phase))
        #Number of config moves
        rigidGroups = getRigidGroups(groupsOfAtoms)
        numberOfRigidGroups = len(rigidGroups)
        filout.write("# Number of config moves\n%d\n# nr_fixes followed by a list\n" % (numberOfRigidGroups))
        for group in rigidGroups:
            filout.write("%d %s\n" % (len(group), ' '.join(str(x.n) for x in group)))
        

def getRigidGroups(groups):
    output = []
    for group in groups:
        if group[0].rigid: output.append(group)
    return output


def writeMixRules(molecules, smxCharge): #Writes force_field_mixing_rules.def (parameters from ffnonbonded)
    for molecule in molecules:
        appendix = getAppendix(molecule)
        if molecule == 'smx':
            consoleprint(" > Writting %s.def for %s in %s%s%s/%s%s" % ('force_field_mixing_rules', molecule.upper(), pathToLibrary, appendix, molecule, molecule, smxCharge))
            pdbFile = "%s%s%s/%s%s/%s.pdb" % (pathToLibrary, appendix, molecule, molecule, smxCharge, molecule)
            itpFile = "%s%s%s/%s%s/%s.itp" % (pathToLibrary, appendix, molecule, molecule, smxCharge, molecule)
            path = "%s%s%s/%s%s/force_field_mixing_rules.def" % (pathToLibrary, appendix, molecule, molecule, smxCharge)
        else:
            consoleprint(" > Writting %s.def for %s in %s%s%s" % ('force_field_mixing_rules', molecule.upper(), pathToLibrary, appendix, molecule))
            pdbFile = "%s%s%s/%s.pdb" % (pathToLibrary, appendix, molecule, molecule)
            itpFile = "%s%s%s/%s.itp" % (pathToLibrary, appendix, molecule, molecule)
            path = "%s%s%s/force_field_mixing_rules.def" % (pathToLibrary, appendix, molecule)
        #Get Atom List
        molAtoms = getAtoms(pdbFile)
        #Add atomtypes to Atom List
        molAtoms = addTypesAndCharges(itpFile, molAtoms)
        #Add Lennard-Jones parameters to Atom List
        molAtoms = addLJ(molAtoms)
        #Write output file
        with open(path, 'w') as filout:
            filout.write("# general rule for shifted vs truncated\n%s\n" % ('truncated'))
            filout.write("# general rule tailcorrection\n%s\n" % ('no'))
            filout.write("# number of defined interactions\n%d\n" % (len(molAtoms)))
            filout.write("# type interaction\n")
            for atom in molAtoms:
                if molecule == 'smx': atomName = "%s%s%s%s" % (molecule.upper(), smxCharge, atom.symbol, atom.ID)
                else: atomName = "%s%s%s" % (molecule.upper(), atom.symbol, atom.ID)
                filout.write("%8s %20s %9.3f %9.3f\n" % (atomName, 'LENNARD_JONES', atom.W, atom.V))
            filout.write("# general mixing rule for Lennard-Jones\n%s" % ('Lorentz-Berthelot'))
        
               
def addTypesAndCharges(itp, atomArray): #Adds atomtypes and charges to atom's list based on itp
    with open(itp, 'r') as filin:
        #Find [ atoms ] header
        for i, line in enumerate(filin):
            elem = line.split()
            if len(elem) >= 2 and elem[1] == 'atoms': break
        next(filin)
        #Parse list of atoms and add atomtypes to the list
        for line in filin:
            elem = line.split()
            if len(elem) == 0: break
            for atom in atomArray:
                if atom.ID == int(elem[0]):
                    atom.addTypeAndCharge(elem[1], float(elem[6]))
    return atomArray
    

def addLJ(atomArray): #Searched ffnonbonded for Lennard-Jones parameters for each atom
    with open(pathToForcefield, 'r') as filin:
        #Find 'gaff parameters' header
        for line in filin:
            elem = line.split()
            if len(elem) >= 1 and elem[1] == 'gaff': break
        #Add parameters to each atom if they match the atomtype being read
        for line in filin:
            elem = line.split()
            for atom in atomArray:
                if elem[0] == atom.atomtype: 
                    atom.addLJ(float(elem[5])*10, float(elem[6])*1000/(1.38*6.022))
    return atomArray


def writePseudoAtoms(molecules, smxCharge): #Writes pseudo_atoms.def
    for molecule in molecules:
        appendix = getAppendix(molecule)
        if molecule == 'smx': 
            pdbFile = "%s%s%s/%s%s/%s.pdb" % (pathToLibrary, appendix, molecule, molecule, smxCharge, molecule)
            itpFile = "%s%s%s/%s%s/%s.itp" % (pathToLibrary, appendix, molecule, molecule, smxCharge, molecule)
            consoleprint(" > Writting %s.def for %s in %s%s%s/%s%s" % ('pseudo_atoms', molecule.upper(), pathToLibrary, appendix, molecule, molecule, smxCharge))
            path = "%s%s%s/%s%s/%s.def" % (pathToLibrary, appendix, molecule, molecule, smxCharge, 'pseudo_atoms')
        else:
            pdbFile = "%s%s%s/%s.pdb" % (pathToLibrary, appendix, molecule, molecule)
            itpFile = "%s%s%s/%s.itp" % (pathToLibrary, appendix, molecule, molecule)
            consoleprint(" > Writting %s.def for %s in %s%s%s" % ('pseudo_atoms', molecule.upper(), pathToLibrary, appendix, molecule))
            path = "%s%s%s/%s.def" % (pathToLibrary, appendix, molecule, 'pseudo_atoms')
        #Get Atom List
        molAtoms = getAtoms(pdbFile)
        #Add charges
        molAtoms = addTypesAndCharges(itpFile, molAtoms)
        #Add masses
        molAtoms = addMassesAndConnectivities(molAtoms)
        #Write output file
        with open(path, 'w') as filout:
            #Number of pseudo atoms
            filout.write("#number of pseudo atoms\n%d\n" % (len(molAtoms)))
            #Parameters
            filout.write("%-10s %10s %4s %7s %10s %10s %12s %14s %8s %5s %13s %11s %16s %13s\n" % ('#type', 'print', 'as', 'chem', 'oxidation', 'mass', 'charge', 'polarization', 'B-factor', 'radii', 'connectivity', 'anisotropic', 'anisotropic-type', 'tinker-type'))
            for atom in molAtoms:
                if molecule == 'smx': atomName = "%s%s%s%s" % (molecule.upper(), smxCharge, atom.symbol, atom.ID)
                else: atomName = "%s%s%s" % (molecule.upper(), atom.symbol, atom.ID)
                filout.write("%-10s %10s %4s %7s %10d %10.4f %12f %14.1f %8.1f %5.1f %13d %11d %16s %13d\n" % (atomName, 'yes', atom.symbol, atom.symbol, 0, atom.mass, atom.charge, 0.0, 1.0, 1.0, atom.connectivity, 0, 'absolute', 0))
        

def addMassesAndConnectivities(atomArray): #Adds masses and connectivities to atom list
    masses = {'C': 12.0107, 'O': 15.9994, 'H': 1.00794, 'N': 14.0067, 'S': 32.065}
    connectivities = {'C': 4, 'O': 2, 'H': 1, 'N': 3, 'S': 6}
    for atom in atomArray:
        atom.addMassAndConnectivity(masses[atom.symbol], connectivities[atom.symbol])
    return atomArray


def getSwapableMolecules(moleculeList, essay):
    swapableMolecules = []
    for molecule in moleculeList:
        if not molecule.fractions[essay] == 0 and molecule.cation == 'no':
            swapableMolecules.append(molecule.n)
    return swapableMolecules


def fixMoleculeN(moleculeList, essay):
    nmole = -1
    for molecule in moleculeList:
        if not molecule.fractions[essay] == 0:
            nmole += 1
            molecule.n = nmole
    return moleculeList


def writeSimulationInput(moleculeList, simulation, essay): #Writes simulation.input for RASPA
    moleculeList = fixMoleculeN(moleculeList, essay)
    swapableMolecules = getSwapableMolecules(moleculeList, essay)
    with open("simulation.input", 'w') as filout:
        #Add simulation parameters
        filout.write("%-37s %s\n" % ('SimulationType', 'MonteCarlo'))
        filout.write("%-37s %d\n" % ('NumberOfCycles', simulation.collectCycles))
        filout.write("%-37s %d\n" % ('NumberOfInitializationCycles', simulation.initCycles))
        filout.write("%-37s %d\n" % ('PrintEvery', simulation.printEvery))
        filout.write("%-37s %d\n" % ('WriteBinaryRestartFileEvery', simulation.printEvery))
        filout.write("%-37s %s\n" % ('PrintForcefieldToOutput', 'no'))
        filout.write("%-37s %s\n" % ('PrintPseudoAtomsToOutput', 'no'))
        filout.write("%-37s %s\n\n" % ('PrintMoleculeDefinitionToOutput', 'no'))
        filout.write("%-37s %s\n\n" % ('Forcefield', 'MC_Forcefield'))
        #Add carbon box
        filout.write("%-37s %d\n" % ('Framework', 0))
        filout.write("%-37s %s\n" % ('FrameworkName', 'carbon'))
        filout.write("%-37s %d %d %d\n" % ('UnitCells', 1, 1, 1))
        filout.write("%-37s %s\n" % ('UseChargesFromCIFFile', 'yes'))
        filout.write("%-37s %.1f\n" % ('ExternalTemperature', simulation.temperature))
        filout.write("%-37s %s\n" % ('ExternalPressure', simulation.pressure))
        filout.write("%-37s %s\n" % ('Movies', 'yes'))
        filout.write("%-37s %d\n\n" % ('WriteMoviesEvery', simulation.printEvery))
        #Add pharmaceuticals
        nmole = -1
        for molecule in moleculeList:
            if not molecule.fractions[essay] == 0:
                nmole += 1
                filout.write("%-9s %d %-25s %s\n" % ('Component', nmole, 'MoleculeName', molecule.name))
                filout.write("%-11s %-25s %s\n" % ('', 'ExtraFrameworkMolecule', molecule.cation))
                if molecule.fractions[essay] > 1000 or molecule.fractions[essay] == 0:
                    filout.write("%-11s %-25s %.0f\n" % ('', 'MolFraction', molecule.fractions[essay]))
                else:
                    filout.write("%-11s %-25s %.4f\n" % ('', 'MolFraction', molecule.fractions[essay]))
                filout.write("%-11s %-25s %.0f\n" % ('', 'TranslationProbability', 1))
                filout.write("%-11s %-25s %.0f\n" % ('', 'SwapProbability', 1))
                if molecule.cation == 'no':
                    filout.write("%-11s %-25s %.0f\n" % ('', 'RotationProbability', 1))
                    filout.write("%-11s %-25s %.0f\n" % ('', 'ReinsertionProbability', 1))
                    filout.write("%-11s %-25s %.0f\n" % ('', 'IdentityChangeProbability', 1))
                    filout.write("%-13s %-23s %.0f\n" % ('', 'NumberOfIdentityChanges', len(swapableMolecules)))
                    filout.write("%-13s %-23s %s\n" % ('', 'IdentityChangesList', ' '.join(str(molecule) for molecule in swapableMolecules)))
                filout.write("%-11s %-25s %d\n\n" % ('', 'CreateNumberOfMolecules', molecule.quantity))


if __name__ == '__main__':
    
    div()
    
    #Global variable definitions
      #Paths
    global gcmcInput
    gcmcInput = 'inputAds.dat'
    global pathToLibrary
    pathToLibrary = '/home/jpereira/carbons/1_models/1_pharma/'
    global pathToForcefield
    pathToForcefield = '/home/jpereira/carbons/1_models/sheet_ga.ff/ffnonbonded.itp'
    global gcmcOutput
    gcmcOutput = 'output.dat'
      #Raspa input variables
    global H2OFRACT
    H2OFRACTION = 1000000000
    global SMX_PKA
    SMX_PKA = 6.2
    

    #Get input from user
    userInput = parseInput()
    moleculesToGather, test, molDef, mix, pseudo, pH = userInput[0], userInput[1], userInput[2], userInput[3], userInput[4], userInput[5]
    #Set smx charge
    if pH > SMX_PKA: smxCharge = '-1'
    else: smxCharge = '0'
    #Write option files for pharmaceutical definition
    if molDef: writeMolDef(moleculesToGather, smxCharge)   
    if mix: writeMixRules(moleculesToGather, smxCharge)
    if pseudo: writePseudoAtoms(moleculesToGather, smxCharge)
    # Don't run raspa if producing definition files
    if molDef or mix or pseudo: exit(1)
    #Parse input.dat file for input
    fileInput = parseInputFile(gcmcInput)
    title, molecules, simulation = fileInput[0], fileInput[1], fileInput[2]
    consoleprint(" > GCMC simulation title: %s" % (title))
    #Refresh smx charge
    if simulation.pH > SMX_PKA: smxCharge = '-1'
    else: smxCharge = '0'
    #Override raspa input parameters for TEST mode
    if test:
        consoleprint(" > NOTE: Running GCMC simulation in TEST mode")
        simulation = Simulation(simulation.pressure, simulation.temperature, simulation.carbonCharge, simulation.pH, 10, 0, 1)
    #Write simulation.input
    for essay in range(len(molecules[0].concentrations)):
        writeSimulationInput(molecules, simulation, essay)
        #Store simulation.input
        dirName = '%s_%s_%s' % (title, molecules[0].concentrations[essay], molecules[1].concentrations[essay])
        if os.path.isdir(dirName): os.system("rm -rf %s" % (dirName))
        os.system("mkdir %s" % (dirName))
        os.system("mv simulation.input %s" % (dirName))
        #Copy necessary files to this dir
        os.system("cp %s%s/%s.def %s" % (pathToLibrary, '1_cbz', 'cbz', dirName))
        os.system("cp %s%s/%s/%s.def %s" % (pathToLibrary, '2_smx', 'smx' + smxCharge, 'smx', dirName))
        os.system("cp %s%s/%s.def %s" % (pathToLibrary, '3_ions', 'h2o', dirName))
        os.system("cp %s%s/%s.def %s" % (pathToLibrary, '3_ions', 'Na', dirName))
        os.system("cp carbon.cif %s" % (dirName))
    div()
    
