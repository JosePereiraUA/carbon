#!/usr/bin/python

from __future__ import print_function
import sys
import os

def consoleprint(*args, **kwargs):
    print(*args, file=sys.stderr, end="\n", **kwargs)

def getInitialConcentrations(folderName):
    charList = list(folderName)
    i = 0
    group = [[]]
    for char in charList:
        if char == '_':
            group.append([])
            i += 1
        else: group[i].append(char)
    return (int(''.join(group[1])), int(''.join(group[2])))      
    

def orderFolderList(folderList):
    folderList = sorted(folderList)
    newFolderList = folderList[:3]
    group2 = list(reversed([x for x in folderList[3:] if getInitialConcentrations(x)[1] == 0]))
    newFolderList = newFolderList + group2
    group3 = list(reversed([x for x in folderList if x not in newFolderList]))
    newFolderList = newFolderList + group3
    return newFolderList


def parseInput():
    o = False
    for arg in sys.argv:
        if arg == '-o': o = True
    return o


class Pharmaceutical:
    def __init__(self, mol=0, mol_std=0, mgg=0, mgg_std=0, done=False):
        self.mgg = mgg
        self.mgg_std = mgg_std
        self.mol = mol
        self.mol_std = mol_std
        self.done = done

class Simulation:
    def __init__(self, cbzI, smxI, components={'cbz': Pharmaceutical(), 'smx': Pharmaceutical(), 'h2o': Pharmaceutical(), 'Na': Pharmaceutical()}):
        self.cbzI = cbzI
        self.smxI = smxI
        self.components = components
    
if __name__ == '__main__':

    o = parseInput()
    
    outputFile = 'output.log'
    os.system("rm %s" % (outputFile))
    folderList = [x for x in os.listdir(os.getcwd()) if os.path.isdir(x)]
    if o: folderList = orderFolderList(folderList) #Order folder list
    simulations = []
    #Parse input in all folders on current directory
    for folder in folderList:
        if o: (cbzI, smxI) = getInitialConcentrations(folder)
        else: (cbzI, smxI) = (int(folder), 0)
        simulation = Simulation(cbzI, smxI, [])
        fileName = "%s/%s" % (folder, 'output.log')
        if not os.path.isfile(fileName): continue
        if os.stat(fileName).st_size == 0: continue
        molecules = {'cbz': Pharmaceutical(), 'smx': Pharmaceutical(), 'h2o': Pharmaceutical(), 'Na': Pharmaceutical()}
        with open(fileName, 'r') as filin:
            for line in filin:
                if line == 'Number of molecules:\n':
                    break
            for line in filin:
                elem = line.split()
                if len(elem) == 3 and elem[0] == 'Component':
                    compName = elem[2][1:-1]
                    for line in filin:
                        elem = line.split()
                        if len(elem) >= 6 and elem[2] == 'absolute' and elem[3] == '[molecules/unit':
                            thisPharma = Pharmaceutical(float(elem[5]), float(elem[7]))
                        if len(elem) >= 6 and elem[2] == 'absolute' and elem[3] == '[milligram/gram':
                            thisPharma.mgg = float(elem[5])
                            thisPharma.mgg_std = float(elem[7])
                            molecules[compName] = thisPharma
                            break
        simulation.components = molecules
        simulations.append(simulation)
    #Output
    with open(outputFile, 'w') as filout:
        filout.write("%9s %9s %-18s %-18s %-18s %-18s %-18s %-18s %-18s %-18s\n" % ('CBZ_[ppb]', 'SMX_[ppb]', '  CBZ_[mg/g]', '  CBZ_[mol]', '  SMX_[mg/g]', '  SMX_[mol]', '  H2O_[mg/g]', '  H2O_[mol]', '  Na_[mg/g]', '  Na_[mol]'))
        for simulation in simulations:
            filout.write("%9d %9d %7.3f +- %-7.3f %7.3f +- %-7.3f %7.3f +- %-7.3f %7.3f +- %-7.3f %7.3f +- %-7.3f %7.3f +- %-7.3f %7.3f +- %-7.3f %7.3f +- %-7.3f\n" % (simulation.cbzI, simulation.smxI, simulation.components['cbz'].mgg, simulation.components['cbz'].mgg_std, simulation.components['cbz'].mol, simulation.components['cbz'].mol_std, simulation.components['smx'].mgg, simulation.components['smx'].mgg_std, simulation.components['smx'].mol, simulation.components['smx'].mol_std, simulation.components['h2o'].mgg, simulation.components['h2o'].mgg_std, simulation.components['h2o'].mol, simulation.components['h2o'].mol_std, simulation.components['Na'].mgg, simulation.components['Na'].mgg_std, simulation.components['Na'].mol, simulation.components['Na'].mol_std))
            
