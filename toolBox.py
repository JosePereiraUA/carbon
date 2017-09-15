from __future__ import print_function
import numpy as np
import random
import math
import copy
import sys
import os

def div():
    consoleprint("-------------------------")

def consoleprint(*args, **kwargs):
    print(*args, file=sys.stderr, end="\n", **kwargs)

def progress(*args, **kwargs):
    print(*args, file=sys.stderr, end="\r", **kwargs)

def mufla(nx, ny):
    r, n = 1.4, 1
    x1 = 0.0
    y1 = r
    x2 = r * math.cos(math.pi / 6.0)
    y2 = r * math.sin(math.pi / 6.0)
    ax, ay = 2 * x2, 0.0
    bx, by = x2, 1.5 * r

    lamina = Lamina()
    for j in range(ny+1):
        for i in range(-j / 2, nx - j / 2):
            dx = i * ax + j * bx
            dy = i * ay + j * by
            atom1 = Atom(0, "C", [x1 + dx, y1 + dy, 0.0], 'ca')
            atom2 = Atom(0, "C", [x2 + dx, y2 + dy, 0.0], 'ca')
            lamina.addAtom(atom1)
            lamina.addAtom(atom2)
    
    if (ny) % 2 != 0:
        del lamina.atoms[(nx * 2 * (ny) + ((ny + 1 / 2) - 1) )]
 
    for j in range(2, ny + 1, 2):
        dx = (-j / 2 - 1) * ax + j * bx
        dy = (-j / 2 - 1) * ay + j * by
        atom1 = Atom(0, "C", [x2 + dx, y2 + dy, 0.0], 'ca')
        lamina.addAtom(atom1)

    for j in range(0, ny, 2):
        dx = (nx - j / 2) * ax + j * bx
        dy = (nx - j / 2) * ay + j * by
        atom1 = Atom(0, "C", [x1 + dx, y1 + dy, 0.0], 'ca')
        lamina.addAtom(atom1)

    for n, atom in enumerate(lamina.atoms, start=1):
        atom.n = n

    return lamina

def axisAngleRotation(w, t, v):
    '''
    Using Rodrigues rotation formula:
    http://en.wikipedia.org/wiki/Axis-angle_representation
      v' = cos(t)v + (w x v)sin(t) + w(w.v)(1-cos(t))

     @param: w - axis
     @param: t - theta
     @param: v - Nx3 matrix to be rotated
    '''
    m = 0
    try:
        n, m = v.shape
    except:
        n, = v.shape

    wn = w / np.sqrt(np.dot(w, w))  # normalize axis
    ct = np.cos(t)

    # case 2D array
    if m:
        return v * ct + np.cross(wn, v) * np.sin(t) + wn * (1 - ct) * np.reshape(np.repeat((wn * v).sum(1), 3), (-1, 3))

    # case 1D array
    return v * ct + np.cross(wn, v) * np.sin(t) + wn * (1 - ct) * (wn * v).sum()

def randomUnitVector(v, w):
    """ v,w - two uniform random numbers from the interval [0,1[
    returns: unit 3-vector randomly distributed
    """
    ctheta = 1.0 - 2.0 * v
    phi = np.pi * 2.0 * w
    stheta = np.sin(np.arccos(ctheta))
    sphi = np.sin(phi)
    cphi = np.cos(phi)
    return np.array([stheta * cphi, stheta * sphi, ctheta])

def copyDic(dic):
    out = {}
    for key, i in dic.items():
        out[key] = i
    return (out)

def recoverParams(library, n):
    params = {}
    tot = 0
    for key, array in library.items():
        this = int(math.floor(array[n]))
        params[key] = this
        tot += this
    dif = 100 - tot
    if not dif == 0 and len(params) > 3:
        for v in range(abs(dif)):
            key = random.sample(params, 1)
            if dif > 0: params[key[0]] += 1
            else: params[key[0]] -= 1
    
    return params

#For itp preparation
def writeSource(source, title):
    with open(source, 'ab') as filout:
        filout.write("source leaprc.gaff\n")
        filout.write("frcmod = loadamberparams ~/carbons/1_models/3_masters/master.frcmod\n")
        filout.write("molecule = loadmol2 %s.mol2\n" % (title))
        filout.write("saveamberparm molecule %s.top %s.crd\n" %(title, title))
        filout.write("quit")

def writeITP(title):
    saved = []
    if int(title) > 99: fix(title)
    top = "R%s_GMX.top" % (title)
    itp = "%s.itp" % (title)
    with open(top, 'r') as fin:
        for line in fin:
            elem = line.split()
            if len(elem) > 0 and elem[1] == 'moleculetype':
                saved.append(line)
                break
        for line in fin:
            elem = line.split()
            if len(elem) > 0 and elem[1] == 'system': break
            saved.append(line)
    with open(itp, 'ab') as filout:
        for line in saved:
            filout.write(line)
    

def fix(title):
    os.system("mv %s_GMX.top R%s_GMX.top" % (title, title))
    filin = "R%s_GMX.top" % (title)
    new_data = []
    with open(filin, 'r') as fin:
        data = fin.readlines()
    for line in data:
        elem = line.split()
        if len(elem) == 2 and elem[0] == title:
            new_data.append("%4s              3\n" % ("R"+title))
        else: new_data.append(line)
    with open(filin, 'w') as filout:
        filout.writelines(new_data)

class Residue: #-----------------------------------------------------------
    def __init__(self):
        self.laminas = []
        self.list = []

    def addLamina(self, lamina, nz):
        for n in range(nz):
            if n == 0:
                lamina_new = lamina
            else: 
                vi = ([1.4 * math.cos(math.pi / 6.0) * n, (1.4 / 2) * n, 3.4 * n])
                lamina_new = lamina.move([vi[0], vi[1], vi[2]], n)
            self.laminas.append(lamina_new)
    
    def checkChrg(self):
        total, totalAtoms = 0, []
        for lamina in self.laminas:
            for atom in lamina.atoms:
                total += atom.charge
        if not total == 0.0:
            for lamina in self.laminas:
                for atom in lamina.atoms:
                    totalAtoms.append(atom)
            theChosenOne = random.sample(totalAtoms, 1)[0]
            theChosenOne.charge -= total
            total = 0
            for lamina in self.laminas:
                for atom in lamina.atoms:
                    total += atom.charge
            if total < 0.0000001: total = 0.0
        if total != 0: consoleprint("Total charge:", total)
        
        
    def writePDB(self, fout):
        self.defineBonds()
        atoms = []
        for lamina in self.laminas:
            for atom in lamina.atoms:
                atoms.append(atom)
        atoms = sorted(atoms, key= lambda atom: atom.n)
        with open(fout, 'ab') as filout:
            form = 'ATOM %6d  %-3s UNK   %3d    %8.3f%8.3f%8.3f ; %-3s %-3.4f\n'
            for atom in atoms:
                filout.write(form % (atom.n, atom.symbol, 1, atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.atmtype, atom.charge))
            filout.write('TER\n')
            for n, lamina in enumerate(self.laminas):
                for i, atom in enumerate(lamina.atoms):
                    filout.write("CONECT%5d" % atom.n)
                    filout.write(''.join('%5d' % (j, ) for j in atom.connects))
                    filout.write("\n")
    
    def writeMol2(self, fout, title):
        resN = "R"+title
        self.defineBonds()
        atoms = []
        for lamina in self.laminas:
            for atom in lamina.atoms:
                atoms.append(atom)
        atoms = sorted(atoms, key= lambda atom: atom.n) #mol2 format requires order
        bonds = self.mol2Bonds()
        with open(fout, 'ab') as filout:
            filout.write("# created with carbFunct2.py by Ze Manel\n@<TRIPOS>MOLECULE\n%s\n%d %d %d\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n" % (resN, self.countAtoms(), len(bonds), 1))
            form  = '%-6d %-4s %6.3f %6.3f %6.3f %-9s %2d %-8s %9.6f\n'
            for atom in atoms:
                filout.write(form % (atom.n, atom.symbol, atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.atmtype, 1, 'CRV1', atom.charge))
            filout.write("@<TRIPOS>BOND\n")
            form = '%d %d %d %d\n'
            for n, bond in enumerate(bonds):
                filout.write(form % (n, bond[0], bond[1], 1))
            filout.write("@<TRIPOS>SUBSTRUCTURE\n%-3d %-7s %-3d %-5s %3d %-5s %s\n" % (1, resN, 1, 'GROUP', 1, '****', resN))
     
    def mol2Bonds(self):
        bonds = []
        for lamina in self.laminas:
            for atom1 in lamina.atoms:
                for atom2_ID in atom1.connects:
                    this = [atom1.n, atom2_ID]
                    if len(bonds) == 0:
                        bonds.append(this)
                    else:
                        exists = False
                        for bond in bonds:
                            if len(set(this)^set(bond)) == 0:
                                exists = True
                        if exists == False:
                            bonds.append(this)
        return bonds
                
    def writeDat(self, fout, sizes):
        tot = self.countAtoms()
        count = self.countFcns()
        elem = self.elementalAnalysis()
        with open(fout, 'ab') as filout:
            filout.write("\n[ Residue details ]\n")
            filout.write("charge : %d\n" % (count['coo'] * -1))
            for key in sizes:
                filout.write("%s : %d\n" % (key, sizes[key]))
            filout.write("\n[ Elemental Analysis (mol)  ]\n")
            for key in elem:
                per = (float(elem[key])/tot)*100
                filout.write("%s : %3.2f%% ( %3d / %3d )\n" % (key, per, elem[key], tot))
            filout.write("\n[ Elemental Analysis (mass) ]\n")
            w, totalW = {'C': 12.0107, 'O': 15.9994, 'H': 1.00794}, 0
            for key in elem: totalW += elem[key]*w[key]
            for key in elem:
                per = ((float(elem[key])*w[key])/totalW)*100
                filout.write("%s : %3.2f%%\n" % (key, per))
            filout.write("\n[ Functional content (mol)  ]\n")
            tot = self.countAtoms('C')
            for key in count:
                per = (float(count[key])/tot)*100
                if not count[key] == 0:
                    filout.write("%-4s : %3.2f%% ( %3d / %3d )\n" % (key, per, count[key], tot))
            filout.write("\n---\n")
    
    def elementalAnalysis(self):
        out = {'C': 0, 'O': 0, 'H': 0}
        for lamina in self.laminas:
            for atom in lamina.atoms:
                out[atom.symbol] += 1
        return out
                
    def defineBonds(self):
        atoms = []
        for lamina in self.laminas:
            for atom in lamina.atoms:
                atoms.append(atom)
        old_atoms = copy.deepcopy(atoms)
        for n, atom in enumerate(atoms, start=1):
            atom.n = n
        table = {}
        for n, atom in enumerate(old_atoms):
            table[atom.n] = atoms[n].n
        self.updateList(table)
        self.updateConnectivities()

    def updateList(self, table):
        for item in self.list:
            item[0] = table[item[0]]
            out = []
            for con in item[2]:
                out.append(table[con])
            item[2] = out
            

    def updateConnectivities(self, threshold=1.55):
        #Clear old connectivities
        for lamina in self.laminas:
            for atom in lamina.atoms:
                del atom.connects[:]
        #Add carbon connects
        dSq = threshold**2
        sizes = self.getSizes()
        for lamina in self.laminas:
            for i, at1 in enumerate(lamina.atoms[:-1]):
                for at2 in lamina.atoms[i+1:]:
                    if at1.symbol == 'C' and at2.symbol == 'C':
                        d = 0
                        for p in range(len(at1.xyz)):
                            d += math.pow(at1.xyz[p]-at2.xyz[p], 2)
                        if d < dSq:
                            at1.connects.append(at2.n)
                            at2.connects.append(at1.n)
        #Add functionalization connects
        for lamina in self.laminas:
            for at1 in lamina.atoms:
                for item in self.list:
                    if at1.n == item[0]:
                        if item[1] == 'ror':
                            for con1 in item[2]:
                                at1.connects.append(con1)
                        else:
                            at1.connects.append(item[2][0])
                            for i, con1 in enumerate(item[2]):
                                for j, con2 in enumerate(item[2]):
                                    if self.getAtom(con1).symbol != self.getAtom(con2).symbol:
                                        if item[1] == 'cooh':
                                            if i == j+1:
                                                self.getAtom(con1).connects.append(con2)
                                            self.getAtom(item[2][0]).connects.append(item[2][2])
                                        else:
                                            self.getAtom(con1).connects.append(con2)
        for lamina in self.laminas:
                for atom in lamina.atoms:
                    atom.connects = sorted(atom.connects)
       
    def distribute(self):
        lanes = {'bot': [], 'mid': [], 'top': [], 'lat': []}
        nz = len(self.laminas)
        for n, lamina in enumerate(self.laminas):
            for atom in lamina.atoms:
                if len(atom.connects) == 2:
                    #atom.symbol = 'S'
                    lanes['lat'].append(atom.n)
                if n == 0 and len(atom.connects) == 3:
                    #atom.symbol = 'N'
                    lanes['bot'].append(atom.n)
                if nz > 1 and n == (nz-1) and len(atom.connects) == 3:
                    #atom.symbol = 'P'
                    lanes['top'].append(atom.n)
                if nz > 1 and not n == 0 and not n == (nz-1) and len(atom.connects) == 3:
                    #atom.symbol = 'He'
                    lanes['mid'].append(atom.n)
        if nz == 1:
            lanes['top'] = random.sample(lanes['bot'], random.randint(0, len(lanes['bot'])))
            lanes['bot'] = [x for x in lanes['bot'] if x not in lanes['top']]
        
        return lanes
    
    def listAtoms(self):
        out = []
        for lamina in self.laminas:
            for atom in lamina.atoms:
                out.append(atom.n)
        return out
    
    def setRORs(self, lanes, objective):
        possibleRORs = {2: 2, 3: 4}
        total = self.countAtoms('C')
        percentage, oxygen = 0, 0
        while percentage < objective:
            verified = False
            while not verified:
                #Pick random lane
                lane = random.choice(lanes.items())
                if len(lane[1]) == 0: continue
                #Pick random atom from that lane
                ID = random.sample(lane[1], len(lane[1]))[0]
                atom0 = self.getAtom(ID)
                lane = lane[0]
                #Pick random ROR type
                if not lane == "lat":
                    ROR = random.choice(possibleRORs.items())
                    loss = ROR[1]
                    ROR = ROR[0]
                if lane == "lat": ROR, loss = 1, 0
                #Verify if placement of this ROR is possible
                (verified, selectedAtom) = self.verifyROR(atom0, ROR, lanes[lane])
            #Add ROR to selectedAtom and remove neighbouring carbons
            atom0.symbol = 'O'
            lanes[lane] = self.addROR(atom0, ROR, lanes[lane], selectedAtom)
            #Calculate current oxygen percentage
            total -= loss
            oxygen += ROR
            percentage = (2 * float(oxygen) / total) * 100 #!
    
    def verifyROR(self, atom0, ROR, lane):
        if ROR == 1:
            for connect in atom0.connects:
                if self.getAtom(connect).symbol == 'O':
                    return (False, atom0)
            return (True, atom0)
        if ROR == 2:
            i = 0
            for connect in atom0.connects:
                if self.getAtom(connect).symbol == 'O':
                    return (False, atom0)
            while i < len(atom0.connects):
                atom2 = self.getAtom(atom0.connects[i])
                if not atom2.n in lane:
                    i += 1
                    continue
                isOk = True
                for connect in atom2.connects:
                    if self.getAtom(connect).symbol == 'O':
                        i += 1
                        isOk = False
                        break
                if isOk: 
                    return (True, atom2)
            return (False, atom0)
        if ROR == 3:
            atom2_IDs = set(atom0.connects)
            for ID in atom2_IDs:
                atom2 = self.getAtom(ID)
                if len(atom2.connects) <= 2:
                    return (False, atom0)
                for connect in atom2.connects:
                    if self.getAtom(connect).symbol == 'O':
                        return (False, atom0)
            return (True, atom0)
    
    def getAtom(self, n):
        for lamina in self.laminas:
            for atom in lamina.atoms:
                if atom.n == n: return atom        
    
    def addROR(self, atom0, ROR, lane, selectedAtom):
        if ROR == 1:
            atom0.atmtype = 'os'
            for _ in range(2): self.list.append([atom0.n, 'ror', atom0.connects])
            return sorted(set(lane) - set(atom0.connects + [atom0.n]))
        if ROR == 2:
            #Set symbols
            selectedAtom.symbol = 'O'
            #Set atomtypes
            atom0.atmtype, selectedAtom.atmtype = 'os', 'os'
            #Append functionalizations to residue list
            for _ in range(2): self.list.append([atom0.n, 'ror', atom0.connects])
            for _ in range(2): self.list.append([selectedAtom.n, 'ror', selectedAtom.connects])
            #Remove unwanted connections
            toRemove = set(atom0.connects + [atom0.n] + selectedAtom.connects + [selectedAtom.n])
            if selectedAtom.n in atom0.connects: atom0.connects.remove(selectedAtom.n)
            if atom0.n in selectedAtom.connects: selectedAtom.connects.remove(atom0.n)
            #Set correct charges
            self.spreadCharge(FunctGroups['ror'+str(ROR)], atom0)
            self.spreadCharge(FunctGroups['ror'+str(ROR)], selectedAtom)
            return sorted(set(lane) - toRemove)
        if ROR == 3:
            toRemove = []
            for atom_ID in atom0.connects:
                atom = self.getAtom(atom_ID)
                #Set symbols
                atom.symbol = 'O'
                #Set atomtypes
                atom.atmtype = 'os'
                #Append functionalizations to residue list
                for _ in range(2): self.list.append([atom.n, 'ror', atom.connects])
                #Remove unwanted connections
                toRemove.append(atom.connects)
                toRemove.append([atom.n])
            toRemove = [x for x in sum(toRemove, []) if not x == atom0.n] + [atom0.n]
            for atom_ID in atom0.connects:
                atom = self.getAtom(atom_ID)
                atom.connects.remove(atom0.n)
                self.spreadCharge(FunctGroups['ror'+str(ROR)], atom)
            self.removeAtom(atom0)
            return sorted(set(lane) - set(toRemove))
    
    def spreadCharge(self, charge, atom):
        atom.charge = charge
        toSpread = (charge * -1)/len(atom.connects)
        for atom_ID in atom.connects:
            atom2 = self.getAtom(atom_ID)
            atom2.charge += toSpread
    
    def removeAtom(self, atomDel):
        for lamina in self.laminas:
            for n, atom in enumerate(lamina.atoms):
                if atomDel.n == atom.n:
                    lamina.atoms.pop(n)
    
    def countAtoms(self, only='NaN'):
        i = 0
        for lamina in self.laminas:
            for atom in lamina.atoms:
                if only == 'NaN': i += 1
                elif atom.symbol == only: i+= 1
        return i
    
    def countFcns(self):
        out = {'sp2': 0, 'ch': 0, 'ror': 0, 'co': 0, 'coo': 0, 'cooh': 0}
        for reg in self.list:
            out[reg[1]] += 1
        tot = self.countAtoms('C')
        out['sp2'] = tot - len(self.list)
        return out
    
    def getLastAtom(self):
        i = 0
        for lamina in self.laminas:
            for atom in lamina.atoms:
                if atom.n > i: i = atom.n
        return i
               
    def getSizes(self):
        sizes = [0]
        for n, lamina in enumerate(self.laminas):
            sizes.append(len(lamina.atoms) + sizes[n])
        return sizes

    def setPos(self, count, params, lanes):
        #Set number of each funct
        total = 0
        for key in params:
            if params[key] <= 0: i = 0
            else: i = int(math.floor(count*(float(params[key])/100)))
            params[key] = i
            total += i
        excess = count - total
        params['sp2'] += excess
        #Randomize layers
        bot = random.sample(lanes['bot'], len(lanes['bot']))
        mid = random.sample(lanes['mid'], len(lanes['mid']))
        top = random.sample(lanes['top'], len(lanes['top']))
        lat = random.sample(lanes['lat'], len(lanes['lat']))
        #Set positions of each funct
        kompensam = 0
        #CO
        co = []
        for i in range(params['co']):
            if len(lat) <= 0: break
            co.append(lat[0])
            if len(lat) > 0: lat.remove(lat[0])
            if len(lat) > 0: lat.remove(lat[0])
        params['co'] = co
        '''kompensam += params['co']
        co = lat[0:params['co']]
        params['co'] = co
        lat = [x for x in lat if x not in co]'''
        #COO
        botTopLat = bot + top + lat
        botTopLat = random.sample(botTopLat, len(botTopLat))
        coo = []
        for i in range(params['coo']):
            coo.append(botTopLat[0])
            if botTopLat[0] in top or botTopLat[0] in bot: kompensam += 1
            botTopLat = self.removeClose(botTopLat)
        params['coo'] = coo
        #CH
        lat = [x for x in lat if x not in coo]
        if kompensam == 0: params['ch'] = lat
        else: params['ch'] = lat[:-kompensam]
        
        return params    
    
    def removeClose(self, array):
        for connect in self.getAtom(array[0]).connects:
            if connect in array: array.remove(connect)
        array.remove(array[0])
        return array
        
    def functionalize(self, params, lanes):
        sizes = self.getSizes()
        for key, array in params.items():
            if key != 'sp2' and key != 'ror':
                for f in array:
                    for n, lamina in enumerate(self.laminas):
                        for atom in lamina.atoms:
                            if atom.n == f:
                                lat, bot = False, False
                                if f in lanes['lat']:
                                    lat = True
                                if not lat:
                                    if f in lanes['bot']:
                                        bot = True
                                #if atom.symbol != 'C' : break
                                tot = self.getLastAtom()
                                lamina.addFunct(key, atom, sizes[n], lat, bot, tot)
                                self.list.append(lamina.list[len(lamina.list)-1])
                                break
                    
    
class Lamina: #-------------------------------------------------------------------

    def __init__(self):
        self.atoms = []
        self.list = []
        
    def addAtom(self, atom):
        self.atoms.append(atom)
        
    def move(self, v, n=0):
        lam = Lamina()
        n = (len(self.atoms) * n)
        for atom in self.atoms:
            lam.addAtom(atom.move(v, n))
        return lam
        
    def extend(self, residue):
        for atom in residue.atoms:
            self.addAtom(atom)
    
    def getBot(self):
        out = Lamina()
        for atom in self.atoms:
            tmpatm = Atom(atom.n, atom.symbol, [atom.xyz[0], atom.xyz[1], -atom.xyz[2]], atom.atmtype, atom.charge)
            out.addAtom(tmpatm)
        return out
    
    def getLat(self, at0, fcn, n):
        for atom in self.atoms:
            if atom.n == at0.connects[0]: at1 = atom
            if atom.n == at0.connects[1]: at2 = atom   
        v1 = np.asarray(at1.xyz) - np.asarray(at0.xyz)
        v2 = np.asarray(at2.xyz) - np.asarray(at0.xyz)
        v3 = v1 + v2
        w = np.cross([0,0,1], v3)
        for n, fcnatom in enumerate(fcn.atoms):
            new_xyz = axisAngleRotation(w, -(np.pi / 2), np.asarray(fcnatom.xyz))
            fcn.atoms[n] = Atom(fcn.atoms[n].n, fcn.atoms[n].symbol, new_xyz, fcn.atoms[n].atmtype, fcn.atoms[n].charge)
        return fcn
    
    def upright(self, lamina, step = 0.05):
        tmplam = Lamina()
        xyz0, rot = [], 0
        conflict = True
        for atom in self.atoms:
            xyz0.append(atom.xyz[0])
            xyz0.append(atom.xyz[1])
            xyz0.append(atom.xyz[2])
        xyz0 = np.array(xyz0)
        xyz0 = xyz0.reshape(len(self.atoms), 3)
        axis = [0.0, 0.0, 1.0]
        xyz2 = axisAngleRotation(axis, 0, xyz0)
        for n, atom in enumerate(self.atoms):
            tmpatm = Atom(atom.n, atom.symbol, xyz2[n], atom.atmtype, atom.charge)
            tmplam.addAtom(tmpatm)
        return tmplam
    
    def addFunct(self, key, at0, n, lat, bot, tot):
        fcn = FunctGroups[key]
        for atom in fcn.atoms:
            tot += 1
            atom.n = tot
        fcn = fcn.upright(self)
        if bot: fcn = fcn.getBot()
        elif lat: fcn = self.getLat(at0, fcn, n)
        new = fcn.move(at0.xyz)
        #Set correct atomtypes
        if lat: at0.atmtype = 'ca'
        else: at0.atmtype = 'c3'
        if key == 'co': at0.atmtype = 'c'
        #Set correct charges
        at0.charge += FunctGroups[key+'_Ccharge']
        split = FunctGroups[key+'_Fcharge'] / len(at0.connects)
        for atom in self.atoms:
            for index in at0.connects:
                if atom.n == index:
                    atom.charge += split
        #Set correct connectivities
        '''new.atoms[0].connects.append(at0.n)
        for atom1 in new.atoms[1:end]:
                new.atoms[0].connects.append(atom1.n)
        if fcn == 'cooh':
            new.atoms[len(new.atoms)].connects = [new.atoms[len(new.atoms)-1]]'''
        self.extend(new)
        new_atoms = []
        for atom in new.atoms:
            new_atoms.append(atom.n)
        self.list.append([at0.n, key, new_atoms])
     
class Atom: #-------------------------------------------------------q
    
    def __init__(self, n="NaN", symbol='NaN', xyz=[0.,0.,0.], atmtype='NaN', charge = 0.000):
        self.n = n
        self.symbol = symbol
        self.xyz = xyz;
        self.connects = []
        self.atmtype = atmtype
        self.charge = charge
        
    def move(self, v, n):
        new = ([self.xyz[0]+v[0], self.xyz[1]+v[1], self.xyz[2]+v[2]])
        return Atom(self.n+n, self.symbol, new, self.atmtype, self.charge)
       
       
FunctGroups = {}

#------CH------------------------------------------------------------
FunctGroups['ch'] = Lamina()
FunctGroups['ch'].addAtom(Atom('NaN', 'H', [0.0, 0.0, 0.8], 'ha', -0.16))
FunctGroups['ch_Ccharge'] = 0.16
FunctGroups['ch_Fcharge'] = 0

#------2RING---------------------------------------------------------
FunctGroups['ror2'] = -0.25

#------3RING---------------------------------------------------------
FunctGroups['ror3'] = -0.29

#------COO-----------------------------------------------------------
FunctGroups['coo'] = Lamina()
FunctGroups['coo'].addAtom(Atom('NaN', 'C',[ 0.00,  0.00, 1.49],'c2', 0.83))
FunctGroups['coo'].addAtom(Atom('NaN', 'O',[ 1.04,  0.17, 2.10],'o', -0.76))
FunctGroups['coo'].addAtom(Atom('NaN', 'O',[-1.15, -0.19, 2.16],'o', -0.76))
FunctGroups['coo_Ccharge'] = 0.03
FunctGroups['coo_Fcharge'] = -0.34

#------COOH----------------------------------------------------------
FunctGroups['cooh'] = Lamina()
FunctGroups['cooh'].addAtom(Atom('NaN', 'C',[ 0.00,  0.00, 1.68],'c2', 0.85))
FunctGroups['cooh'].addAtom(Atom('NaN', 'O',[ 0.47, -0.90, 2.41],'o', -0.63))
FunctGroups['cooh'].addAtom(Atom('NaN', 'O',[-0.58,  1.14, 2.19],'oh', -0.63))
FunctGroups['cooh'].addAtom(Atom('NaN', 'H',[-0.51,  1.01, 3.15],'ho', 0.42))
FunctGroups['cooh_Ccharge'] = 0.07
FunctGroups['cooh_Fcharge'] = -0.08

#------CO------------------------------------------------------------
FunctGroups['co'] = Lamina()
FunctGroups['co'].addAtom(Atom('NaN', 'O', [0.0, 0.0, 1.43], 'o', -0.56))
FunctGroups['co_Ccharge'] = 0.54
FunctGroups['co_Fcharge'] = 0.02

