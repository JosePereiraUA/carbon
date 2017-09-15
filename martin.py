#!/usr/bin/python

from toolBox import *
import time

def help():
    consoleprint("\n  > \033[92mMartin\033[0m: Draws activated carbon residue from scratch, using input data from experiemental essays\n  > Last updated: 26 May 2017 by Ze Manel\n")
    consoleprint("  \033[93mUSAGE\033[0m:")
    consoleprint("  martin.py [-\033[92mm\033[0m] [-\033[92ma\033[0m] [-\033[92mr\033[0m] [-\033[92mt\033[0m] [-\033[92mi\033[0m] [-\033[92mn\033[0m] [-e]\n\n  \033[93mPARAMETERS\033[0m:\n  -h (--help)        : Displays this help screen\n  -\033[92mm\033[0m (--mol2)        : Write complementary mol2 file\n  -\033[92ma\033[0m (--all)         : Write complementary mol2 file, itp file and itp_rest file\n  -\033[92mr\033[0m (--rest)        : Write rest_itp file\n  -\033[92mt\033[0m (--title) [str] : Manually select output files names\n  -\033[92mi\033[0m (--itp)         : Write itp file (requires -m option for valid mol2 file)\n  -\033[92mn\033[0m (--num)   [int] : Perform the whole process a number of times, separate folders of each residue\n  -e (--extra)       : Output dat file with elem. analysis and functionalization content information\n  -v (--verb)        : Run script in verbose mode\n")
    consoleprint("  \033[93mOUTPUT\033[0m:")
    consoleprint("  + title.pdb\n [+ title.mol2]\n [+ title.itp]\n [+ title_rest.itp]\n [+ title.dat]\n")
    div()

if __name__ == '__main__':

    #Benchmark process
    startTime = time.clock()
    
    div()

    #Input parsing
    m, a, r, i, e = False, False, False, False, False
    t = "NaN"
    v = ">/dev/null 2>&1"
    n = 1
    for ind, arg in enumerate(sys.argv):
        if arg == '-h' or arg == '--help':
            help()
            exit(1)
        elif arg == '-m' or arg == '--mol2': m = True
        elif arg == '-a' or arg == '--all': a = True
        elif arg == '-r' or arg == '--rest': r = True
        elif arg == '-t' or arg == '--title': t = sys.argv[ind+1]
        elif arg == '-i' or arg == '--itp': i = True
        elif arg == '-n' or arg == '--num': n = int(sys.argv[ind+1])
        elif arg == '-e' or arg == '--extra': e = True
        elif arg == '-v' or arg == '--verb': v = ""
    pdbF = ("/home/jpereira/carbons/1_models/3_masters/master_RL.dat")
    
    #Input verification
    if a:
        m, i, r, e = True, True, True, True
    if i and not m:
        consoleprint("No -m option detected. Turning --mol2 on for correct .mol2 and .itp definition.")
        m = True
        div()
    
    #Obtain parameter from input master file
    params, total, scale = {'sp2': 0, 'ch': 0, 'ror': 0, 'co': 0, 'coo': 0}, 0, 0
    sizes = {'nx': (0, 0), 'ny': (0, 0), 'nz': (0, 0)}
    with open(pdbF, 'r') as fin:
        for line in fin:
            elem = line.split()
            if elem[0] == 'scale':
                scale = np.sqrt(float(elem[1]))
                continue
            if elem[0] == 'charge':
                charge = elem[1]
                continue
            for key in params:
                if elem[0] == key:
                    params[key] = float(elem[1])
                    total += float(elem[1])
                    continue
            for key in sizes:
                if elem[0] == key:
                    sizes[key] = (int(elem[1]), int(elem[2]))
                    continue

                    
    #Check integrity of obtained parameters
    if total != 100:
        consoleprint("Input values don't add to 100%. Please check master.dat.")
        div()
        exit(1)
    for key in sizes:
        if sizes[key] == 0:
            consoleprint("No value found for %s. Please check master.dat." % key)
            div()
            exit(1)

    #Set library for multiple files
    libraryF = {'sp2': 0, 'ch': 0, 'ror': 0, 'co': 0, 'coo': 0}
    libraryS = {'nx': 0, 'ny': 0, 'nz': 0}
    for key in libraryF:
        libraryF[key] = np.random.normal(params[key], scale, n).tolist()
    for key in libraryS:
        libraryS[key] = [random.randint(sizes[key][0], sizes[key][1]) for p in range(n)]

    #Delete all folders in current working directory
    if n > 1:
        for index in range(n):
            
            title = str(index + 1)
            if not t == "NaN": title = t + title
            os.system("rm -r %s >/dev/null 2>&1" % (title))
    
    #Display header
    consoleprint("Initiating library assembly.\nUsing experimental data stored in: %s\nCurrent working directory: %s\n" % (pdbF, os.getcwd()))
    
    #Build and write pdb (and mol2) files
    for index in range(n):
    
        #Show progress
        progress("Assembling library. Progress:", index+1, "/", n, "(%2.f%%)" % ((float(index+1)/n)*100))
    
        #Title atribution
        if n == 1 and not t == "NaN": title = t
        else: 
            title = str(index + 1)
            if not t == "NaN": title = t + title
        if not os.path.isdir("%s" % (title)) and n > 1: os.system("mkdir %s" % (title))
        
        #Clean previous files and set output files names (in case of n = 1)
        foutPDB = "%s.pdb" % (title)
        with open(foutPDB, 'w'): pass
        if m:
            foutMol2 = "%s.mol2" % (title)
            with open(foutMol2, 'w'): pass
        if i:
            foutITP = "%s.itp" % (title)
            source = "source"
            with open(source, 'w'): pass
            with open(foutITP, 'w'): pass
        if e:
            foutDAT = "%s.dat" % (title)
            with open(foutDAT, 'w'): pass
    
        #Create residue
        params = recoverParams(libraryF, index)
        sizes = recoverParams(libraryS, index)
        res = Residue()
        lam = mufla(sizes['nx'], sizes['ny'])
        #lam = mufla(4, 4) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        res.addLamina(lam, sizes['nz'])
        #res.addLamina(lam, 5) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        res.updateConnectivities()
        lanes = res.distribute()
        lanes_saved = copyDic(lanes)
        res.setRORs(lanes, params['ror'])
        #Fix lanes because of removed atoms and added oxygens
        fullAtoms = res.listAtoms()
        for key, lane in lanes_saved.items():
            for atom_ID in sorted(lane):
                if not atom_ID in fullAtoms or res.getAtom(atom_ID).symbol == 'O':
                    lanes_saved[key].remove(atom_ID)
        #Functionalize
        count = res.countAtoms()
        params = res.setPos(count, params, lanes_saved)
        #Check PCZ
        if charge == '+':
            params['cooh'] = params['coo']
            del params['coo']
        res.functionalize(params, lanes_saved)
        #res.checkChrg()
        res.writePDB(foutPDB)
        if m: res.writeMol2(foutMol2, title)
        
        #Write dat file
        if e: res.writeDat(foutDAT, sizes)
    
        #Write itp file
        if i:
            writeSource(source, title)
            os.system("tleap -f %s %s" % (source, v))
            os.system("acpype -p %s.top -x %s.crd %s" % (title, title, v))
            writeITP(title)
            
        #Write rest_itp file
        if r:
            os.system("yes 0 2>/dev/null 2>&1 | gmx genrestr -f %s.pdb -o %s_rest.itp %s" % (title, title, v))
        
        #Save required files, delete the rest
        if n > 1: os.system("mv -t %s/ *.itp *.mol2 *.pdb *.dat >/dev/null 2>&1 && rm * >/dev/null 2>&1" % (title))
        else:
            os.system("mkdir temp && mv -t temp/ *.itp *.mol2 *.pdb *.dat *.top >/dev/null 2>&1 && rm * >/dev/null 2>&1")
            os.system("mv temp/* . && rm -r temp")
    
    #Display and exit
    consoleprint("\n\nLibrary assembly sucessefull. %d out of %d requested residues built in %f s" % (n, n, time.clock() - startTime))
    div()
    
