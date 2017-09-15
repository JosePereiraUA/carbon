from __future__ import print_function
import sys
import os

def consoleprint(*args, **kwargs):
    print(*args, file=sys.stderr, end="\n", **kwargs)

if __name__ == '__main__':
    
    #Update run.sh with simulation name
    run = '/home/jose.manuel.pereira@UA.PT/run.sh'
    if len(sys.argv) == 2: title = sys.argv[1]
    else: title = 'default'
    data = []
    with open(run, 'r') as filin:
        for line in filin:
            data.append(line[:-1])
    data[1] = '#SBATCH --job-name="%s"' % (title)
    with open(run, 'w') as filout:
        for line in data:
            filout.write(line + "\n")
    
    #Copy run.sh to simulation folders and start RASPA
    main = os.getcwd()
    l = [x for x in os.listdir(main)]
    for folder in l:
        os.system("cp %s %s" % (run, folder))
        os.chdir(folder)
        os.system("sbatch run.sh")
        os.chdir(main)
