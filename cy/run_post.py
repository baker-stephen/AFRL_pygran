import sys

import cy.pdn as pdn
import cy.ultp as ultp
import cy.plotme_sin as plotme_sin

sys.path.insert(1,'../')

from param_defn import PD

def go(params: PD, name_mod = ""):
        
    params.out_dir = params.out_dir[:params.out_dir.rfind('/')+1]

    print("new write directory:", params.out_dir)

    # sys.stdout = open(params.out_dir + 'stdout.txt', 'w')
    # sys.stderr = open(params.out_dir + 'stderr.txt', 'w')

    pdn.go(params,name_mod=name_mod)
    ultp.go(params,name_mod=name_mod)
    plotme_sin.go(params,name_mod=name_mod)