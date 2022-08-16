import sys

import cy.pdn
import cy.ultp
import cy.plotme_sin

sys.path.insert(1,'../')

from param_defn import PD

def go(params: PD):
        
    # params.out_dir = '../'+params.out_dir[:params.out_dir.rfind('/')+1]
    #
    # print("write directory:", params.out_dir)
    #
    # sys.stdout = open(params.out_dir + 'stdout.txt', 'w')
    # sys.stderr = open(params.out_dir + 'stderr.txt', 'w')

    pdn.go(params)
    ultp.go(params)
    plotme_sin.go(params)