import sys

import convert_to_csvs
import cy.run_post as rp

from param_defn import PD

ID_DP_dict_all = {
    '0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
    '0.602': ['1/4', '1/8', '1/16', '3/16', '7/16', '2.8/25.4', '3/32'],
    '1.029': ['1/4', '1/8', '3/16', '3/32', '5/16', '7/16', '15/32', '5/32', ],
    '1.59': ['1/4', '1/8', '3/16', '5/16', '5/32', '7/16', '7/32', '9/32', '9/64', '15/64']}

ID_DP_dict = {
    '0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
    '0.602': ['3/32', '1/4', '1/8', '3/16', '7/16', '2.8/25.4', ], #1/16 still running
    '1.029': ['5/32', '1/4', '1/8', '3/16', '3/32', '5/16', '7/16', '15/32',],
    '1.59': ['1/4', '1/8', '3/16', '5/16', '7/16', '9/32']} #'5/32', '7/32', '9/64', '15/64' still running

if __name__ == "__main__":
    num_insertions = 10
    # Define the parameters of the simulation
    sys.stdout = open('outputs/stdout_post_mutli.txt', 'w')
    sys.stderr = open('outputs/sterr_post_mutli.txt', 'w')
    for Dstr in ID_DP_dict.keys():
        for Dpstr in ID_DP_dict[Dstr]:
            params = PD(Dstr, Dpstr, num_inserts=num_insertions)

            rp.go(params,name_mod="n10")

    sys.stdout.close()
    sys.stderr.close()


