import pdn
import ultp
import plotme_sin

import sys

if __name__ == "__main__":

    # ID_DP_dict = {'0.26': ['0.7/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
    #               '1.029': ['1/4', '1/8', '3/16', '7/16'],
    #               '1.59': ['1/4', '1/8', '3/16', '5/16', '5/32', '7/16', '7/32', '9/32', '9/64', '15/64']}


    ID_DP_dict = {'0.602': ['1/4', '1/8', '3/16',],}



    for ID in ID_DP_dict:
        file_name = 'in_files/hi_res-'
        file_name += ID.replace('.', 'pt') + '-'
        wrt_dir = '../all_data/'
        wrt_dir += ID.replace('.', 'pt') + '/'
        for DP in ID_DP_dict[ID]:
            this_file = file_name + DP.replace('.', 'pt').replace('/', '_') + '.in'
            this_wrt_dir = wrt_dir+DP.replace('.', 'pt').replace('/', '_') + '/'
            print("write directory:", this_wrt_dir)
            print("input file:",this_file)
            sys.stdout = open(this_wrt_dir + 'stdout.txt', 'w')
            sys.stderr = open(this_wrt_dir + 'stderr.txt', 'w')
            pdn.go(this_file)
            ultp.go(this_file)
            plotme_sin.go(this_file)