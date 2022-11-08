import numpy as np
from param_defn import PD
import convert_to_csvs as ccsv
from glob import glob

"""
Running still: 
0.602 3/32
0.602 2.8/25.4
0.602 1/8
1.029 1/8
1.59 1/8
1.029 3/32
1.029 5/32
1.59 1/4
1.59 3/16
0.26 0.7/25.4
Last 1.59 3/16
Missed 0.602 1/16
"""

ID_DP_dict_all = {'0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
              '0.602': ['1/4', '1/8', '1/16', '3/16', '7/16', '2.8/25.4', '3/32'],
              '1.029': ['1/4', '1/8', '3/16', '3/32', '5/16', '7/16', '15/32', '5/32',],
              '1.59': ['1/4', '1/8', '3/16', '5/16', '5/32', '7/16', '7/32', '9/32', '9/64', '15/64']}

ID_DP_dict = {'0.26': ['0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
              '0.602': ['1/4', '3/16', '7/16',],
              '1.029': ['1/4', '3/16', '5/16', '7/16', '15/32',],}

if __name__ == "__main__":
    base_dir = 'outputs/'
    for ID_str in ID_DP_dict.keys():
        print("ID:",ID_str)
        ID_dir = base_dir+ID_str.replace('.','pt')+'/'
        for Dp_str in ID_DP_dict[ID_str]:
            params = PD(ID_str,Dp_str,num_inserts=10)
            print("Dp:",Dp_str)
            # Dp_str = Dp_str.replace('.', 'pt').replace('/', '_')
            out_dir = ID_dir + Dp_str.replace('.', 'pt').replace('/', '_')+'/'
            dirs = glob(out_dir+'*/')
            for dir in dirs:
                dir_check = dir[len(out_dir):len(dir)-1]
                print(dir_check)
                if dir_check.__contains__('7-11-2022') or dir_check.__contains__('6-11-2022'):
                    ccsv.go(params.N_spheres(),params.final_step(),dir[:len(dir)-1],name_mod='_n10')
                    print("done")
                    break
            # print(dirs)