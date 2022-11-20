from fluent import poros_Guo, poros_Sato, poros_Cheng, poros_Foumeny
import numpy as np
import matplotlib.pyplot as plt
from param_defn import PD

Guo_Data = \
    "1.183 0.5120 0.5039 0.4961 0.5040 ± 0.0080\
    1.263 0.5632 0.5650 0.5574 0.5619 ± 0.0040\
    1.391 0.6385 0.6310 0.6195 0.6297 ± 0.0096\
    1.667 0.6735 0.6720 0.6885 0.6780 ± 0.0091\
    1.718 0.6740 0.6750 0.6891 0.6794 ± 0.0084\
    1.775 0.6422 0.6363 0.6490 0.6425 ± 0.0064\
    2.525 0.4971 0.4927 0.4861 0.4920 ± 0.0055\
    2.727 0.4750 0.4811 0.4814 0.4792 ± 0.0036\
    2.880 0.4920 0.5007 0.4989 0.4972 ± 0.0046\
    3.117 0.4198 0.4235 0.4280 0.4238 ± 0.0041\
    3.367 0.4640 0.4691 0.4721 0.4684 ± 0.0041\
    3.408 0.4601 0.4592 0.4555 0.4583 ± 0.0024\
    3.556 0.4467 0.4415 0.4405 0.4429 ± 0.0033\
    3.752 0.4613 0.4650 0.4669 0.4644 ± 0.0028\
    4.090 0.4285 0.4290 0.4307 0.4294 ± 0.0012\
    4.208 0.4309 0.4321 0.4351 0.4327 ± 0.0022\
    4.383 0.4256 0.4271 0.4214 0.4247 ± 0.0030\
    4.633 0.4267 0.4251 0.4211 0.4243 ± 0.0029\
    4.733 0.4153 0.4140 0.4122 0.4138 ± 0.0016\
    5.000 0.4047 0.4056 0.4014 0.4039 ± 0.0022\
    5.050 0.4079 0.4073 0.4067 0.4073 ± 0.0006\
    5.113 0.4305 0.4319 0.4325 0.4316 ± 0.0010\
    5.917 0.3935 0.3907 0.3924 0.3922 ± 0.0014\
    6.313 0.4073 0.4041 0.4051 0.4055 ± 0.0016\
    6.514 0.3896 0.3855 0.3862 0.3871 ± 0.0022\
    7.100 0.4009 0.4052 0.4038 0.4033 ± 0.0022"

# Pygran_N_strs = \
#     '9.434285714\
#     9.434285714\
#     9.434285714\
#     8.255\
#     8.255\
#     8.255\
#     7.769411765\
#     7.769411765\
#     7.769411765\
#     7.337777778\
#     7.337777778\
#     7.337777778\
#     6.604\
#     6.604\
#     6.604\
#     6.003636364\
#     6.003636364\
#     6.003636364\
#     5.503333333\
#     5.503333333\
#     5.503333333\
#     4.402666667\
#     4.402666667\
#     4.402666667\
#     3.884705882\
#     3.884705882\
#     3.884705882\
#     4.816\
#     4.816\
#     4.816\
#     3.210666667\
#     3.210666667\
#     3.210666667\
#     2.408\
#     2.408\
#     2.408\
#     10.976\
#     10.976\
#     10.976\
#     8.232\
#     8.232\
#     8.232\
#     5.488\
#     5.488\
#     5.488\
#     4.116\
#     4.116\
#     4.116\
#     3.2928\
#     3.2928\
#     3.2928\
#     2.352\
#     2.352\
#     2.352\
#     12.72\
#     12.72\
#     12.72\
#     11.30666667\
#     11.30666667\
#     11.30666667\
#     10.176\
#     10.176\
#     10.176\
#     8.48\
#     8.48\
#     8.48\
#     7.268571429\
#     7.268571429\
#     7.268571429\
#     6.784\
#     6.784\
#     6.784\
#     6.36\
#     6.36\
#     6.36\
#     5.653333333\
#     5.653333333\
#     5.653333333\
#     5.088\
#     5.088\
#     5.088\
#     3.634285714\
#     3.634285714\
#     3.634285714\
#     9.632\
#     9.632\
#     9.632\
#     1.376\
#     1.376\
#     1.376'
#
# Pygran_totals_str = \
#     '0.417075328\
#     0.417124131\
#     0.423372709\
#     0.433506542\
#     0.42204091\
#     0.418303156\
#     0.447253535\
#     0.444293468\
#     0.46172313\
#     0.413872601\
#     0.453256365\
#     0.469135768\
#     0.387944454\
#     0.459424899\
#     0.472640754\
#     0.48452612\
#     0.466998421\
#     0.49125576\
#     0.368853806\
#     0.373714098\
#     0.376171953\
#     0.383219076\
#     0.387200962\
#     0.393487057\
#     0.392108968\
#     0.40074237\
#     0.411979343\
#     0.447279093\
#     0.389611210039648\
#     0.614629929191254'
#
# Pygran_vol_avg_str = \
#     '0.401253773\
#     0.401253773\
#     0.401253773\
#     0.40286297\
#     0.40286297\
#     0.40286297\
#     0.406926964\
#     0.406926964\
#     0.406926964\
#     0.413968424\
#     0.413968424\
#     0.413968424\
#     0.416128576\
#     0.416128576\
#     0.416128576\
#     0.399367731\
#     0.399367731\
#     0.399367731\
#     0.429277902\
#     0.429277902\
#     0.429277902\
#     0.43062107\
#     0.43062107\
#     0.43062107\
#     0.44720406\
#     0.44720406\
#     0.44720406\
#     0.397493138\
#     0.397493138\
#     0.397493138\
#     0.441156077\
#     0.441156077\
#     0.441156077\
#     0.466470357\
#     0.466470357\
#     0.466470357\
#     0.368251872\
#     0.368251872\
#     0.368251872\
#     0.375201072\
#     0.375201072\
#     0.375201072\
#     0.392447128\
#     0.392447128\
#     0.392447128\
#     0.4706465\
#     0.4706465\
#     0.4706465\
#     0.45224866\
#     0.45224866\
#     0.45224866\
#     0.470458047\
#     0.470458047\
#     0.470458047\
#     0.349582882\
#     0.349582882\
#     0.349582882\
#     0.354139427\
#     0.354139427\
#     0.354139427\
#     0.357693398\
#     0.357693398\
#     0.357693398\
#     0.366362826\
#     0.366362826\
#     0.366362826\
#     0.371240033\
#     0.371240033\
#     0.371240033\
#     0.378325404\
#     0.378325404\
#     0.378325404\
#     0.376659882\
#     0.376659882\
#     0.376659882\
#     0.386272807\
#     0.386272807\
#     0.386272807\
#     0.398651618\
#     0.398651618\
#     0.398651618\
#     0.435471127\
#     0.435471127\
#     0.435471127\
#     0.372539998\
#     0.372539998\
#     0.372539998\
#     0.61016592\
#     0.61016592\
#     0.61016592'

ID_DP_dict = {
    '0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
    '0.602': ['1/4', '1/8', '1/16', '3/16', '7/16', '2.8/25.4', '3/32'],
    '1.029': ['1/4', '1/8', '3/16', '3/32', '5/16', '7/16', '15/32', '5/32', ],
    '1.59': ['1/4', '1/8', '3/16', '5/16', '5/32', '7/16', '7/32', '9/32', '9/64', '15/64']}

def get_vol_avg(file)->float:
    key = "volume averaged porosity: "
    por = -1
    for line in file:
        if line.startswith(key):
            por_str = line[len(key)+1:].strip()
            try:
                por = float(por_str)
            except Exception as err:
                continue

    return por

if __name__ == "__main__":

    fig, ax = plt.subplots()

    res = 500
    N_max = 12.8

    N_min_sato = 2.5
    N_satos = np.linspace(N_min_sato,N_max,res)
    ax.plot(N_satos,[poros_Sato(N) for N in N_satos], label="Sato")

    N_min_Foumeny = 1.866
    N_Foumenys = np.linspace(N_min_Foumeny,N_max,res)
    ax.plot(N_Foumenys, [poros_Foumeny(N) for N in N_Foumenys], label="Foumeny")

    N_min_Cheng = 1.01
    N_Chengs = np.linspace(N_min_Cheng, N_max,res)
    ax.plot(N_Chengs, [poros_Cheng(N) for N in N_Chengs], label="Cheng")

    N_min_Guo = 1.0
    N_Guos = np.linspace(N_min_Guo,N_max,res)
    guos = [poros_Guo(N) for N in N_Guos]
    guos_bools = [guo!=-1 for guo in guos]
    N_Guos_plottable = N_Guos[guos_bools]
    guo_plottable = np.array(guos)[guos_bools]
    ax.scatter(N_Guos_plottable, guo_plottable, label="Guo", s=5, c='y', marker="_")

    guo_lines = Guo_Data.split("    ")
    Guo_data_Ns = []
    Guo_data_avgs = []
    Guo_data_stdevs = []
    for line in guo_lines:
        vals = line.split(" ")
        Guo_data_Ns.append(float(vals[0]))
        Guo_data_avgs.append(float(vals[4]))
        Guo_data_stdevs.append(float(vals[6]))

    ax.errorbar(Guo_data_Ns,Guo_data_avgs, yerr=Guo_data_stdevs, label="Guo experimental data", fmt="o", c="black")

    pg_Ns_single = []
    pg_Ns_10 = []
    pg_single = []
    pg_10 = []
    for Dstr in ID_DP_dict.keys():
        for Dpstr in ID_DP_dict[Dstr]:
            p = PD(Dstr,Dpstr)
            dirr = '../'+p.out_dir

            try:
                outsingle = open(dirr+'outputs.txt')
                va_por = get_vol_avg(outsingle)
                if va_por != -1:
                    pg_Ns_single.append(p.ID / p.DP)
                    pg_single.append(va_por)
            except Exception as err:
                print("No single drop file for "+Dstr+", "+Dpstr+". Exc:",err)

            try:
                out10 = open(dirr+'outputs_n10.txt')
                va_por = get_vol_avg(out10)
                if va_por != -1:
                    pg_Ns_10.append(p.ID / p.DP)
                    pg_10.append(va_por)
            except Exception as err:
                print("No 10 drop file for "+Dstr+", "+Dpstr+". Exc:",err)

    ax.scatter(pg_Ns_single, pg_single, label="Pygran single drop", c="red")
    ax.scatter(pg_Ns_10, pg_10, label="Pygran 10 drop", c="purple")


    ax.legend()
    ax.set_title("Porosity Comparisons")
    ax.set_xlabel("Inner Diameter / Sphere Diameter")
    ax.set_ylabel("Porosity")
    plt.savefig('porosity_comparisons.png')
    plt.show()



