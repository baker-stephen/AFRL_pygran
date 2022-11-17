import numpy as np
from param_defn import PD
import matplotlib.pyplot as plt

ID_DP_dict_all = {
    '0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
    '0.602': ['1/4', '1/8', '1/16', '3/16', '7/16', '2.8/25.4', '3/32'],
    '1.029': ['1/4', '1/8', '3/16', '3/32', '5/16', '7/16', '15/32', '5/32', ],
    '1.59': ['1/4', '1/8', '3/16', '5/16', '5/32', '7/16', '7/32', '9/32', '9/64', '15/64']}

# ID_DP_dict = {'0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
#               '0.602': ['1/4', '1/8', '1/16', '3/16', '7/16', '2.8/25.4','3/32'], #Not ready yet: '3/32'
#               '1.029': ['1/4', '1/8', '3/16', '3/32', '5/16', '7/16', '15/32'], #Not ready yet: '5/32',
#               '1.59': ['1/4', '1/8', '3/16', '5/16', '5/32', '7/16', '7/32', '9/32', '9/64', '15/64']}

ID_DP_dict = {
    '0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
    '0.602': ['1/4', '1/8', '3/16', '7/16', '2.8/25.4', ], #1/16 still running, '3/32', needs post
    '1.029': ['1/4', '1/8', '3/16', '3/32', '5/16', '7/16', '15/32', ], #'5/32', needs post
    '1.59': ['1/4', '1/8', '3/16', '5/16', '7/16', '9/32']} #'5/32', '7/32', '9/64', '15/64' still running

# ID_DP_dict = {'0.26': ['0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
#               '0.602': ['1/4', '3/16', '7/16',],
#               '1.029': ['1/4', '3/16', '5/16', '7/16', '15/32',],}

if __name__ == "__main__":
    csv = open('guo_layers.csv', 'w')
    csv.write("D,Dp,N,vol avg por,n,l,k1,k2\n")

    Ns = []
    k1s = []

    base_dir = 'outputs/'
    for D_str in ID_DP_dict.keys():
        ID_dir = base_dir + D_str.replace('.', 'pt') + '/'
        for Dp_str in ID_DP_dict[D_str]:

            print("D: "+D_str+", Dp: "+Dp_str)

            params = PD(D_str,Dp_str,num_inserts=10)

            out_dir = ID_dir + Dp_str.replace('.', 'pt').replace('/', '_')+'/'

            # params.output_dir()

            # wrt_dir = params.out_dir[:params.out_dir.rfind('/') + 1]

            pipe_in_rad = params.ID / 2
            D = params.ID
            Dp = params.DP

            res = params.res()
            x_res = res[0]
            x0 = -pipe_in_rad
            xM = pipe_in_rad
            dx = (xM - x0) / (x_res + 1)

            y_res = res[1]
            y0 = -pipe_in_rad
            yM = pipe_in_rad
            dy = (yM - y0) / (y_res + 1)

            z_res = res[2]
            z0 = params.z0_cutoff()
            zM = params.zM_cutoff()
            dz = (zM - z0) / (z_res + 1)

            positions = []
            final_csv = params.final_step()-params.final_step()%50000
            with open(out_dir+'csvs_n10/particles'+str(final_csv)+'.csv', 'r') as r:
            # with open(out_dir + 'particles.csv', 'r') as r:
                r.readline()
                for line in r:
                    xyz = [float(p) for p in line.split(",")]
                    if xyz[2]<z0 or xyz[2]>zM:
                        continue
                    positions.append(xyz)
                r.close()

            positions = sorted(positions, key=lambda x: x[2])

            layers = []
            layer_pos = []
            z_diffs = [0.0]

            for i,xyz in enumerate(positions):
                z = xyz[2]
                for j,xyz_2 in enumerate(positions[i+1:]):
                    z2 = xyz_2[2]
                    diff = z2-z
                    if diff >= Dp:
                        break
                    elif len(z_diffs) <= j:
                        z_diffs.append(diff)
                    else:
                        z_diffs[j] = (z_diffs[j]+diff)/2

            z_diffs = [zd / Dp for zd in z_diffs]

            # L2s = []
            n_fromz = 0
            for i in range(len(z_diffs)-1,0,-1):
                z_model = np.array(z_diffs[i:])
                zmavg = np.mean(z_model)
                L2 = 0
                for z in z_model:
                    L2 += (z-zmavg)**2
                if L2/len(z_model) > 0.001:
                    n_fromz = i
                    break
                # L2s.append(L2/len(z_model))

                # print("linear regression model slope: ", model.coef_[0])
                # print("fit quality: ", model.score(xs_model, press_model))
            print("n_fromz:", n_fromz)
            if n_fromz==0:
                continue
            else:
                Ns.append(D/Dp)
            k1 = n_fromz*(Dp**2)/(D**2)
            k1s.append(k1)
            print("k1:",k1)

            vol_avg_por = -1
            with open(out_dir+'outputs.txt', 'r') as out_txt:
                key_str = "volume averaged porosity: "
                for line in out_txt:
                    if line.__contains__(key_str):
                        vol_avg_por_str = line[len(key_str)+1:]

                vol_avg_por = float(vol_avg_por_str)

                out_txt.close()

            # z_k1 = 1.5*(1-vol_avg_por)/z_k2
            # print("z_k1:",z_k1)


            # fig, ax = plt.subplots()
            # ax.plot([zd/Dp for zd in z_diffs])
            # ax.set_xlabel("spheres away")
            # ax.set_ylabel("z diff / Dp")
            # ax.set_title("D: "+D_str+", Dp: "+Dp_str)
            # plt.show()

            # fig, ax = plt.subplots()
            # ax.plot(L2s)
            # ax.set_xlabel("spheres away")
            # ax.set_ylabel("L2s")
            # ax.set_title("D: " + D_str + ", Dp: " + Dp_str)
            # plt.show()

            # continue

            #
            # for xyz in positions:
            #     z = xyz[2]
            #     if len(layers) == 0:
            #         layers.append([xyz])
            #         layer_pos.append([z,z,z])
            #     else:
            #         did_add = False
            #         for li,lp in enumerate(layer_pos):
            #             if abs(z-lp[0]) < np.sqrt(3)*Dp/2:
            #                 did_add=True
            #                 layers[li].append(xyz)
            #                 if z<lp[0]:
            #                     lp[0] = z
            #                 elif z>lp[1]:
            #                     lp[1] = z
            #                 avg = 0
            #                 for pos in layers[li]:
            #                     avg += pos[2]
            #                 lp[2] = avg/len(layers[li])
            #                 break
            #
            #         if not did_add:
            #             layers.append([xyz])
            #             layer_pos.append([z, z, z])
            # ns = [len(layer) for layer in layers]
            # print("ns:",ns)
            # print("zs:",layer_pos)
            # print("ns sort:",sorted(ns))
            # n_avg = np.mean(ns)
            # # with open(out_dir + 'pos_layers.csv', 'w') as w:
            # #     w.write("x,y,z,layer\n")
            # #     for i,layer in enumerate(layers):
            # #         for pos in layer:
            # #             w.write(str(pos[0])+","+str(pos[1])+","+str(pos[2])+","+str(i)+"\n")
            # #     w.close()
            # # print("n avg:",n_avg)
            #
            # #Alternate n calc method:
            # # n_alts_counter = 0
            # # min_zi = 0
            # # for pi,xyz in enumerate(positions):
            # #     z = xyz[2]
            # #     zi = min_zi
            # #     min_zi_set = False
            # #     while zi<len(positions):
            # #         z_other = positions[zi][2]
            # #         if z-z_other < Dp/4:
            # #             zi += 1
            # #             continue
            # #         elif z_other-z > Dp/4:
            # #             break
            # #         elif not min_zi_set:
            # #             min_zi = zi
            # #             min_zi_set = True
            # #
            # #         n_alts_counter += 1
            # #         zi+=1
            # n_std = np.std(ns)
            # print("std:",n_std)
            # n_alts = []
            # for n in ns:
            #     if abs(n - n_avg) <= n_std:
            #         n_alts.append(n)
            # n_alt = np.mean(n_alts)
            # print("n_alt:",n_alt)
            #
            # ls = []
            # for i in range(1,len(layer_pos)-1):
            #     prev = layer_pos[i][2] - layer_pos[i-1][2]
            #     next = layer_pos[i+1][2] - layer_pos[i][2]
            #     ls.append((prev+next)/2)
            #
            # # print("ls:",ls)
            # l_avg = np.mean(ls)
            # print("l avg:",l_avg)

            # k1 = n_avg*Dp**2 / D**2
            # k2 = Dp/l_avg

            # porosity = 1 - (2/3)*k1*k2

            N = D/Dp


            # l_alt = 2*n_alt*Dp**3/(3*D**2*(1-vol_avg_por))

            # rel_err = 2*abs(vol_avg_por - porosity)/(vol_avg_por + porosity)

            # k1 = n_alt * (Dp**2)/(D**2)
            # print("k1:",k1)

            k2 = 1.5 * (1 - vol_avg_por) / k1
            l = Dp/k2
            csv.write(D_str+','+str(Dp)+','+str(N)+','+str(vol_avg_por)+','+str(n_fromz)+','+str(l)+','+str(k1)+','+str(k2)+'\n')
            print()

    csv.close()
    fig, ax = plt.subplots()
    ax.scatter(Ns,k1s)
    ax.set_xlabel("N")
    ax.set_ylabel("k1")
    # ax.set_title("D: " + D_str + ", Dp: " + Dp_str)
    plt.show()

