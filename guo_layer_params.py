import numpy as np
from param_defn import PD

ID_DP_dict = {'0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
              '0.602': ['1/4', '1/8', '1/16', '3/16', '7/16'], #Not ready yet: '3/32', '2.8/25.4',
              '1.029': ['1/4', '1/8', '3/16', '3/32', '5/16', '7/16', '15/32'], #Not ready yet: '5/32',
              '1.59': ['1/4', '1/8', '3/16', '5/16', '5/32', '7/16', '7/32', '9/32', '9/64', '15/64']}

if __name__ == "__main__":
    csv = open('guo_layers.csv', 'w')
    csv.write("D,Dp,N,por pg,por nl,rel err,n,l\n")

    for D_str in ID_DP_dict.keys():
        for Dp_str in ID_DP_dict[D_str]:

            print("D: "+D_str+", Dp: "+Dp_str)

            params = PD(D_str, Dp_str)

            params.output_dir()

            wrt_dir = params.out_dir[:params.out_dir.rfind('/') + 1]

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
            with open(wrt_dir+'particles.csv', 'r') as r:
                r.readline()
                for line in r:
                    xyz = [float(p) for p in line.split(",")]
                    positions.append(xyz)
                r.close()

            positions = sorted(positions, key=lambda x: x[2])

            layers = []
            layer_pos = []
            for xyz in positions:
                z = xyz[2]
                if len(layers) == 0:
                    layers.append([xyz])
                    layer_pos.append([z,z,z])
                else:
                    did_add = False
                    for li,lp in enumerate(layer_pos):
                        if abs(z-lp[2]) < Dp/4:
                            did_add=True
                            layers[li].append(xyz)
                            if z<lp[0]:
                                lp[0] = z
                            elif z>lp[1]:
                                lp[1] = z
                            avg = 0
                            for pos in layers[li]:
                                avg += pos[2]
                            lp[2] = avg/len(layers[li])
                            break

                    if not did_add:
                        layers.append([xyz])
                        layer_pos.append([z, z, z])

            ns = [len(layer) for layer in layers]
            print("ns:",ns)
            n_avg = np.mean(ns)
            print("n avg:",n_avg)

            ls = []
            for i in range(1,len(layer_pos)-1):
                prev = layer_pos[i][2] - layer_pos[i-1][2]
                next = layer_pos[i+1][2] - layer_pos[i][2]
                ls.append((prev+next)/2)

            print("ls:",ls)
            l_avg = np.mean(ls)
            print("l avg:",l_avg)

            k1 = n_avg*Dp**2 / D**2
            k2 = Dp/l_avg

            porosity = 1 - (2/3)*k1*k2

            N = D/Dp
            vol_avg_por = -1
            with open(wrt_dir+'outputs.txt', 'r') as out_txt:
                key_str = "volume averaged porosity: "
                for line in out_txt:
                    if line.__contains__(key_str):
                        vol_avg_por_str = line[len(key_str)+1:]

                vol_avg_por = float(vol_avg_por_str)

                out_txt.close()

            rel_err = 2*abs(vol_avg_por - porosity)/(vol_avg_por + porosity)

            csv.write(D_str+','+Dp_str+','+str(N)+','+str(vol_avg_por)+','+str(porosity)+','+str(rel_err)+','+str(n_avg)+','+str(l_avg)+'\n')
            print()

    csv.close()
