import numpy as np

def sin_cos_list(x,params):
    ret = params[0]
    for i in range(1,len(params)):
        if i % 2 != 0:
            ret += params[i]*np.cos(2 * np.pi * (i//2 + 1) * x)
        else:
            ret += params[i]*np.sin(2 * np.pi * (i//2) * x)
    return ret

if __name__ == "__main__":

    # ID_DP_dict = {
    #     '0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4',
    #              '1.7/25.4',],
    #     '0.602': ['1/4', '1/8', '2.8/25.4', '3/16', '3/32', '7/16',], #Post processing missing n10 '1/16'
    #     '1.029': ['1/4', '1/8', '3/16', '3/32', '5/16', '5/32', '7/16', '15/32', ],
    #     '1.59': ['1/4', '1/8', '3/16', '5/16', '7/16', '7/32', '9/32', '9/64', '15/64',]} #n10 pygran missing '5/32', '7/32', '9/64', '15/64'

    ID_DP_dict = {
        '1.59': ['9/32']}

    for ID in ID_DP_dict:
        # wrt_dir = '../all_data/'
        wrt_dir = ID.replace('.', 'pt') + '/'
        for DP in ID_DP_dict[ID]:
            this_wrt_dir = wrt_dir+DP.replace('.', 'pt').replace('/', '_') + '/'
            print("read/write directory:", this_wrt_dir)

            fourier_fit_str = ""
            save_next = False
            with open(this_wrt_dir+'outputs_n10.txt', 'r') as op_file:
                for line_str in op_file:
                    line_str = str(line_str).strip()
                    if save_next:
                        fourier_fit_str = line_str
                        save_next = False
                    elif line_str.lower().__contains__("fourier fit"):
                        save_next = True

                op_file.close()

            if fourier_fit_str == "":
                print("no fourier fit found for DP = " + DP + ", ID = " + ID + ", skipping.")
            else:
                term_strs = fourier_fit_str.split("+")
                params = [float(term_strs[0])]
                for i,term in enumerate(term_strs):
                    if i==0:
                        continue
                    else:
                        splt_term = term.split("*")
                        params.append(float(splt_term[0]))

                print("got fourier fit")

                n = 10000
                r_fits = np.linspace(0, float(ID)/2, n)
                p_fits = []
                for i, r in enumerate(r_fits):
                    p_fits.append(sin_cos_list(r, params))

                pfmean = np.mean(p_fits)
                print("mean from generated pts:", pfmean)
                mean_weight = np.average(p_fits,weights=np.linspace(0,1,n))
                print("weighted mean:", mean_weight)

                wght_vari = 0
                weight_sum = 0
                for i, p in enumerate(p_fits):
                    w = i/(n+1)
                    weight_sum += w
                    wght_vari += w*(p-mean_weight)**2
                wght_vari /= ((n-1)/n)*weight_sum
                print("weighted variance:",wght_vari)
                wght_std = np.sqrt(wght_vari)
                print("weighted std dev:",wght_std)


                # out.write("mean from generated pts: " + str(np.mean(p_fits)) + '\n')
                print("Normal std dev from generated pts:", np.std(p_fits))
                # out.write("std dev from generated pts: " + str(np.std(p_fits)) + '\n')

                with open(this_wrt_dir + 'outputs_n10.txt', 'a') as op_file:
                    op_file.write("\n\nweighted std dev: "+str(wght_std)+'\n')
                    op_file.close()
