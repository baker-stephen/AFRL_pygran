import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

if __name__ =="__main__":
    xs = []
    press = []
    with open('axial_pressure_data/0pt26_ID_1pt7mm_DP_1_mdot_r_new.xy', 'r') as read:
        for i,line in enumerate(read):
            if i<4:
                # print("non-data line:",line)
                continue
            elif line == ")\n":
                break
            else:
                data = line.split()
                xs.append(float(data[0]))
                press.append(float(data[1]))
        read.close()

    # print(xs)
    # print(press)
    # Sort pressure based on x location
    press_sort = [p for x,p in sorted(zip(xs,press))]

    # Sort x location
    xs_sort = np.array(sorted(xs))
    # print("press sort")
    # print(press_sort)
    # print("x sort")
    # print(xs_sort)
    dpdx = (press_sort[len(press_sort)-1]-press_sort[0])/(xs_sort[-1]-xs_sort[0])

    print("dpdx: ", dpdx)

    #linear fit:
    model = LinearRegression()

    xs_model = xs_sort.reshape((-1,1))
    press_model = np.array(press_sort)
    model.fit(xs_model,press_model)
    print("linear regression model slope: ",model.coef_[0])
    print("fit quality: ",model.score(xs_model,press_model))

    #convert pressure to psi:
    # press_sort = [p/6894.76 for p in press_sort]

    #convert position to in, set x0 = 0:
    xs_sort_in = np.array([x*.0254 for x in xs_sort])
    fig, ax = plt.subplots()
    ax.plot(xs_sort_in,press_sort,label="data")
    ax.plot(xs_sort_in,model.intercept_+xs_sort*model.coef_[0],label="linear fit")
    ax.set_xlabel("radius (in)")
    ax.set_ylabel("pressure (Pa)")
    ax.legend()
    plt.show()