
def output_vti(postions:list,tstep:int):
    with open('csvs_new/particles'+str(tstep)+'.csv','w') as write:
        write.write("x,y,z\n")
        for part in range(len(postions)):
            if len(postions[part]) < 3:
                continue
            write.write(str(postions[part][0]) + ',' + str(postions[part][1]) + ',' + str(postions[part][2]) + '\n')

        write.close()


if __name__ == "__main__":
    atom_count = 100
    for step in range(1000,2000001,1000):
        positions = []
        if step%10000==0:
            print("step: ",step)
        with open('out-SpringDashpot-12:53:55-14.5.2022/traj/particles'+str(step)+'.vtk', 'r') as read:

            this_atom_count = 0
            for i in range(9):
                if i == 3:
                    this_atom_count = int(read.readline())
                    # print("this_atom_count:",this_atom_count)
                else:
                    read.readline()
            for i in range(atom_count):
                if i<this_atom_count:
                    positions.append([])
                    data = read.readline().split()
                    positions[i].append(float(data[2]))
                    positions[i].append(float(data[3]))
                    positions[i].append(float(data[4]))
                else:
                    positions.append([-1, -1, 0])
            output_vti(positions,step)
            read.close()

    print("done")