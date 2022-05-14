
def output_vti(postions:list,tstep:int):
    with open('csvs_new/particles'+str(tstep)+'.csv','w') as write:
        write.write("x,y,z\n")
        for part in range(len(postions)):
            if len(postions[part]) < 3:
                continue
            write.write(str(postions[part][0]) + ',' + str(postions[part][1]) + ',' + str(postions[part][2]) + '\n')

        write.close()


if __name__ == "__main__":
    for step in range(1000,100001,1000):
        positions = []
        if step%10000==0:
            print("step: ",step)
        with open('out-SpringDashpot-20:20:20-13.5.2022/traj/particles'+str(step)+'.vtk', 'r') as read:
            atom_count = 6
            for i in range(9):
                if i == 3:
                    atom_count = int(read.readline())
                else:
                    read.readline()
            for i in range(atom_count):
                positions.append([])
                data = read.readline().split()
                positions[i].append(float(data[2]))
                positions[i].append(float(data[3]))
                positions[i].append(float(data[4]))
            output_vti(positions,step)
            read.close()

    print("done")