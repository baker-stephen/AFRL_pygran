import os

#TODO: Comment this file

def output_csv(postions: list, tstep: int):
    with open('csvs/particles' + str(tstep) + '.csv', 'w') as write:
        write.write("x,y,z\n")
        for part in range(len(postions)):
            if len(postions[part]) < 3:
                continue
            write.write(str(postions[part][0]) + ',' + str(postions[part][1]) + ',' + str(postions[part][2]) + '\n')

        write.close()


def go(atom_count: int, final_step: int, source_dir: str):
    try:
        os.makedirs(source_dir + '/csvs/', exist_ok=True)
        print("Directory '/csvs/' created successfully\n")
    except OSError as error:
        print("Directory '/csvs/' can not be created. Error: %s\n" % str(error))
    box_bounds = []
    for step in range(50000, final_step, 50000):
        positions = []
        if step % 100000 == 0:
            print("step: ", step)
        with open(source_dir + '/traj/particles' + str(step) + '.vtk', 'r') as read:
            this_atom_count = 0
            for i in range(9):
                if i == 3:
                    this_atom_count = int(read.readline())
                elif 5 <= i <= 7 and len(box_bounds) < 3:
                    box_bounds.append([float(b) for b in read.readline().split()])
                else:
                    read.readline()
            for i in range(atom_count):
                if i < this_atom_count:
                    positions.append([])
                    data = read.readline().split()
                    positions[i].append(float(data[2]))
                    positions[i].append(float(data[3]))
                    positions[i].append(float(data[4]))
                else:
                    positions.append([box_bounds[0][0], box_bounds[1][0], box_bounds[2][0]])
            output_csv(positions, step)
            read.close()

    print("Done converting to csvs")
