from pygran import simulation
from pygran.params import steel, glass
import numpy as np
import convert_to_csvs

import sys

from param_defn import PD

def go(sim_params: PD):
    parts_per_insert = sim_params.N_spheres() // sim_params.num_inserts
    print('Running ' + sim_params.ID_str + '" ID with atom count ' + str(
        sim_params.N_spheres()) + ', and ' + sim_params.DP_str + '" DP.\n')
    box = sim_params.bounds()
    min = box[0]
    max = box[5]
    all_params = (min,max)
    center = (all_params[0]+all_params[1])/2
    c = str(center)
    box_l = [min,max,min,max,min,max,]
    # Create a dictionary of physical parameters
    params = {

        # Define the system
        'boundary': ('f', 'f', 'f'),  # fixed BCs

        'box': box_l,  # simulation box size in inches

        # Define component(s)
        'species': (
            {'material': glass, 'style': 'sphere', 'radius': sim_params.DP / 2},),

        # Set skin distance to be 1/4 particle diameter
        'nns_skin': sim_params.DP / 4,

        # Timestep
        'dt': 2.5e-7,

        # Apply gravitional force in the negative direction along the z-axis. In inches per second squared
        'gravity': (385.827, 0, 0, -1),

        # Setup I/O
        'traj': {'pfile': 'particles*.vtk', 'mfile': 'pipe*.vtk', 'freq': 50000},

        'output': sim_params.out_dir,

        # Stage runs
        'stages': {'insertion': 2.3e6 / sim_params.num_inserts},

        # Define mesh
        'mesh': {
            'pipe': {'file': sim_params.mesh(), 'mtype': 'mesh/surface/stress', 'material': steel,
                     'args': {'curvature_tolerant': 'yes'}
                     },
        },

        # 'restart' : [0,'outputs/0pt602/2pt8_25pt4/sim_out_12:12:16_6-11-2022/restart/restart.binary.48955000']
    }

    print("create sim")
    # Create an instance of the DEM class
    sim = simulation.DEM(**params)

    print("read restart")

    sim.command('read_restart outputs/0pt602/2pt8_25pt4/sim_out_12:12:16_6-11-2022/restart/restart.binary.48955000')

    print("read done.")
    # #Add a dissipative force
    # air_resistance = sim.addViscous(species=1, gamma=0.1)
    #
    # for i in range(sim_params.num_inserts):
    #     # Insert a group of particles
    #     insert = sim.insert(species=1, value=((i+1)*parts_per_insert - sim.get_natoms()), region=sim_params.insert(),
    #                         args={'orientation': 'random'})
    #
    #     # Run insertion stage, let the particles settle into the cylinder
    #     sim.run(params['stages']['insertion'] * 2, params['dt'])
    #     sim.remove(insert)
    #
    #     # Setup shaking:
    #     freq = 10 * 2 * np.pi
    #     nTaps = 30
    #     period = 1 / freq
    #     nSteps = period / params['dt']
    #     ampz = sim_params.ampz()
    #     ampxy = sim_params.ampxy()
    #
    #     for j in range(nTaps // 2):
    #         # vibrate x
    #         mm = sim.moveMesh('pipe', viblin=(
    #             'axis 1 0 0', 'order 1', 'amplitude {}'.format(ampxy), 'phase 0', 'period {}'.format(period)))
    #         sim.run(nSteps, params['dt'])
    #         sim.remove(mm)
    #         # vibrate y
    #         mm = sim.moveMesh('pipe', viblin=(
    #             'axis 0 1 0', 'order 1', 'amplitude {}'.format(ampxy), 'phase 0', 'period {}'.format(period)))
    #         sim.run(nSteps, params['dt'])
    #         sim.remove(mm)
    #
    #     # Allow for some settling
    #     sim.run(params['stages']['insertion'] / 2, params['dt'])
    #
    #     for k in range(nTaps):
    #         # vibrate z
    #         mm = sim.moveMesh('pipe', viblin=(
    #             'axis 0 0 1', 'order 1', 'amplitude {}'.format(ampz), 'phase 0', 'period {}'.format(period)))
    #         sim.run(nSteps, params['dt'])
    #         sim.remove(mm)
    #
    #     # Let simulation settle before next insertion
    #
    #     sim.run(params['stages']['insertion'] * 2, params['dt'])
    #
    # sim.run(params['stages']['insertion'] * 2, params['dt'])

    rotate = sim.moveMesh('pipe', rotate=(
        'origin '+c+" "+c+" "+c,'axis 1 0 0', 'period 1'
    ))
    sim.run((1.0/4)/params['dt'],params['dt'])
    sim.remove(rotate)
    sim.run(params['stages']['insertion'], params['dt'])




if __name__ == "__main__":
    # Retrieve command-line arguments. First element is always file name, we can skip that.
    args = sys.argv[1:]
    print('args:', args)

    if len(args) != 2 and len(args) != 3:
        raise Exception("Please specify both ID and DP")

    # Define the parameters of the simulation
    if len(args) == 3:
        num_insertions = int(str(args[2]).strip())
        params = PD(str(args[0]).strip(), str(args[1]).strip(), num_inserts=num_insertions)
    else:
        params = PD(str(args[0]).strip(), str(args[1]).strip())

    # Initialize the output directory
    params.output_dir()

    print("out_dir: %s" % params.out_dir)

    go(params)

    # convert_to_csvs.go(params.N_spheres(),params.final_step(),params.out_dir,'_n'+str(params.num_inserts))