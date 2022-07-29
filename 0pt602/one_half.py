from pygran import simulation
from pygran.params import steel, organic, glass
import numpy as np
from datetime import datetime as dt

import sys

if __name__ == "__main__":

    args = sys.argv[1:]
    print('args:', args)
    if len(args) != 1:
        raise Exception("Please specify only DP")

    DP_N_dict = {'1/4':220,'1/8':2200,'1/16':16500,'3/16':600,'7/16':50}

    dp_str = str(args[0]).strip()
    total_parts = DP_N_dict[dp_str]
    num_insertions = 5
    parts_per_insert = total_parts // num_insertions
    dp_frac = [float(x) for x in dp_str.split('/')]
    dp_in = dp_frac[0]/dp_frac[1]
    print('Running 0.602" ID with atom count '+str(total_parts)+', and '+str(dp_in)+'" DP.\n')
    time = dt.now()
    # Create a dictionary of physical parameters
    params = {

        # Define the system
        'boundary': ('f', 'f', 'f'),  # fixed BCs
        # z bound given by funnel height (2"/sqrt(2) = 35.921mm) + pipe height (1.2"=30.48mm) + extra insertion room
        # x and y bounds given by funnel OR (funnel height) + pipe OR (.375"/2 = 4.7625mm) = 40.68mm
        #'box': (-22, 22, -22, 22, -3, 75),  # simulation box size mm
        'box': (-6, 6, -6, 6, -2, 25),  # simulation box size in inches
        # Define component(s)
        # Dp mini = 1mm, r = .5mm = .0198505"
        'species': (
            {'material': glass, 'style': 'sphere', 'radius': dp_in/2},),

        # Set skin distance to be 1/4 particle diameter
        # 'nns_skin': .25e-3,
        'nns_skin': dp_in/4,

        # Timestep
        # Needs to be reduced to satisfy rayleigh time constraint, apparently dependent on particle size
        'dt': 2.5e-7,

        # Apply gravitional force in the negative direction along the z-axis
        #'gravity': (9.81e3, 0, 0, -1),
         'gravity': (385.827, 0, 0, -1),

        # Setup I/O
        'traj': {'pfile': 'particles*.vtk', 'mfile': 'pipe*.vtk', 'freq': 10000},

        'output': dp_str.replace('/','_').replace('.','pt')+'/{}-{}-{}_{}-{}-{}'.format(
                    time.hour,
                    time.minute,
                    time.second,
                    time.day,
                    time.month,
                    time.year,
                ),

        # Stage runs [optional]
        'stages': {'insertion': 2.3e6/num_insertions},

        # Define mesh for rotating mesh (tumbler)
        # Scale since stl is in inches
        # Used meshlab to reduce mesh count by .9, still needed curvature tolerance
        # TODO: define PVC material
        'mesh': {
            'pipe': {'file': '../mesh/half_pg_taller_blended_funnel_reduced.stl', 'mtype': 'mesh/surface/stress', 'material': steel,
                     'args': {'curvature_tolerant': 'yes'}
                     },
        },
    }

    print("sim")

    # Create an instance of the DEM class
    sim = simulation.DEM(**params)


    air_resistance = sim.addViscous(species=1, gamma=0.1)


    for i in range(num_insertions):
        # insert = sim.insert(species=1, value=parts_per_insert, region=('cylinder', 'z', 0, 0, 12.2e-3, 45e-3, 85e-3),
        #                     args={'orientation': 'random'})
        insert = sim.insert(species=1, value=parts_per_insert, region=('cylinder', 'z', 0, 0, 5.5, 15, 20),
                            args={'orientation': 'random'})


        # Run insertion stage
        sim.run(params['stages']['insertion']*2, params['dt'])
        sim.remove(insert)


        #Setup shaking:
        freq = 10*2*np.pi
        nTaps = 30
        period = 1/freq
        nSteps = period / params['dt']
        ampz = .022
        ampxy = .017

        for i in range(nTaps//2):
            #vibrate x
            mm = sim.moveMesh('pipe', viblin=(
                'axis 1 0 0', 'order 1', 'amplitude {}'.format(ampxy), 'phase 0', 'period {}'.format(period)))
            sim.run(nSteps, params['dt'])
            sim.remove(mm)
            #vibrate y
            mm = sim.moveMesh('pipe', viblin=(
                'axis 0 1 0', 'order 1', 'amplitude {}'.format(ampxy), 'phase 0', 'period {}'.format(period)))
            sim.run(nSteps, params['dt'])
            sim.remove(mm)

        sim.run(params['stages']['insertion']/2, params['dt'])

        for i in range(nTaps):
            #vibrate z
            mm = sim.moveMesh('pipe', viblin=(
                'axis 0 0 1', 'order 1', 'amplitude {}'.format(ampz), 'phase 0', 'period {}'.format(period)))
            sim.run(nSteps, params['dt'])
            sim.remove(mm)

        #let settle

        sim.run(params['stages']['insertion'], params['dt'])
