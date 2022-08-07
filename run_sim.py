from pygran import simulation
from pygran.params import steel, glass
import numpy as np

from init_sim import PD


def go(sim_params: PD):
    parts_per_insert = sim_params.N_spheres() // sim_params.num_inserts
    print('Running ' + sim_params.ID_str + '" ID with atom count ' + str(
        sim_params.N_spheres()) + ', and ' + sim_params.DP_str + '" DP_str.\n')

    # Create a dictionary of physical parameters
    params = {

        # Define the system
        'boundary': ('f', 'f', 'f'),  # fixed BCs

        'box': sim_params.bounds(),  # simulation box size in inches

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
        'traj': {'pfile': 'particles*.vtk', 'mfile': 'pipe*.vtk', 'freq': 10000},

        'output': sim_params.output_dir(),

        # Stage runs
        'stages': {'insertion': 2.3e6 / sim_params.num_inserts},

        # Define mesh
        'mesh': {
            'pipe': {'file': sim_params.mesh(), 'mtype': 'mesh/surface/stress', 'material': steel,
                     'args': {'curvature_tolerant': 'yes'}
                     },
        },
    }

    # Create an instance of the DEM class
    sim = simulation.DEM(**params)

    #Add a dissipative force
    air_resistance = sim.addViscous(species=1, gamma=0.1)

    for i in range(sim_params.num_inserts):
        # Insert a group of particles
        insert = sim.insert(species=1, value=parts_per_insert, region=sim_params.insert(),
                            args={'orientation': 'random'})

        # Run insertion stage, let the particles settle into the cylinder
        sim.run(params['stages']['insertion'] * 2, params['dt'])
        sim.remove(insert)

        # Setup shaking:
        freq = 10 * 2 * np.pi
        nTaps = 30
        period = 1 / freq
        nSteps = period / params['dt']
        ampz = sim_params.ampz()
        ampxy = sim_params.ampxy()

        for i in range(nTaps // 2):
            # vibrate x
            mm = sim.moveMesh('pipe', viblin=(
                'axis 1 0 0', 'order 1', 'amplitude {}'.format(ampxy), 'phase 0', 'period {}'.format(period)))
            sim.run(nSteps, params['dt'])
            sim.remove(mm)
            # vibrate y
            mm = sim.moveMesh('pipe', viblin=(
                'axis 0 1 0', 'order 1', 'amplitude {}'.format(ampxy), 'phase 0', 'period {}'.format(period)))
            sim.run(nSteps, params['dt'])
            sim.remove(mm)

        # Allow for some settling
        sim.run(params['stages']['insertion'] / 2, params['dt'])

        for i in range(nTaps):
            # vibrate z
            mm = sim.moveMesh('pipe', viblin=(
                'axis 0 0 1', 'order 1', 'amplitude {}'.format(ampz), 'phase 0', 'period {}'.format(period)))
            sim.run(nSteps, params['dt'])
            sim.remove(mm)

        # Let simulation settle before next insertion

        sim.run(params['stages']['insertion'], params['dt'])
