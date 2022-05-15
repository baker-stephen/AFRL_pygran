from pygran import simulation
from pygran.params import steel, organic, glass
import numpy as np


if __name__ == "__main__":

    total_parts = 77
    num_insertions = 7
    parts_per_insert = total_parts//num_insertions
    # Create a dictionary of physical parameters
    params = {

        # Define the system
        'boundary': ('f', 'f', 'f'),  # fixed BCs
        'box': (-1, 1, -1, 1, -.2, 8),  # simulation box size

        # Define component(s)
        # Dp small = .25" = .00635m
        'species': (
            {'material': glass, 'style': 'sphere', 'radius': .125},),

        # Set skin distance to be 1/4 particle diameter
        # For dpsmall, 1/4 is .0015875
        'nns_skin': .03125,

        # Timestep
        'dt': 5e-6,

        # Apply gravitional force in the negative direction along the z-axis
        'gravity': (385.827, 0, 0, -1),

        # Setup I/O
        'traj': {'pfile': 'particles*.vtk', 'mfile': 'pipe*.vtk'},

        # Stage runs [optional]
        'stages': {'insertion': 1e5},

        # Define mesh for rotating mesh (tumbler)
        # TODO: define PVC material
        'mesh': {
            'pipe': {'file': 'mesh/pipe_pygran_half_bottom_1.stl', 'mtype': 'mesh/surface/stress', 'material': steel,
                     # 'args': {'scale': .0254}
                     },
        }
    }

    # Create an instance of the DEM class
    sim = simulation.DEM(**params)

    #Setup shaking:
    # freq = 40*2*np.pi
    # nTaps = 100
    # period = 1/freq
    # nSteps = period / params['dt']
    # amp = .1
    #
    # for i in range(nTaps):
    #     sim.moveMesh('pipe', viblin=('axis 0 0 1', 'order 1', 'amplitude {}'.format(amp), 'phase 0', 'period {}'.format(period)))
    #     sim.run(nSteps, params['dt'])
    #     sim.remove('moveMesh')

    # Insert 800 particles once in a cylinder
    # My best guess for cylinder numbers: x0, y0, r, z_min, z_max
    # pipe ID: 5/8" = .015875m, r = .0079375m, subtract sphere radius: .0015875, height: 4" = .1016m

    for i in range(num_insertions):
        insert = sim.insert(species=1, value=parts_per_insert, region=('cylinder', 'z', 0, 0, 0.1875, 4.25, 7.75),
                            args={'orientation': 'random'})
        # Add dissipative force proprtional to tablet velocity
        air_resistance = sim.addViscous(species=1, gamma=0.1)

        # Run insertion stage
        sim.run(params['stages']['insertion'], params['dt'])
        sim.remove(insert)

    #Setup shaking:
    freq = 40*2*np.pi
    nTaps = 100
    period = 1/freq
    nSteps = period / params['dt']
    ampz = .01
    ampxy = .005

    for i in range(nTaps):
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

        #vibrate z
        mm = sim.moveMesh('pipe', viblin=(
            'axis 0 0 1', 'order 1', 'amplitude {}'.format(ampz), 'phase 0', 'period {}'.format(period)))
        sim.run(nSteps, params['dt'])
        sim.remove(mm)

    #let settle

    sim.run(params['stages']['insertion'], params['dt'])