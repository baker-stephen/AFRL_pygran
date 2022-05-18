from pygran import simulation
from pygran.params import steel, organic, glass
import numpy as np

#TODO: mesh does not like something about mesh angel and curvature. Try a chamfer.

if __name__ == "__main__":

    total_parts = 1000
    num_insertions = 10
    parts_per_insert = total_parts//num_insertions
    # Create a dictionary of physical parameters
    params = {

        # Define the system
        'boundary': ('f', 'f', 'f'),  # fixed BCs
        # z bound given by funnel height (2"/sqrt(2) = 35.921mm) + pipe height (1.2"=30.48mm) + extra insertion room
        # x and y bounds given by funnel OR (funnel height) + pipe OR (.375"/2 = 4.7625mm) = 40.68mm
        'box': (-50e-3, 50e-3, -50e-3, 50e-3, -10e-3, 90e-3),  # simulation box size
        # 'box': (-2, 2, -2, 2, -.2, 4),  # simulation box size in inches
        # Define component(s)
        # Dp mini = 1mm, r = .5mm = .0198505"
        'species': (
            {'material': glass, 'style': 'sphere', 'radius': .5e-3},),

        # Set skin distance to be 1/4 particle diameter
        # 'nns_skin': .25e-3,
        'nns_skin': (.5e-3)/2,

        # Timestep
        'dt': 5e-6,

        # Apply gravitional force in the negative direction along the z-axis
        'gravity': (9.81, 0, 0, -1),
        # 'gravity': (385.827, 0, 0, -1),

        # Setup I/O
        'traj': {'pfile': 'particles*.vtk', 'mfile': 'pipe*.vtk'},

        # Stage runs [optional]
        'stages': {'insertion': 1e5},

        # Define mesh for rotating mesh (tumbler)
        # Scale since stl is in inches
        # TODO: define PVC material
        'mesh': {
            'pipe': {'file': 'mesh/one_quarter_funnel_blend.stl', 'mtype': 'mesh/surface/stress', 'material': steel,
                     'args': {'scale': .0254,}
                     },
        },
    }

    print("sim")

    # Create an instance of the DEM class
    sim = simulation.DEM(**params)

    print("command")

    sim.command("fix cad all mesh/surface file /home/stephen/pg/DEM_tablet_ex/mesh/one_quarter_funnel_blend.stl type 1")

    print("after")

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

    # particle insertion
    # My best guess for cylinder numbers: x0, y0, r, z_min, z_max
    # Minimum height of insert cylinder of rad .5": 12.2mm above top of pipe, or z = 43.18
    for i in range(num_insertions):
        insert = sim.insert(species=1, value=parts_per_insert, region=('cylinder', 'z', 0, 0, 12.2e-3, 45e-3, 85e-3),
                            args={'orientation': 'random'})
        # insert = sim.insert(species=1, value=parts_per_insert, region=('cylinder', 'z', 0, 0, 12.2/25.4, 45/25.4, 85/25.4),
        #                     args={'orientation': 'random'})
        # Add dissipative force proprtional to tablet velocity
        air_resistance = sim.addViscous(species=1, gamma=0.1)

        # Run insertion stage
        sim.run(params['stages']['insertion'], params['dt'])
        sim.remove(insert)

    #Setup shaking:
    # freq = 40*2*np.pi
    # nTaps = 100
    # period = 1/freq
    # nSteps = period / params['dt']
    # ampz = .01
    # ampxy = .005
    #
    # for i in range(nTaps):
    #     #vibrate x
    #     mm = sim.moveMesh('pipe', viblin=(
    #         'axis 1 0 0', 'order 1', 'amplitude {}'.format(ampxy), 'phase 0', 'period {}'.format(period)))
    #     sim.run(nSteps, params['dt'])
    #     sim.remove(mm)
    #
    #     #vibrate y
    #     mm = sim.moveMesh('pipe', viblin=(
    #         'axis 0 1 0', 'order 1', 'amplitude {}'.format(ampxy), 'phase 0', 'period {}'.format(period)))
    #     sim.run(nSteps, params['dt'])
    #     sim.remove(mm)
    #
    #     #vibrate z
    #     mm = sim.moveMesh('pipe', viblin=(
    #         'axis 0 0 1', 'order 1', 'amplitude {}'.format(ampz), 'phase 0', 'period {}'.format(period)))
    #     sim.run(nSteps, params['dt'])
    #     sim.remove(mm)
    #
    # #let settle
    #
    # sim.run(params['stages']['insertion'], params['dt'])