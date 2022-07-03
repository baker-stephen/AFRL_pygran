from pygran import simulation
from pygran.params import steel, organic, glass
import numpy as np

#TODO: THIS IS SET UP FOR 4MM ID 1MM DP
#TODO: EVERYTHING IN MM

if __name__ == "__main__":

    total_parts = 300
    num_insertions = 1
    parts_per_insert = total_parts//num_insertions
    dp_in = 7/16
    # Create a dictionary of physical parameters
    params = {

        # Define the system
        'boundary': ('f', 'f', 'f'),  # fixed BCs
        # z bound given by funnel height (2"/sqrt(2) = 35.921mm) + pipe height (1.2"=30.48mm) + extra insertion room
        # x and y bounds given by funnel OR (funnel height) + pipe OR (.375"/2 = 4.7625mm) = 40.68mm
        #'box': (-22, 22, -22, 22, -3, 75),  # simulation box size mm
         'box': (-4, 4, -4, 4, -1.1, 17),  # simulation box size in inches
        # Define component(s)
        # Dp mini = 1mm, r = .5mm = .0198505"
        'species': (
            {'material': glass, 'style': 'sphere', 'radius': dp_in/2},), #TODO: change this for different Dp

        # Set skin distance to be 1/4 particle diameter
        # 'nns_skin': .25e-3,
        'nns_skin': dp_in/4, #TODO: change this for different Dp

        # Timestep
        # Needs to be reduced to satisfy rayleigh time constraint, apparently dependent on particle size
        'dt': 2.5e-7,

        # Apply gravitional force in the negative direction along the z-axis
        #'gravity': (9.81e3, 0, 0, -1),
         'gravity': (385.827, 0, 0, -1),

        # Setup I/O
        'traj': {'pfile': 'particles*.vtk', 'mfile': 'pipe*.vtk', 'freq': 50000},

        # Stage runs [optional]
        'stages': {'insertion': 2.3e6},

        # Define mesh for rotating mesh (tumbler)
        # Scale since stl is in inches
        # Used meshlab to reduce mesh count by .9, still needed curvature tolerance
        # TODO: define PVC material
        'mesh': {
            'pipe': {'file': 'mesh/three_half_in_pg.stl', 'mtype': 'mesh/surface/stress', 'material': steel,
                     'args': {'curvature_tolerant': 'yes'}
                     },
        },
    }

    print("sim")

    # Create an instance of the DEM class
    sim = simulation.DEM(**params)

    air_resistance = sim.addViscous(species=1, gamma=0.1)

    # particle insertion
    # My best guess for cylinder numbers: x0, y0, r, z_min, z_max
    # Minimum height of insert cylinder of rad .5": 12.2mm above top of pipe, or z = 43.18
    #TODO: change range back from 1 to num_insertions and value from 10 to parts_per_insert
    for i in range(num_insertions):
        # insert = sim.insert(species=1, value=parts_per_insert, region=('cylinder', 'z', 0, 0, 12.2e-3, 45e-3, 85e-3),
        #                     args={'orientation': 'random'})
        insert = sim.insert(species=1, value=parts_per_insert, region=('cylinder', 'z', 0, 0, 3, 12, 15),
                            args={'orientation': 'random'})

        # Run insertion stage
        sim.run(params['stages']['insertion']*2, params['dt'])
        sim.remove(insert)

    #allow for more settling after restart
    # sim.run(2e-7, params['dt'])

    #Setup shaking:
    freq = 10*2*np.pi
    nTaps = 20
    period = 1/freq
    nSteps = period / params['dt']
    ampz = .02
    ampxy = .015

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

    sim.run(params['stages']['insertion']*2, params['dt'])
