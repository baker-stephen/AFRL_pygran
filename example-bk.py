from pygran import simulation
from pygran.params import organic


if __name__ == "__main__":
    # Create a dictionary of physical parameters
    params = {

        # Define the system
        'boundary': ('f', 'f', 'f'),  # fixed BCs
        'box': (-1, 1, -1, 1, -1, 1),  # simulation box size

        # Define component(s)
        'species': (
            {'material': organic, 'style': 'multisphere/tablet', 'radius': 2e-2, 'nspheres': 12, 'length': 1e-1},),

        # Set skin distance to be 1/4 particle diameter
        'nns_skin': 5e-3,

        # Timestep
        'dt': 5e-6,

        # Apply gravitional force in the negative direction along the z-axis
        'gravity': (9.81, 0, 0, -1),

        # Setup I/O
        'traj': {'pfile': 'particles*.vtk', 'mfile': 'tumbler*.vtk'},

        # Stage runs [optional]
        'stages': {'insertion': 1e6, 'rotation': 1e6},

        # Define mesh for rotating mesh (tumbler)
        'mesh': {
            'tumbler': {'file': 'mesh/tumbler.stl', 'mtype': 'mesh/surface/stress', 'material': organic,
                        'args': {'scale': 1e-3}},
        }
    }

    # Create an instance of the DEM class
    sim = simulation.DEM(**params)

    # Insert 800 particles once in a cylinder
    insert = sim.insert(species=1, value=800, region=('cylinder', 'y', 0, 0, 0.7, -0.4, 0.4),
                        args={'orientation': 'random'})

    # Add dissipative force proprtional to tablet velocity
    air_resistance = sim.addViscous(species=1, gamma=0.1)

    # Run insertion stage
    sim.run(params['stages']['insertion'], params['dt'])

    # Delete insertion fix
    sim.remove(insert)

    # Rotate mesh (tumbler) along the xoz plane
    sim.moveMesh(name='tumbler', rotate=('origin', 0, 0, 0), axis=(0, 1, 0), period=5e-1)

    # Run rotation stage
    sim.run(params['stages']['rotation'], params['dt'])