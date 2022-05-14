from pygran import simulation
from pygran.params import steel, organic


if __name__ == "__main__":
    # Create a dictionary of physical parameters
    params = {

        # Define the system
        'boundary': ('f', 'f', 'f'),  # fixed BCs
        'box': (-1, 1, -1, 1, -.2, 5),  # simulation box size

        # Define component(s)
        # Dp small = .25" = .00635m
        'species': (
            {'material': organic, 'style': 'sphere', 'radius': .125},),

        # Set skin distance to be 1/4 particle diameter
        # For dpsmall, 1/4 is .0015875
        'nns_skin': .03125,

        # Timestep
        'dt': 1e-6,

        # Apply gravitional force in the negative direction along the z-axis
        'gravity': (385.827, 0, -1, 0),

        # Setup I/O
        'traj': {'pfile': 'particles*.vtk', 'mfile': 'pipe*.vtk'},

        # Stage runs [optional]
        'stages': {'insertion': 2e6},

        # Define mesh for rotating mesh (tumbler)
        # TODO: define PVC material
        'mesh': {
            'pipe': {'file': 'mesh/pipe_pygran_half_bottom.stl', 'mtype': 'mesh/surface/stress', 'material': steel,
                     # 'args': {'scale': .0254}
                     },
        }
    }

    # Create an instance of the DEM class
    sim = simulation.DEM(**params)

    # Insert 800 particles once in a cylinder
    # My best guess for cylinder numbers: x0, y0, r, z_min, z_max
    # pipe ID: 5/8" = .015875m, r = .0079375m, subtract sphere radius: .0015875, height: 4" = .1016m
    insert = sim.insert(species=1, value=15, region=('cylinder', 'z', 0, 0, 0.2875, .3, 4.7),
                        args={'orientation': 'random'})

    # Add dissipative force proprtional to tablet velocity
    air_resistance = sim.addViscous(species=1, gamma=0.1)

    # Run insertion stage
    sim.run(params['stages']['insertion'], params['dt'])