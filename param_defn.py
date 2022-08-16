import numpy as np
from datetime import datetime as dt
import os


class PD:
    """
    PD is short for 'Parameter Dictionaries'.
    Initialize a PD object with an inner diameter, and call the various methods to retrieve the proper parameter.
    Make sure to add every parameter when adding a new ID.
    """
    # Length of the cylinder to fill
    lens = {'0.26': 1.0, '0.602': 12.0, '1.029': 5.0, '1.59': 10.0}
    # Boundaries of the simulation: (-x,+x,-y,+y,-z,+z). Be sure to provide extra space for shaking.
    boxes = {'0.26': (-0.6, 0.6, -0.6, 0.6, -1.5, 6),
             '0.602': (-6, 6, -6, 6, -2, 25),
             '1.029': (-4.5, 4.5, -4.5, 4.5, -.5, 15),
             '1.59': (-4, 4, -4, 4, -1.1, 17)}
    # The name of the mesh file
    meshes = {'0.26': 'one_quarter_funnel_closed.stl',
              '0.602': 'half_pg_taller_blended_funnel_reduced.stl',
              '1.029': 'one_in_pg_bl_reduced.stl',
              '1.59': 'three_half_in_pg_reduced_smooth.stl'}
    # The shape and dimension of the group of inserted particles.
    # Format: (shape, axial direction, x position of center axis, y position of center axis, radius, z minimum, z max)
    # Be sure to give a bit of room between the mesh walls and this cylinder.
    inserts = {'0.26': ('cylinder', 'z', 0, 0, .4, 2, 3),
               '0.602': ('cylinder', 'z', 0, 0, 5.5, 15, 20),
               '1.029': ('cylinder', 'z', 0, 0, 3.8, 9, 13.5),
               '1.59': ('cylinder', 'z', 0, 0, 3, 12, 15)}
    # The amplitude of the shaking in the z direction.
    # Not an exact science, generally want enough shaking so there's a good amount of movement,
    # but not so much that the particles violently fly upward.
    ampzs = {'0.26': 0.02, '0.602': 0.022, '1.029': 0.025, '1.59': 0.025}
    # The amplitude of the shaking in the x and y directions.
    # Same applies here.
    ampxys = {'0.26': 0.015, '0.602': 0.017, '1.029': 0.02, '1.59': 0.02}
    # The first two values are the minimum and maximum z values where spheres could exist with a valid inner diameter.
    # The second two values are the min and max z values where the spheres have settled into a more ordered state.
    post_zs = {'0.26': [0.15, 1.0, 0.3, 0.95],
               '1.029': [0.5, 4.5, 1.0, 4.5],
               '1.59': [0.0, 9.5, 1.0, 9.5],
               '0.602': [0.0, 10.0, 2.0, 10.0]}

    def __init__(self, ID: str, DP: str, num_inserts=1):
        self.ID_str = ID
        self.DP_str = DP
        # Inner diameter is given simply as a decimal
        self.ID = float(ID)
        # All sphere diameters are provided as a fraction in inches.
        Dp_frac = [float(d) for d in DP.split('/')]
        self.DP = Dp_frac[0] / Dp_frac[1]
        # Number of insertion phases to run.
        # Increasing this number greatly increases the runtime, but can yield more optimal packing.
        self.num_inserts = num_inserts
        self.out_dir = "not_created"

    def length(self) -> float:
        return self.lens[self.ID_str]

    def bounds(self) -> tuple:
        return self.boxes[self.ID_str]

    def mesh(self) -> str:
        return 'mesh/' + self.meshes[self.ID_str]

    def insert(self) -> tuple:
        return self.inserts[self.ID_str]

    def ampz(self) -> float:
        return self.ampzs[self.ID_str]

    def ampxy(self) -> float:
        return self.ampxys[self.ID_str]

    def output_dir(self) -> str:
        """
        The output directory structure follows: outputs/[ID]/[DP_str],
        where '.' are replaced with 'pt', and fractions (or '/') are replaced with '_'
        :return: A time-stamped directory at which the simulation outputs are saved
        """

        # If the output directory has already been created, simply return that
        if self.out_dir != "not_created":
            return self.out_dir
        out_dir = 'outputs/'
        out_dir += self.ID_str.replace('.', 'pt') + '/'
        out_dir += self.DP_str.replace('.', 'pt').replace('/', '_') + '/'
        # Create this directory if it does not exist
        try:
            os.makedirs(out_dir, exist_ok=True)
            print("Directory '%s' created successfully\n" % out_dir)
        except OSError as error:
            print("Directory '%s' can not be created. Error: %s\n" % (out_dir, str(error)))
        time = dt.now()
        out_dir += 'sim_out_{}:{}:{}_{}-{}-{}'.format(
            time.hour,
            time.minute,
            time.second,
            time.day,
            time.month,
            time.year,
        )
        self.out_dir = out_dir
        return out_dir

    def porosity_Guo(self) -> float:
        """
        :param N: Ratio of ID to DP_str
        :return: porosity analytically determined in paper by Guo
        """
        N = self.ID / self.DP
        if N >= 3.0 or 1.866 < N < 2.0:
            # print("out of range for Guo")
            return -1

        k1 = 1
        k2 = 1

        if N <= 1.866:
            k1 = N ** (-2)
            k2 = 1 / (np.sqrt(1 - (N - 1) ** 2))
        elif 2 <= N <= 3:
            n = 2
            alpha = np.pi / 4
            if 2.1547 <= N < 2.4142:
                n = 3
                alpha = np.pi / 6
            elif 2.4142 <= N < 2.7013:
                n = 4
                alpha = 22.5 * (np.pi / 180)
            elif 2.7013 <= N:
                n = 5
                alpha = 18 * (np.pi / 180)
            k1 = n / (N ** 2)
            k2 = 1 / np.sqrt(1 - ((N - 1) * np.sin(alpha)) ** 2)

        return 1 - (2 / 3) * k1 * k2

    def poros_Foumeny(self) -> float:
        """
        :param N: Ratio of ID to DP_str
        :return: porosity empirically determined by Foumeny
        """
        N = self.ID / self.DP
        return .383 + .254 * np.power(N, -0.923) * np.power(.723 * N - 1, -0.5)

    def poros_Cheng(self) -> float:
        """
        :param N: Ratio of ID to DP_str
        :return: porosity empirically determined by Cheng
        """
        N = self.ID / self.DP
        eps_1 = .8 * (N - 1) ** .27
        eps_2 = .38 * (1 + (1 / (N - 1)) ** 1.9)
        return (eps_1 ** (-3) + eps_2 ** (-3)) ** (-1 / 3)

    def N_spheres(self, buffer=.1) -> int:
        """
        :param ID: Inner diameter of the cylinder, inches
        :param DP: Diameter of the spheres/pellets, inches
        :param len: Length of the cylinder to fill
        :param buffer: Percentage above the estimated number of spheres that will exactly fill the cylinder
        :return: The number of spheres needed to fill the cylinder
        """
        N = self.ID / self.DP
        # Guo is the most accurate predictor, but only applicable in a small range
        if N < 3.0 and not (1.866 < N < 2.0):
            porosity = self.porosity_Guo()
        # Both Foumeny and Cheng have similar predictions, they are almost interchangeable
        else:
            porosity = self.poros_Foumeny()
        vol_to_fill = (np.pi * (self.ID / 2) ** 2 * self.length()) * (1 - porosity)
        vol_sphere = (4 / 3) * np.pi * (self.DP / 2) ** 3
        N_spheres = vol_to_fill / vol_sphere
        return int(N_spheres * (1 + buffer))

    def final_step(self) -> int:
        # Total number of steps
        settle_steps = 2.3e6 / self.num_inserts
        dt = 2.5e-7
        freq = 10 * 2 * np.pi
        period = 1 / freq
        shake_steps = period / dt
        # shake_steps = (1 / 20 * np.pi) / 2.5e-7

        # nsteps = period / dt

        # print((2* settle_steps + 15*shake_steps*2 + .5*settle_steps + 30*shake_steps + settle_steps)*self.num_inserts)

        return int((3.5 * settle_steps + 60 * shake_steps) * self.num_inserts)

    def z0(self) -> float:
        return self.post_zs[self.DP_str][0]

    def zM(self) -> float:
        return self.post_zs[self.DP_str][1]

    def z0_cutoff(self) -> float:
        return self.post_zs[self.DP_str][2]

    def zM_cutoff(self) -> float:
        return self.post_zs[self.DP_str][3]

    def res(self) -> (int,int,int):
        """
        :return: x, y, z, and resolution to make into a 3D grid, and then 1D radial average.
        I have found these values to work best at all ID and DP in my testing.
        """
        return (800,800,2500,200)
