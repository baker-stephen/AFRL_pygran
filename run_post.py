import sys

# import run_sim
import convert_to_csvs
import cy.run_post as rp

from param_defn import PD

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

    params.out_dir = 'outputs/1pt029/15_32/sim_out_9:15:15_16-8-2022'

    print("out_dir: %s" % params.out_dir)

    # Run the simulation
    # run_sim.go(params)

    wrt_dir = params.out_dir[:params.out_dir.rfind('/') + 1]

    print("write directory:", wrt_dir)

    sys.stdout = open(wrt_dir + 'stdout.txt', 'w')
    sys.stderr = open(wrt_dir + 'stderr.txt', 'w')

    print("write directory:", wrt_dir)

    print("\n\nSimulation Complete!\n\n")

    # Convert the output to Paraview-compatible and readable csvs
    convert_to_csvs.go(params.N_spheres(), params.final_step(), params.out_dir)  # TODO: test new csvs method

    print("\nAll Pygran steps complete")

    rp.go(params)

    print("\n\nRun Completed!!")
