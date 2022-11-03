from param_defn import PD

if __name__=="__main__":
    Dpstr = "1/4"
    IDstr = "0.602"
    params = PD(IDstr,Dpstr)
    params.output_dir()
    params.out_dir = params.out_dir[:params.output_dir().rfind('/')+1]
    import_csv = params.out_dir+'particles.csv'
    print(import_csv)

    print(params.z0_cutoff())
    print(params.zM_cutoff())
    slice_mid = (params.zM_cutoff()+params.z0_cutoff())/2
    slice_range = (slice_mid-params.DP, slice_mid+params.DP)
    print(slice_range)
    positions = []

    with open(import_csv, 'r') as r:
        r.readline()
        for line in r:
            xyz = [float(p) for p in line.split(",")]
            if xyz[2]<slice_range[0] or xyz[2]>slice_range[1]:
                #Outside our slice
                continue
            positions.append(xyz)
        r.close()
    positions = sorted(positions, key=lambda x: x[2])
    print(len(positions))
    print(positions)
    low = positions[0][2]
    for pos in positions:
        pos[2] -= low
        pos[2] = round(pos[2],7)
        print(pos)