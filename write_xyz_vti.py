import numpy as np

if __name__ == "__main__":

    pipe_in_rad = 1.59 / 2

    x_res = 100
    x0 = -pipe_in_rad
    xM = pipe_in_rad
    dx = (xM - x0) / (x_res + 1)
    y_res = 100
    y0 = -pipe_in_rad
    yM = pipe_in_rad
    dy = (yM - y0) / (y_res + 1)
    z_res = 305
    z0 = 0.0
    zM = 9.5
    dz = (zM - z0) / (z_res + 1)
    loaded = []
    print("opening filled arr")
    with open('outputs/1pt59/1_4/filled_array_' + str(x_res) + '-' + str(y_res) + '-' + str(z_res) + '.npy', 'rb') as f:
        loaded = np.load(f)
        f.close()

    print("writing vti")
    with open('outputs/1pt59/1_4/xyz_vol_' + str(x_res) + '-' + str(y_res) + '-' + str(z_res) + '.vti', 'w') as f:
        f.write('<VTKFile type="ImageData">\n')
        f.write('<ImageData WholeExtent="'+str(0)+' '+str(x_res)+' '+str(0)+' '+str(y_res)+' '+str(0)+' '+str(z_res))
        f.write('" Origin="'+str(x0)+' '+str(y0)+' '+str(z0)+'" Spacing="'+str(dx)+' '+str(dy)+' '+str(dz)+'">\n')
        f.write('<Piece Extent="'+str(0)+' '+str(x_res)+' '+str(0)+' '+str(y_res)+' '+str(0)+' '+str(z_res)+'">\n')
        f.write('<PointData>\n')
        f.write('<DataArray Name="vol" NumberOfComponents="1" format="ascii" type="Float64">\n')
        for z in range(z_res+1):
            for y in range(y_res+1):
                for x in range(x_res+1):
                    f.write(str(loaded[x][y][z])+' ')
        f.write('\n</DataArray>\n')
        f.write('</PointData>\n')
        f.write('</Piece>\n')
        f.write('</ImageData>\n')
        f.write('</VTKFile>\n')
        f.close()
