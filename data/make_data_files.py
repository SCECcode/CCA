#!/usr/bin/env python

##
#  Retrieves the iteration number from intensity and builds the data
#  files in the expected format.
#

import getopt
import sys
import subprocess
import struct
import array

model = "CCA"
dimension_x = 1024
dimension_y = 896
dimension_z = 100

def usage():
    print("\n./make_data_files.py -i [iteration number]\n\n")
    print("-i - The iteration number to retrieve from intensity.\n\n")
    sys.exit(0)

def main():

    # Set our variable defaults.
    iteration = -1
    username = ""
    path = "/home/scec-01/enjuilee/work/" + model + "_ASCII"

    # Get the iteration number.
    try:
        opts, args = getopt.getopt(sys.argv[1:], "u:i:", ["user=", "iteration="])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)

    for o, a in opts:
        if o in ("-i", "--iteration"):
            iteration = str(a)
        if o in ("-u", "--user"):
            username = str(a) + "@"

    # If the iteration number was not provided, display the usage.
    if iteration == -1:
        usage()
        sys.exit(1)

    print("\nDownloading model file\n")

    subprocess.check_call(["scp", username +
                           "intensity.usc.edu:" + path + "/" + model + iteration.zfill(2) + ".ascii",
                           "."])

    # Now we need to go through the data files and put them in the correct
    # format for CCA. More specifically, we need a Vp.dat, Vs.dat, and Density.dat

    subprocess.check_call(["mkdir", "-p", "./i" + iteration.zfill(2)])

    print("\nWriting out CCA data files\n")

    f = open("./" + model + iteration.zfill(2) + ".ascii")

    f_vp = open("./i" + iteration.zfill(2) + "/vp.dat", "wb")
    f_vs = open("./i" + iteration.zfill(2) + "/vs.dat", "wb")
    f_density = open("./i" + iteration.zfill(2) + "/density.dat", "wb")

    vp_arr = array.array('f', (0,) * (dimension_x * dimension_y * dimension_z))
    vs_arr = array.array('f', (0,) * (dimension_x * dimension_y * dimension_z))
    density_arr = array.array('f', (0,) * (dimension_x * dimension_y * dimension_z))

    for line in f:
        arr = line.split()
        x_pos = int(arr[0]) - 1
        y_pos = int(arr[1]) - 1
        z_pos = int(arr[2]) - 1
        vp = float(arr[3])
        vs = float(arr[4])
        density = float(arr[5])

        vp_arr[z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos] = vp
        vs_arr[z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos] = vs
        density_arr[z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos] = density

    vp_arr.tofile(f_vp)
    vs_arr.tofile(f_vs)
    density_arr.tofile(f_density)

    f.close()
    f_vp.close()
    f_vs.close()
    f_density.close()

    print("\nDone!")

if __name__ == "__main__":
    main()
