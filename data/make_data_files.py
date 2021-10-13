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

if sys.version_info.major >= (3) :
  from urllib.request import urlopen
else:
  from urllib2 import urlopen


model = "CCA"
dimension_x = 0
dimension_y = 0
dimension_z = 0

def usage():
    print("\n./make_data_files.py -i [6]\n\n")
    print("-i - The iteration number to retrieve\n")
    sys.exit(0)

def download_urlfile(url,fname):
  try:
    response = urlopen(url)
    CHUNK = 16 * 1024
    with open(fname, 'wb') as f:
      while True:
        chunk = response.read(CHUNK)
        if not chunk:
          break
        f.write(chunk)
  except:
    e = sys.exc_info()[0]
    print("Exception retrieving and saving model datafiles:",e)
    raise
  return True

def main():

    # Set our variable defaults.
    iteration = -1
    path = ""

    try:
        fp = open('./config','r')
    except:
        print("ERROR: failed to open config file")
        sys.exit(1)
    
    ## look for model_data_path and other varaibles
    lines = fp.readlines()
    for line in lines :
        if line[0] == '#' :
          continue
        parts = line.split('=')
        if len(parts) < 2 :
          continue;
        variable=parts[0].strip()
        val=parts[1].strip()
       
        if (variable == 'model_data_path') : 
            path = val + '/' + model
            continue
        if (variable == 'nx') : 
            dimension_x = int(val)
            continue
        if (variable == 'ny') : 
            dimension_y = int(val)
            continue
        if (variable == 'nz') : 
            dimension_z = int(val)
            continue
        
        continue
    if path == "" :
        print("ERROR: failed to find variables in config file")
        sys.exit(1)

    fp.close()

    # Get the iteration number.
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:", ["iteration="])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)

    for o, a in opts:
        if o in ("-i", "--iteration"):
            iteration = str(a)

    # If the iteration number was not provided, display the usage.
    if iteration == -1:
        usage()
        sys.exit(1)

    print("\nDownloading model file\n")

    fname= model + iteration.zfill(2) + ".ascii"
    url = path + "/" + fname
    download_urlfile(url,fname)

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
