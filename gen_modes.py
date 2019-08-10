#!/usr/bin/env python

from __future__ import division
import os
import re
import sys
import numpy as np
import argparse as arg
from elements import ELEMENTS

au2ang = 0.5291771


def skiplines(openfile, nlines=0):
    '''Skips nlines + 1 lines in openfile. In other words, if nlines=0 it will
    go to the next line.'''

    for i in range(nlines):
        next(openfile)

    return next(openfile)


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''Generates .xyz trajectories for
                                each Normal Mode from a QM frequency calculation''',
                                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input files
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument('-f', '--filename',
                     default=None, type=str, dest="File", required=True,
                     help='''Excited State file''')

    inp.add_argument('-s', '--sel', default=None, nargs='+', type=str,
                     dest='AtomSel', help='''Atom Selection.''')

    inp.add_argument('--scale', default=5, type=int,
                     dest='Scale', help='''Scaling factor for displacement
                     vectors in the VMD script.''')

    #
    # Output files
    #
    out = parser.add_argument_group("Output Data")

    out.add_argument('-o', '--outdir', default="normal_modes",
                     type=str, dest="OutDir", help='''Output folder''')

    #
    # Parse and create the Options Dictionary
    #
    args = parser.parse_args()
    Opts = vars(args)

    if Opts['AtomSel']:
        Opts['AtomSel'] = read_sel(Opts['AtomSel'])

    return Opts


def extend_compact_list(idxs):

    extended = []

    # Uncomment this line if idxs is a string and not a list
    idxs = idxs.split()

    for idx in idxs:

        to_extend = idx.split('-')

        if len(to_extend) > 1:

            sel =  list(map(int, to_extend))
            extended += range(sel[0],sel[1]+1,1)

        else:
        
            extended.append(int(idx))
    
    return extended


def read_sel(string):

    string =  ','.join(string).replace(',,',',')

    try:
        f = open(string, 'r')
        string = f.readlines()
        f.close()
        string =  ','.join(string).replace(',,',',')
        string = string.replace(',', ' ')
        string = list(map(lambda x: x - 1, extend_compact_list(string)))

    except IOError:
        string = string.replace(',', ' ')
        string = list(map(lambda x: x - 1, extend_compact_list(string)))

    return string


def save_visdisps(coords, disps, sel=None, filename="structure"):

    header = ("#\n"
              "# VMD script to draw vectors\n"
              "#\n"
              "menu main on\n"
              "display projection orthographic\n"
              "display depthcue off\n"
              "display nearclip set 0.01\n"
              "axes location lowerleft\n"
              "\n"
              "#\n"
              "# VMD functions to draw a vector\n"
              "#\n"
              "proc vmd_draw_arrow {mol start end} {\n"
              "    set length [veclength [vecsub $end $start]]\n"
              "    set conelen [expr max(0.4,0.2*$length) ]\n"
              "    set scale [expr max(0.5,(1.0-$conelen/$length))]\n"
              "\n"
              "    set middle [vecadd $start [vecscale $scale [vecsub $end $start]]]\n"
              "    graphics $mol cylinder $start $middle radius 0.05\n"
              "    puts [list cone $middle $end radius 0.15]\n"
              "    graphics $mol cone $middle $end radius 0.15\n"
              "}\n"
              "\n"
              "proc vmd_draw_vector { mol pos val } {\n"
              "    set end   [ vecadd $pos [ vecscale +1 $val ] ]\n"
              "    vmd_draw_arrow $mol $pos $end\n"
              "}\n"
              "\n"
              "\n")

    #
    # Save the geometry and write the VMD script file with transition dips
    #
    with open('%s.vmd' % filename, 'w') as f:
        f.write(header)
        f.write("mol new equilibrium.xyz type xyz\n")
        f.write("mol new %s.xyz type xyz\n" % os.path.split(filename)[-1])

        if sel:
            f.write("mol rep Licorice 0.2 10\n")
            f.write("mol selection index %d to %d\n" % (sel[0], sel[-1]))
            f.write("mol addrep top\n")

        f.write("mol showrep 0 0 off\n")
        f.write("\n")

        for N, coord in enumerate(coords[sel][0]):
            disp = disps[N] * Opts['Scale']
            if np.linalg.norm(disp) > 0.00001 * Opts['Scale']:
                cmd = "graphics 0 color green; vmd_draw_vector 0 {%8.4f %8.4f %8.4f} {%8.4f %8.4f %8.4f}\n"
                data = [coord[0], coord[1], coord[2], disp[0], disp[1], disp[2]]
                f.write(cmd % tuple(data))

    return


def parse_molden(filename):

    with open(filename) as f:

        structure = []
        freqs = []
        modes = []
        forcecns = None
        redmasses = None

        for line in f:

            if "[FREQ]" in line:

                line = next(f)
                data = line.strip()
                while not data.startswith("["):
                    data = line.split()
                    freq = float(data[0])
                    freqs.append(freq)
                    line = next(f)
                    data = line.strip()

            if "[FR-COORD]" in line:

                line = next(f)
                data = line.strip()
                while not data.startswith("["):
                    data = line.strip().split()
                    data[1:] = list(map(float, data[1:]))
                    structure.append(data)
                    line = next(f)
                    data = line.strip()

            if "[FR-NORM-COORD]" in line:

                nat = len(structure)
                line = next(f)
                data = line.strip()
                i = 0
                while not data.startswith("["):

                    try:
                        if "vibration" in line:
                            i = 0
                            mode = []
                            line = next(f)

                        while i < nat:
                            data = line.strip().split()
                            data = list(map(float, data))
                            mode.append(data)
                            i += 1
                            line = next(f)
                            data = line.strip()

                        modes.append(mode)

                    except StopIteration:
                        while i < nat:
                            data = line.strip().split()
                            data = list(map(float, data))
                            mode.append(data)
                            i += 1
                            line = next(f)
                            data = line.strip()

                        modes.append(mode)

                        break


    Z_atoms = [ ELEMENTS[x[0]].number for x in structure ]
    masses = np.array([ ELEMENTS[x[0]].mass for x in structure ])
    atoms = [ x[0] for x in structure ]
    coords = np.array([x[1:] for x in structure]) #* au2ang
    freqs = np.array(freqs)
    modes = np.array(modes) #* au2ang

    return Z_atoms, masses, atoms, coords, freqs, forcecns, redmasses, modes


if __name__ == '__main__':

    Opts = options()
    Z_atoms, masses, atoms, coords, freqs, forcecns, redmasses, xdisps = parse_molden(Opts["File"])

    # #
    # # Choose the parser depending on the method and on the input files
    # #
    # parsertype = guess(Opts["File"])

    # if parsertype == "G09":
    #     freqparser = parsefreqs_G09

    # elif parsertype == "QChem":
    #     freqparser = parsefreqs_QChem

    # Z_atoms, masses, atoms, coords, freqs, forcecns, redmasses, xdisps = freqparser(Opts['File'])

    #
    # Create OutDir
    #
    if not os.path.exists(Opts['OutDir']):
        os.makedirs(Opts['OutDir'])

    #
    # Save Equilibrium Geometry
    #
    eqfile = os.path.join(Opts['OutDir'], "equilibrium")
    with open("%s.xyz" % eqfile, "w") as f:
        f.write("%d\n\n" % len(Z_atoms))
        np.savetxt(f, np.c_[Z_atoms, coords], fmt="%3d %12.6f %12.6f %12.6f")

    #
    # Cycle over Normal Modes
    #
    for N, xdisp in enumerate(xdisps, start=1):

        mode_anim = []

        #
        # Displace the eq geometry
        #
        for i in np.linspace(-1, 1, 20):

            coor_disp = coords + xdisp * i
            disp_geom = np.c_[Z_atoms, coor_disp]
            mode_anim.append(disp_geom)

        mode_anim = np.array(mode_anim)
        outfile = os.path.join(Opts['OutDir'], "mode_%03d_%d" % (N, freqs[N-1]))

        #
        # Save movie
        #
        with open("%s.xyz" % outfile, "w") as f:

            for step in mode_anim:
                f.write("%d\n\n" % len(Z_atoms))
                np.savetxt(f, step, fmt="%3d %12.6f %12.6f %12.6f")

            for step in mode_anim[::-1]:
                f.write("%d\n\n" % len(Z_atoms))
                np.savetxt(f, step, fmt="%3d %12.6f %12.6f %12.6f")

        save_visdisps(coords, xdisp, sel=Opts['AtomSel'], filename="%s" % outfile)
