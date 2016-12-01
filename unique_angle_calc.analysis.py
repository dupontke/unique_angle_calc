#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
# ----------------------------------------
# USAGE:

# ./pi-pi_stacking_calc.analysis.py pi-pi_stacking_calc.config 

# ----------------------------------------
# PREAMBLE:

import scipy
import sys
import os
import numpy as np
import math
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *
from sel_list import *

# constant to convert radians to degrees
radians_to_degrees = 180.0/3.1415926535

flush = sys.stdout.flush

config_file = sys.argv[1]

# ----------------------------------------
# FUNCTIONS:

def ffprint(string):
    print '%s' %(string)
    flush()

necessary_parameters = ['pdb','traj_loc','start','end','Wrapped','outputname_ang','selection_file']
all_parameters = ['pdb','traj_loc','start','end','Wrapped','outputname_ang','selection_file','important','write_summary','summary_filename','selection_output']
def config_parser(config_file):# Function to take config file and create/fill the parameter dictionary
    for i in range(len(necessary_parameters)):
        parameters[necessary_parameters[i]] = ''

        # SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['important'] = None
        parameters['write_summary'] = False
        parameters['summary_filename'] = 'unique_angle_calc.summary'
        parameters['summary_filename'] = 'selections.txt'

        # GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
        execfile(config_file,parameters)
        for key, value in parameters.iteritems():
            if value == '':
                print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
                sys.exit()

def summary():
    with open('%s' %(parameters['summary_filename']),'w') as f:
        f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
        f.write('To recreate this analysis, run this line:\n')
        for i in range(len(sys.argv)):
            f.write('%s ' %(sys.argv[i]))
        f.write('\n\n')
        f.write('Parameters used:\n')
        for i in all_parameters:
            f.write('%s = %s \n' %(i,parameters[i]))
        f.write('\n\n')
        f.write('output is written to:\n')
        f.write('%s\n' %(parameters['outputname_ang']))
        f.write('\nNumber of steps analyzed: %d\n' %(nSteps))
        f.write('\nAtom selections analyzed have been written out to %s\n' %(parameters['selection_output']))

# ----------------------------------------
# MAIN PROGRAM:
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# initiate coordinate universe
u = MDAnalysis.Universe(parameters['pdb'])
u_important = u.select_atoms(parameters['important'])
if not parameters['Wrapped']:     # Test to see if the 'Wrapped' key is equal to False
    rest = u.select_atoms('not (resname WAT or resname Na+ or resname Cl- or protein)')
    rest_nRes = rest.n_residues
    print rest_nRes


# atom selections for unique angle calculation
nSel = len(sel)
u_sel = []
for i in range(nSel):
    sel1 = sel[i][0]
    sel2 = sel[i][1]
    sel3 = sel[i][2]
    u_sel.append([u.select_atoms(sel1),u.select_atoms(sel2),u.select_atoms(sel3)])

nSteps = 0
start = int(parameters['start'])
end = int(parameters['end'])
with open(parameters['outputname_ang'],'w') as f:
    ffprint('Beginning trajectory analysis')
    while start <= end:
        ffprint('Loading trajectory analysis')
        u.load_new('%sproduction.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))
        nSteps += len(u.trajectory)
        # Loop through trajectory
        for ts in u.trajectory:
            # unnecessary calculations if the trajectory is wrapped.
            if not parameters['Wrapped']: 
                dims = u.dimensions[:3]
                dims2 = dims/2.0

                for i in range(rest_nRes):
                    COM = rest.residues[i].center_of_mass()
                    t = wrapping(COM,dims,dims2)
                    rest.residues[i].atoms.translate(t)

            # get the three atom positions
            for i in range(nSel):
                atom1_pos = u_sel[i][0][0].position
                atom2_pos = u_sel[i][1][0].position
                atom3_pos = u_sel[i][2][0].position
                # determine two vectors using the three atom positions
                r12 = atom1_pos-atom2_pos
                r12 /= np.linalg.norm(r12)
                r32 = atom3_pos-atom2_pos
                r32 /= np.linalg.norm(r32)
                r12r32 = np.dot(r12,r32)
                angle = np.arccos(r12r32)*radians_to_degrees
                f.write('%10.6f     ' %(angle))
            f.write('\n')
        start +=1

if parameters['write_summary']:
    summary()

