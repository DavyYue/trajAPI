
# coding: utf-8

# # testAPI_automate
#
# Created by Davy Yue 2017-06-26

## 2017_06_27 1:00 PM CST
## **** DEBUGGING/TESTING NOTES ****
## - think about inputs and outputs
## - think about where variables can be passed as parameters
## - need test cases - maybe with n-dodecane
## - need more filetypes to read instead of just hoomdxml
## - check with generated rdfs from testAPI_propane to compare results
## - maybe use more functions to call
## - consider adaptation for other elements (any periodic ones)
##       use element dictionary in msibi repo


# ### Imports
import itertools
import string
import os
from xml.etree import cElementTree as ET

import parmed as pmd
from periodic import element
import numpy as np
import matplotlib.pyplot as plt
# get_ipython().magic('matplotlib inline')

from msibi import MSIBI, State, Pair, mie
import mdtraj as md

## check possible other packages to ease process


# ### Functions
def read_search_mapping(search_mapping_filename, user_mapping_filename):
    # parser = ET.XMLParser(encoding="utf-8")
    # root = ET.fromstring(search_mapping_filename, parser=parser)
    root = ET.fromstring(open(search_mapping_filename).read())
    # root = ET.parse(search_mapping_filename).getroot()
    searchstr = [] # list containing all search values ordered by priority
    for value in root.findall('value'):
        searchstr.append(value.attrib['searchstr'])
        print(searchstr)

# read_search_mapping('propane_search_mapping.xml', 'propane_user_mapping.xml')

# SMARTS string
# Current supported SMARTS string: https://github.com/mosdef-hub/foyer/issues/63
# https://pubchem.ncbi.nlm.nih.gov/search/help_search.html

def read_user_mapping(user_mapping_filename):
    # parser = ET.XMLParser(encoding="utf-8")
    # root = ET.fromstring(user_mapping_filename, parser=parser)
    root = ET.fromstring(open(user_mapping_filename).read())
    # root = ET.parse(user_mapping_filename).getroot()

    # Get molecule_name
    bead = root.find('bead')
    molecule_name = bead.attrib["molecule_name"]

    # Get element_names, n_unitsPerSection
    element_names = [] # only need for atom definition
    n_sections_BEAD = 0
    n_unitsPerSection = 0 # number units per section being combined
    for section in bead.findall('section'):
        # allow for different number of units in each section per bead
        # Element identification using Periodic python package
        # https://github.com/luisnaranjo733/periodic
        # use foyer here to find indices

        n_unitsPerSection = elements_symb.count("C")
        # need to modify for different center elements
        # counts number of carbon center atoms since each section is organic


    print(n_unitsPerSection)
    # check later for more different elements with loop

    return n_unitsPerSection, molecule_name, element_names
# read_user_mapping(user_mapping_filename='propane_user_mapping.xml')

def read_system_info(struct_filename):
    # parser = ET.XMLParser(encoding="utf-8")
    # root = ET.fromstring(struct_filename, parser=parser)
    root = ET.fromstring(open(struct_filename).read())
    # root = ET.parse(struct_filename).getroot()
    n_unitsTotal = int(root.find('configuration').find('position').attrib['num'])
    print(n_unitsTotal)

    return n_unitsTotal

# read_system_info(struct_filename='start_aa.hoomdxml')

def create_system_mapping(element_names, n_sections_TOTAL, t):
    # Initialize atoms with elements
    ## for loop to traverse element_names array for elements
    ## first test - only use for carbon (one element)
    ## then expand later
    ## need to allow for different types of elements together
    ## maybe use element library from msibi repo
    from mdtraj.core import element
    list(t.top.atoms)[0].element = element.carbon # check element
    list(t.top.atoms)[0].element.mass
    for atom in t.top.atoms:
        atom.element = element.carbon # check element

    # Map the beads accordingly
    cg_idx = 0
    start_idx = 0
    propane_map = {0: [0, 1, 2]} ## mapping definition needs to be created
                                    # from search and user files
    ## TEST CODE
    ######################################################################
    ######################################################################



    ######################################################################
    ######################################################################
    ## TEST CODE


    system_mapping = {}
    for n in range(n_sections_TOTAL):
        for bead, atoms in propane_map.items():
            system_mapping[cg_idx] = [x + start_idx for x in atoms]
            start_idx += len(atoms) # understand this part
            cg_idx += 1

    # Apply mapping for XYZ coordinates
    cg_xyz = np.empty((t.n_frames, len(system_mapping), 3))
    for cg_bead, aa_indices in system_mapping.items():
        cg_xyz[:, cg_bead, :] = md.compute_center_of_mass(t.atom_slice(aa_indices))

    # Apply mapping for Topology object
    cg_top = md.Topology()
    for cg_bead in system_mapping.keys(): #i got the keys keys keys
        cg_top.add_atom('carbon', element.virtual_site, cg_top.add_residue('A',
                            cg_top.add_chain()))
        ## Check element and name for items 'A'
        ## Possible interface with mbuild for better UI and aesthetics


    return cg_xyz, cg_top

def compute_files(cg_xyz, cg_top, t, molecule_name):
    # Create Trajectory object and save to .dcd file
    cg_traj = md.Trajectory(cg_xyz, cg_top, time=None,
                            unitcell_lengths=t.unitcell_lengths,
                            unitcell_angles=t.unitcell_angles)
    ## need better file naming convention and rule guideline
    cg_traj.save_dcd('cg_traj_{0:s}.dcd'.format(molecule_name)) ## need check statements to prevent file overwrite
                                    ## rename old/new files accordingly

    # Create rdfs file from pairs
    ## might need for loop if more elements (later implementation)
    ## need some way to recognize pairs of units - not just carbon elements if expanded
    pairs = cg_traj.top.select_pairs(selection1='name {0:s}'.format(element_names[0]), ## Check element
                                     selection2='name {0:s}'.format(element_names[0])) ## Check element
    r, g_r = md.compute_rdf(cg_traj, pairs=pairs,
                            r_range=(0, 1.2), bin_width=0.005)
                    ## identify end of range with data pairs
                    ## maybe something with read-file function - compare data values next to each other
                    ## See where data drop-off occurs and plot respectively
                    ## maybe use slope - negative less than some number, set as cutoff
                    ## record cutoff point somewhere for debugging purposes
    np.savetxt('rdfs_aa.txt', np.transpose([r, g_r])) # need check statements to prevent file overwrite

    plot_output(r, g_r, molecule_name)


# In[8]:

def plot_output(x, y, molecule_name):
    ## modify figsize according to drop-off point
    plt.figure(num=None, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(x, y, label="CG {0:s}".format(molecule_name)) ## Check label use string format method

    ## look up more aesthetic matplotlib functions
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.legend()

def convert_Traj_RDF():
    ## add parameters to function calls - maybe add other functions
    ## other potential functions:
    ##     - check_file_overwrite()
    ##     - check_g_r_dropoff() - integrate with plot function
    ##     - manage_filetypes() - read in files

    traj_filename = 'traj_unwrapped.dcd'
    struct_filename = 'start_aa.hoomdxml'
    search_mapping_filename = 'propane_search_mapping.xml'
    user_mapping_filename = 'propane_user_mapping.xml'

    t = md.load(traj_filename, top=struct_filename)
    # t.save('propane.mol2')
    # struct_parmed = pmd.load_file('propane.mol2')

    n_units_TOTAL = read_system_info(struct_filename)
    n_unitsPerSection, molecule_name, element_names = read_user_mapping(user_mapping_filename)
    n_sections_TOTAL = n_units_TOTAL // n_unitsPerSection

    cg_xyz, cg_top = create_system_mapping(element_names, n_sections_TOTAL, t)
    compute_files(cg_xyz, cg_top, t, molecule_name)

# Execute functions
## maybe use dictionary for element_names?
## incorporate the arguments in initialization for map_beads()
## maybe initialize element_names array from read-in file with bonds
## bonds recorded in structure indicate elements involved in rdf
convert_Traj_RDF()
