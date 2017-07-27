# # testAPI_automate
# Created by Davy Yue 2017-06-26

import itertools
import string
import os
from xml.etree import cElementTree as ET

from foyer.smarts_graph import SMARTSGraph
from foyer.smarts import SMARTS as SMARTSParser
import parmed as pmd
from periodic import element
import numpy as np
import matplotlib.pyplot as plt
# get_ipython().magic('matplotlib inline') # unknown error need to check this
import time

from msibi import MSIBI, State, Pair, mie
import mdtraj as md


def read_search_mapping(search_mapping_filename, user_mapping_filename, structure):
    # parser = ET.XMLParser(encoding="utf-8")
    # root = ET.fromstring(search_mapping_filename, parser=parser)
    root = ET.fromstring(open(search_mapping_filename).read())
    # root = ET.parse(search_mapping_filename).getroot()

    searchlist = [] # list containing all search values ordered by priority
    for value in root.findall('value'):
        searchlist.append(value.attrib['searchstr'])
        print(searchlist)

    parser = SMARTSParser()
    matches = []

    for searchstr in searchlist:
        graph = SMARTSGRAPH(searchstr, parser=parser)
        matches.append(graph.find_matches(structure))

    return matches

# SMARTS string
# Current supported SMARTS string: https://github.com/mosdef-hub/foyer/issues/63

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
    element_names.append('carbon')
    n_sections_BEAD = 0
    n_unitsPerSection = 0 # number units per section being combined
    for section in bead.findall('section'):
        # allow for different number of units in each section per bead
        # Element identification using Periodic python package
        # https://github.com/luisnaranjo733/periodic
        # use foyer here to find indices

        n_unitsPerSection = section.attrib['elements'].count("C")
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
    # SLOWEST PART OF CODE IS THIS FUNCTION
    # Initialize atoms with elements
    ## for loop to traverse element_names array for elements
    ## need to expand from just carbon to more/different elements
    ## maybe use elements from periodic package
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
    for n in range(n_sections_TOTAL): # what does sections mean in this particular context
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

def compute_files(cg_xyz, cg_top, t, molecule_name, element_names):
    # Create Trajectory object and save to .dcd file
    cg_traj = md.Trajectory(cg_xyz, cg_top, time=None,
                            unitcell_lengths=t.unitcell_lengths,
                            unitcell_angles=t.unitcell_angles)
    ## need better file naming convention and rule guideline
    filepath = os.path.join(os.getcwd(), 'data/cg_traj_{0:s}.dcd'.format(molecule_name))
    cg_traj.save_dcd(filepath) ## need check statements to prevent file overwrite
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
    filepath = os.path.join(os.getcwd(), 'data/rdfs_aa.txt')
    np.savetxt(filepath, np.transpose([r, g_r])) # need check statements to prevent file overwrite
    print("Saved rdfs to file")
    plot_output(r, g_r, molecule_name)

def plot_output(x, y, molecule_name):
    ## modify figsize according to drop-off point
    plt.figure(num=None, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(x, y, label="CG {0:s}".format(molecule_name)) ## Check label use string format method

    ## look up more aesthetic matplotlib functions
    plt.title('Propane CG Graph at {0:s}'.format(time.strftime("%m/%d/%Y %I:%M:%S")))
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.legend()
    filepath = os.path.join(os.getcwd(), 'data/trajAPI_plot_{0:s}.pdf'.format(molecule_name))
    plt.savefig(filepath)
    print("Figure should be saved to data folder")

def convert_Traj_RDF():
    ## add parameters to function calls - maybe add other functions
    ## other potential functions:
    ##     - check_file_overwrite()
    ##     - check_g_r_dropoff() - integrate with plot function
    ##     - manage_filetypes() - read in files (maybe for mdtraj flexibility)

    traj_filename = os.path.join(os.getcwd(), 'data/traj_unwrapped.dcd')
    struct_filename = os.path.join(os.getcwd(), 'data/start_aa.hoomdxml')
    search_mapping_filename = os.path.join(os.getcwd(),
                                'data/propane_search_mapping.xml')
    user_mapping_filename = os.path.join(os.getcwd(),
                                'data/propane_user_mapping.xml')

    t = md.load(traj_filename, top=struct_filename)
    print("Loaded struct & traj files")

    top = t.top.to_openmm(traj=t) #openmm topology

    print(read_search_mapping(search_mapping_filename, user_mapping_filename, structure))

    # n_units_TOTAL = read_system_info(struct_filename)
    # print("Read in system info from struct file")
    # n_unitsPerSection, molecule_name, element_names = read_user_mapping(user_mapping_filename)
    # print("Read in user_mapping file")
    # n_sections_TOTAL = n_units_TOTAL // n_unitsPerSection

    # cg_xyz, cg_top = create_system_mapping(element_names, n_sections_TOTAL, t)
    # print("Created system mapping")
    # compute_files(cg_xyz, cg_top, t, molecule_name, element_names)

# Execute functions
## maybe initialize element_names array from read-in file with bonds
## bonds recorded in structure indicate elements involved in rdf
convert_Traj_RDF()
