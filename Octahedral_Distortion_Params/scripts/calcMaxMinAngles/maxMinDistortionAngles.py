#!/usr/bin/env python3
import heapq
import pymatgen as mg
import os, sys, getopt
from pymatgen.io import vasp
import numpy as np
from numpy import linalg as linalg
import warnings
import csv
import sys

warnings.simplefilter("error", RuntimeWarning)
#a = species_arr[0].symbol
#print(type(a))
# We want to get the indexes of all the atoms in the species_array of the
# same atom type, basically. 

# The function below gets all the B_indexes and the X_indexes so that
# we may do something with their indices to calculate relevant info.
def getSpeciesInfo(species_arr, B_atom, X_atom):
    B_indexes = []   
    X_indexes = []
    B_nearest_neighbors = []
    B_coords = []
    X_coords = []
    for atom_index in range(0, len(species_arr)):
        if B_atom in species_arr[atom_index].symbol:
            B_indexes.append(atom_index)
            B_coords.append(struct[atom_index].coords)
        if X_atom in species_arr[atom_index].symbol:
            X_indexes.append(atom_index)
            X_coords.append(struct[atom_index].coords)
    return B_indexes, X_indexes, B_coords, X_coords
#---------------------------------------------------------------------------
# Coordinates are stored in an array, in which each element is also an
# array that contains the coordinates information.
# Tried to calculate distance between two atoms that were not in the same
# unit cell, but still got closest distance. Therefore, no problem, probably.
# This function gets an octahedral array. The octahedral array is an array
# containing arrays that have:
# 1. The center B-index as the first element. 
# 2. The 6 X-indexes as the next six elements
def get_octahedrals(B_indexes, X_indexes, coord_num):
    octahedron_array = []
    for B_index in B_indexes:
        X_indexes_distances = {}
        for X_index in X_indexes:
            X_indexes_distances[X_index] = struct.get_distance(B_index, X_index)

        octahedral_indices = [B_index]
        # There are 6 Xodines that are the closest to one B
        closest_X_distances = heapq.nsmallest(coord_num, \
                X_indexes_distances.values())
        for value in closest_X_distances:
            closest_X_index = [k for k in X_indexes_distances.keys() \
                    if X_indexes_distances[k] == value and k not in \
                    octahedral_indices]
            octahedral_indices.extend(closest_X_index)
        octahedron_array.append(octahedral_indices)
    return octahedron_array

def calc_DistortionAngles(B_coords, shared_X_coord, center_B, struct, \
        index, B_indexes, shared_X_indexes, centerBIndex):
    # To calculate the in-plane distortion angle, we must first calculate
    # the projection of the B-X vector onto the B-plane. Then, we must
    # calculate the angle from the dot product of B-B and the projected
    # B-X vector onto the plane. 
    def BPlane(B_coords):
        B_v1 = np.subtract(B_coords[index], B_coords[(index + 1)%len(B_coords)])
        B_v2 = np.subtract(B_coords[index], B_coords[(index + 2)%len(B_coords)])
        cross_BVectors = np.cross(B_v1, B_v2)
        nVect_BPlane = cross_BVectors / linalg.norm(cross_BVectors)
        return nVect_BPlane

    def inPlaneAngle(B_coords, X_coord):
        nVector_BPlane = BPlane(B_coords)
        BX_v = np.subtract(X_coord, center_B)
        # Need to project BX-vector onto the normal to the BPlane first
        proj_BX_norm = nVector_BPlane * np.dot(BX_v, nVector_BPlane)
        # Next, need to subtract the above projection from B-X vector.
        # This lets us obtain the part of the B-X vector that's on the
        # B-plane, basically.
        proj_BX_BPlane = np.subtract(BX_v, proj_BX_norm)
        # Then, we have to get the position of the projection of the
        # iodine atom onto the BPlane!
        proj_X_BPlane = np.add(proj_BX_BPlane, center_B)
        # Finally, let's get the two vectors pointing away from
        # proj_X_BPlane. Then, find the inPlaneAngle from the dot prod.
        v1 = np.subtract(center_B, proj_X_BPlane)
        v2 = np.subtract(B_coords[index], proj_X_BPlane)
        inPlaneAngle = np.arccos(np.dot(v1, v2) / \
                (linalg.norm(v1) * linalg.norm(v2))) * 180/np.pi
        return inPlaneAngle

    def outPlaneAngle(B_coords, X_coord):
        # Get vector projection of X atom onto the BPlane norm to get vertical
        # This is the second method that I'm trying.
        nVector_BPlane = BPlane(B_coords)
        BX_v = np.subtract(X_coord, center_B)
        proj_BX_norm = nVector_BPlane * np.dot(BX_v, nVector_BPlane)

        B_vector = np.subtract(B_coords[index], center_B)
        cross_Bvect_Nvect = np.cross(B_vector, nVector_BPlane)
        unit_cross_Bvect_Nvect = cross_Bvect_Nvect / linalg.norm(cross_Bvect_Nvect)
        proj_BX_unit_cross_Bvect_Nvect = unit_cross_Bvect_Nvect \
                * np.dot(BX_v, unit_cross_Bvect_Nvect)
        vector_on_normal_plane = BX_v - proj_BX_unit_cross_Bvect_Nvect
        X_atom_on_normal = np.add(vector_on_normal_plane, center_B)
        v1 = np.subtract(center_B , X_atom_on_normal)
        v2 = np.subtract(B_coords[index], X_atom_on_normal)

        dot_product = (np.dot(v1, v2) / (linalg.norm(v1)*linalg.norm(v2)))
        try:
            outPlaneAngle = np.arccos(dot_product) * 180 / np.pi
        except RuntimeWarning:
            outPlaneAngle = 180
        return outPlaneAngle
    
    # Calculates bond angles that matches up to VESTA bond angles
    # The inputs should be the indices that correspond to each atom.
    def tiltAngle(center_B_index, other_B_index, shared_X_index):
        return struct.get_angle(center_B_index, shared_X_index, other_B_index)

    #Main method for calc_distortion_angles goes here 
    inPlaneAngle = inPlaneAngle(B_coords, shared_X_coord)
    outPlaneAngle = outPlaneAngle(B_coords, shared_X_coord)
    tiltAngle = tiltAngle(centerBIndex , B_indexes[index], shared_X_indexes[index])
    inPlaneDistortion = 180 - inPlaneAngle
    outPlaneDistortion = 180 - outPlaneAngle
    tiltDistortion = 180 - tiltAngle
    return inPlaneDistortion, outPlaneDistortion, tiltDistortion


# -----------------------------------------------------------------------
# Given a 3 x 1 x 3 supercell, return the octahedral arrays that are in
# the middle of the supercell.
# Hypothesis: The fifth cell is the supercell. 
def get_center_octahedrals(octahedral_array, numUnitCells):
    num_unitcells = len(octahedral_array)/4
    num_unitcells = int(num_unitcells)
    center_octahedrals = []
    center_B = int(numUnitCells/2)
    for i in range(0, 4):
        center_octahedrals.append(octahedral_array[center_B])
        center_B += num_unitcells
    return center_octahedrals

# Return an array of B_indexes and the shared X_indexes
# select_octahedral input is the octahedral of interest.
# Want to find all of its shared iodines (a total of 4) also all the
# B_indexes that it shares iodines with. 
# If the B_index is the same as the center, we ignore it and continue with the loop
def get_shared_X(select_octahedral, octahedral_array):
    B_indexes = []
    shared_X_indexes = []
    for i in range(0, len(octahedral_array)):
        for j in range(1, len(octahedral_array[i])):
            if octahedral_array[i][j] in select_octahedral:
                if (octahedral_array[i][0] == select_octahedral[0]):
                    continue
                else:
                    shared_X_indexes.append(octahedral_array[i][j])
                    B_indexes.append(octahedral_array[i][0])
    return B_indexes, shared_X_indexes

#------------------------------------------------------------------------

# The function below is the main function call. Call this function to get
# all the relevant distortion parameters!
# This function gets the coordinates of all the B atoms that share an
# iodine with the center B atom of choice in center_octahedrals
def get_distortion_info(center_octahedrals, octahedral_array, struct):
    uniqueInPlane = []
    uniqueOutPlane = []
    uniqueTilt = []
    for i in range(0, len(center_octahedrals)):
        inPlaneArr = []
        outPlaneArr = []
        tiltingArr = []

        #print("Information for octahedral", i + 1)
        B_indexes, shared_X_indexes = get_shared_X(center_octahedrals[i], octahedral_array)
        center_BIndex = center_octahedrals[i][0]
        center_B = struct[center_BIndex].coords
        shared_X_coords = []
        B_coords = []
        for k in range(0, len(B_indexes)):
            B_coords.append(struct[B_indexes[k]].coords)
            shared_X_coords.append(struct[shared_X_indexes[k]].coords)
        for j in range(0, len(B_coords)):
            inPlaneDistortion, outPlaneDistortion, tiltDistortion =\
            calc_DistortionAngles(B_coords, shared_X_coords[j], center_B,\
            struct, j, B_indexes, shared_X_indexes, center_BIndex)
            if inPlaneDistortion not in uniqueInPlane:
                uniqueInPlane.append(inPlaneDistortion)
            if outPlaneDistortion not in uniqueOutPlane:
                uniqueOutPlane.append(outPlaneDistortion)
            if tiltDistortion not in uniqueTilt:
                uniqueTilt.append(tiltDistortion)
        bond_distortion = bond_length_distortion(center_octahedrals[i], struct)
        

    with open('data.csv', 'a', newline = '') as csv_file:
        data_writer = csv.writer(csv_file, delimiter = ",", quotechar = '"', \
                    quoting = csv.QUOTE_MINIMAL)
        data_writer.writerow([filename, max(uniqueInPlane), min(uniqueInPlane), \
                max(uniqueOutPlane), min(uniqueOutPlane), max(uniqueTilt), min(uniqueTilt)])
        
        #*** HAVE TO FINISH UP THIS METHOD!
        #octahedral_elongation = octahedral_elongation(center_octahedrals[i], \
                #sharedXXndexes, struct)

def bond_length_distortion(octahedral_array, struct):
    d_m = 0
    bond_length_distortion = 0
    for i in range(1, len(octahedral_array)):
        d_m += struct.get_distance(octahedral_array[0], octahedral_array[i]) \
                / (len(octahedral_array) - 1)

    for i in range(1, len(octahedral_array)):
        d_i = struct.get_distance(octahedral_array[0], octahedral_array[i])
        bond_length_distortion += (((d_i - d_m) / d_m)**2)/6
    return bond_length_distortion

# Have to finish defining this method. Cannot average bond lengths to calculate
# volume. http://mathcentral.uregina.ca/QQ/database/QQ.09.06/haivan1.html
# An irregular octahedron can be dissected into four irregular tetrahedra. 
# The volume of any tetrahedron is a third the base area times the perpendicular
# height. Or use formula from this website:
# https://iopscience.iop.org/article/10.1070/RM1995v050n05ABEH002620/pdf
def octahedral_elongation(center_octahedral_array, sharedXXndexes, struct):
    # np.matrix([[5, 6, 7], [4, 6]]) --> ex. how to create numpy array
    return

def oct_angle_variance(center_octahedrals, struct):
    # Take the first octahedral only for testing
    all_angles = []
    PbIndex = center_octahedrals[0][0]
    for i in range(1, len(center_octahedrals[0]) - 1):
        for j in range(i + 1, len(center_octahedrals[0])):
            all_angles.append(struct.get_angle(center_octahedrals[0][i], PbIndex, \
                    center_octahedrals[0][j]))
    oct_angles = heapq.nsmallest(12, all_angles)
    oct_angles_sum = 0
    for i in range(0, oct_angles):
        oct_angles_sum = (oct_angles[i] - 90)**2
    return oct_angles_sum / 11

    #for i in range(0, len(center_octahedrals)):
        #oct_angles = []
        #PbIndex = center_octahedrals[i][0]
        #for j in range(0, len(center_octahedrals[i])):
            #for k in range(j, len(center_octahedrals[i])):
                #oct_angles.append(struct.get_angle(center_octahedrals[i][j], \
                        #PbIndex, center_octahedrals[i][k]))


#---------------------------------------------------------------------------

# Set the command line arguments to read in B atom and X atom.
filename = sys.argv[1]
B_atom = sys.argv[2]
X_atom = sys.argv[3]

struct = vasp.inputs.Poscar.from_file(filename).structure
B_coordination = 6
# print("Information for:", filename)

# Supercell the structure. 
# Xn order to calculate the distortions of all four bonds. You have to
# supercell the structure and impose periodic boundary conditions.
# You can choose how to supercell, just vary the numbers a, b, and c.
a = 3
b = 3
c = 3
struct.make_supercell([a, b, c])
numUnitCells = a*b*c
supercell_lattice_vectors = struct.lattice
supercell_species_arr = struct.species
B_indexes, X_indexes, B_coords, X_coords = getSpeciesInfo(supercell_species_arr, \
        B_atom, X_atom)
octahedral_array = get_octahedrals(B_indexes, X_indexes, B_coordination)
center_octahedrals = get_center_octahedrals(octahedral_array, numUnitCells)

# Function call that outputs things to the terminal
# The below function will print out all the octahedral distortion
# parameters of interest in an organized fashion on the terminal.

# get_distortion_info(center_octahedrals, octahedral_array, struct)
print(oct_angle_variance(center_octahedrals))

