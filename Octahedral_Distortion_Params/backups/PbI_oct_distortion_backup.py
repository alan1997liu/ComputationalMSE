#!/usr/bin/env python3
import heapq
import pymatgen as mg
import os, sys, getopt
from pymatgen.io import vasp
import numpy as np
from numpy import linalg as linalg
import warnings

warnings.simplefilter("error", RuntimeWarning)
#a = species_arr[0].symbol
#print(type(a))
# We want to get the indexes of all the atoms in the species_array of the
# same atom type, basically. 

# The function below gets all the Pb_indexes and the I_indexes so that
# we may do something with their indices to calculate relevant info.
def getSpeciesInfo(species_arr):
    Pb_indexes = []   
    I_indexes = []
    Pb_nearest_neighbors = []
    Pb_coords = []
    I_coords = []
    for atom_index in range(0, len(species_arr)):
        if "Pb" in species_arr[atom_index].symbol:
            Pb_indexes.append(atom_index)
            Pb_coords.append(struct[atom_index].coords)
        if "I" in species_arr[atom_index].symbol:
            I_indexes.append(atom_index)
            I_coords.append(struct[atom_index].coords)
    return Pb_indexes, I_indexes, Pb_coords, I_coords
#---------------------------------------------------------------------------
# Coordinates are stored in an array, in which each element is also an
# array that contains the coordinates information.
# Tried to calculate distance between two atoms that were not in the same
# unit cell, but still got closest distance. Therefore, I don't see
# a problem, but I may be wrong? Below are **tests** of Pb 
## print(struct.get_distance(6, 9))
## print(struct.get_distance(1, 4))
## print(struct.get_distance(2, 4))
# This function gets an octahedral array. The octahedral array is an array
# containing arrays that have:
# 1. The center Pb-index as the first element. 
# 2. The 6 I-indexes as the next six elements
def get_octahedrals(Pb_indexes, I_indexes, coord_num):
    octahedron_array = []
    for Pb_index in Pb_indexes:
        I_indexes_distances = {}
        for I_index in I_indexes:
            I_indexes_distances[I_index] = struct.get_distance(Pb_index, I_index)

        octahedral_indices = [Pb_index]
        # There are 6 Iodines that are the closest to one Pb
        closest_I_distances = heapq.nsmallest(coord_num, \
                I_indexes_distances.values())
        for value in closest_I_distances:
            closest_I_index = [k for k in I_indexes_distances.keys() \
                    if I_indexes_distances[k] == value and k not in \
                    octahedral_indices]
            octahedral_indices.extend(closest_I_index)
        octahedron_array.append(octahedral_indices)
    return octahedron_array

def calc_DistortionAngles(Pb_coords, shared_I_coord, center_Pb, struct, \
        index, Pb_indexes, shared_I_indexes, centerPbIndex):
    # To calculate the in-plane distortion angle, we must first calculate
    # the projection of the Pb-I vector onto the Pb-plane. Then, we must
    # calculate the angle from the dot product of Pb-Pb and the projected
    # Pb-I vector onto the plane. 
    def PbPlane(Pb_coords):
        Pb_v1 = np.subtract(Pb_coords[index], center_Pb)
        Pb_v2 = np.subtract(Pb_coords[(index + 1)%len(Pb_coords)], center_Pb)
        cross_PbVectors = np.cross(Pb_v1, Pb_v2)
        nVect_PbPlane = cross_PbVectors / linalg.norm(cross_PbVectors)
        return nVect_PbPlane

    def inPlaneAngle(Pb_coords, I_coord):
        nVector_PbPlane = PbPlane(Pb_coords)
        PbI_v = np.subtract(I_coord, center_Pb)
        # Need to project PbI-vector onto the normal to the PbPlane first
        proj_PbI_norm = nVector_PbPlane * np.dot(PbI_v, nVector_PbPlane)
        # Next, need to subtract the above projection from Pb-I vector.
        # This lets us obtain the part of the Pb-I vector that's on the
        # Pb-plane, basically.
        proj_PbI_PbPlane = np.subtract(PbI_v, proj_PbI_norm)
        # Then, we have to get the position of the projection of the
        # iodine atom onto the PbPlane!
        proj_I_PbPlane = np.add(proj_PbI_PbPlane, center_Pb)
        # Finally, let's get the two vectors pointing away from
        # proj_I_PbPlane. Then, find the inPlaneAngle from the dot prod.
        v1 = np.subtract(center_Pb, proj_I_PbPlane)
        v2 = np.subtract(Pb_coords[index], proj_I_PbPlane)
        inPlaneAngle = np.arccos(np.dot(v1, v2) / \
                (linalg.norm(v1) * linalg.norm(v2))) * 180/np.pi
        return inPlaneAngle

    def outPlaneAngle(Pb_coords, I_coord):
        # Get projection of Iodine atom onto the PbPlane norm
        nVector_PbPlane = PbPlane(Pb_coords)
        PbI_v = np.subtract(I_coord, center_Pb)
        proj_PbI_norm = nVector_PbPlane*np.dot(PbI_v, nVector_PbPlane)
        # Now get projection of PbI vector (PbI_v) onto Pb-Pb vector
        Pb_vector = np.subtract(Pb_coords[index], center_Pb)
        proj_PbI_PbPbvector = (np.dot(Pb_vector, PbI_v) / \
                (linalg.norm(Pb_vector)**2)) * Pb_vector
        # Now get position of Iodine that is in the plane that is
        # perpendicular to the Pb-Pb plane
        proj_I_outofPlane = np.add(center_Pb, proj_PbI_PbPbvector, \
                proj_PbI_norm)
        # Now get the angles!
        v1 = np.subtract(center_Pb, proj_I_outofPlane)
        v2 = np.subtract(Pb_coords[index], proj_I_outofPlane)
        dot_product = (np.dot(v1, v2) / (linalg.norm(v1)*linalg.norm(v2)))
        try:
            outPlaneAngle = np.arccos(dot_product) * 180 / np.pi
        except RuntimeWarning:
            print("Cannot compute arccos for outPlaneAngle.", \
                    "The dot product of interest is", dot_product)
            outPlaneAngle = 0
            pass
        return outPlaneAngle
    
    # Calculates bond angles that matches up to VESTA bond angles
    # The inputs should be the indices that correspond to each atom.
    def tiltAngle(center_Pb_index, other_Pb_index, shared_I_index):
        return struct.get_angle(center_Pb_index, shared_I_index, other_Pb_index)

    #Main method for calc_distortion_angles goes here 
    inPlaneAngle = inPlaneAngle(Pb_coords, shared_I_coord)
    outPlaneAngle = outPlaneAngle(Pb_coords, shared_I_coord)
    tiltAngle = tiltAngle(centerPbIndex , Pb_indexes[index], shared_I_indexes[index])
    inPlaneDistortion = 180 - inPlaneAngle
    outPlaneDistortion = 180 - outPlaneAngle
    tiltDistortion = 180 - tiltAngle
    return inPlaneDistortion, outPlaneDistortion, tiltDistortion


# -----------------------------------------------------------------------
# Given a 3 x 1 x 3 supercell, return the octahedral arrays that are in
# the middle of the supercell.
# Hypothesis: The fifth cell is the supercell. 
def get_center_octahedrals(octahedral_array):
    num_unitcells = len(octahedral_array)/4
    num_unitcells = int(num_unitcells)
    center_octahedrals = []
    center_Pb = 4
    for i in range(0, 4):
        center_octahedrals.append(octahedral_array[center_Pb])
        center_Pb += num_unitcells
    return center_octahedrals

# Return an array of Pb_indexes and the shared I_indexes
# select_octahedral input is the octahedral of interest.
# Want to find all of its shared iodines (a total of 4) also all the
# Pb_indexes that it shares iodines with
def get_shared_I(select_octahedral, octahedral_array):
    Pb_indexes = []
    shared_I_indexes = []
    for i in range(0, len(octahedral_array)):
        for j in range(1, len(octahedral_array[i])):
            if octahedral_array[i][j] in select_octahedral:
                if (octahedral_array[i][0] == select_octahedral[0]):
                    continue
                else:
                    shared_I_indexes.append(octahedral_array[i][j])
                    Pb_indexes.append(octahedral_array[i][0])
    return Pb_indexes, shared_I_indexes

#------------------------------------------------------------------------

# The function below is the main function call. Call this function to get
# all the relevant distortion parameters!
# This function gets the coordinates of all the Pb atoms that share an
# iodine with the center Pb atom of choice in center_octahedrals
def get_distortion_info(center_octahedrals, octahedral_array, struct):
    for i in range(0, len(center_octahedrals)):
        print("Information for octahedral", i + 1)
        Pb_indexes, shared_I_indexes = get_shared_I(center_octahedrals[i], octahedral_array)
        print(shared_I_indexes)
        print(center_octahedrals[i])
        center_PbIndex = center_octahedrals[i][0]
        center_Pb = struct[center_PbIndex].coords
        shared_I_coords = []
        Pb_coords = []
        for k in range(0, len(Pb_indexes)):
            Pb_coords.append(struct[Pb_indexes[k]].coords)
            shared_I_coords.append(struct[shared_I_indexes[k]].coords)
        for j in range(0, len(Pb_coords)):
            inPlaneDistortion, outPlaneDistortion, tiltDistortion =\
            calc_DistortionAngles(Pb_coords, shared_I_coords[j], center_Pb,\
            struct, j, Pb_indexes, shared_I_indexes, center_PbIndex)
            print("Bond", j + 1)
            print("In plane distortion:", inPlaneDistortion)
            print("Out of plane distortion:", outPlaneDistortion)
            print("Tilting angle:", tiltDistortion)
            print()
        
        bond_distortion = bond_length_distortion(center_octahedrals[i], struct)
        print("Bond length distortion:", bond_distortion)
        print()
        print()
        #octahedral_elongation = octahedral_elongation(center_octahedrals[i], \
                #sharedIIndexes, struct)

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

def octahedral_elongation(center_octahedral_array, sharedIIndexes, struct):
    return
#---------------------------------------------------------------------------

struct = vasp.inputs.Poscar.from_file("PBE_CONTCAR").structure
Pb_coordination = 6
#species_arr = struct.species
#Pb_indexes, I_indexes, Pb_coords, I_coords = getSpeciesInfo(species_arr)
#octahedron_array = get_octahedrals(Pb_indexes, I_indexes, Pb_coordination)
#print(octahedron_array)
#tilting_distortion, inPlane_distortion, outPlane_distortion = \
        #calc_DistortionAngles(octahedron_array, struct)


# Supercell the structure. 
# In order to calculate the distortions of all four bonds. You have to
# supercell the structure and impose periodic boundary conditions.
struct.make_supercell([3, 1, 3])
supercell_lattice_vectors = struct.lattice
supercell_species_arr = struct.species
Pb_indexes, I_indexes, Pb_coords, I_coords = getSpeciesInfo(supercell_species_arr)
octahedral_array = get_octahedrals(Pb_indexes, I_indexes, Pb_coordination)
center_octahedrals = get_center_octahedrals(octahedral_array)
# Function call that outputs things to the terminal
# The below function will print out all the octahedral distortion
# parameters of interest in an organized fashion on the terminal.
get_distortion_info(center_octahedrals, octahedral_array, struct)
print("End of distortion information")
print()
