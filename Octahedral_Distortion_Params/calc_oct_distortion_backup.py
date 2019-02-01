#!/usr/bin/env python3
import heapq
import pymatgen as mg
import os, sys, getopt
from pymatgen.io import vasp
import numpy as np
from numpy import linalg as linalg

#a = species_arr[0].symbol
#print(type(a))

# We want to get the indexes of all the atoms in the species_array of the
# same atom type, basically. 

def getSpeciesInfo(species_array):
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
    octahedral_array = []
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
        octahedral_array.append(octahedral_indices)
    return octahedral_array

def calc_DistortionAngles(octahedral_array, struct):
    # To calculate the in-plane distortion angle, we must first calculate
    # the projection of the Pb-I vector onto the Pb-plane. Then, we must
    # calculate the angle from the dot product of Pb-Pb and the projected
    # Pb-I vector onto the plane. 
    def PbPlane(Pb_coords):
        Pb_v1 = np.subtract(Pb_coords[3], Pb_coords[0])
        Pb_v2 = np.subtract(Pb_coords[2], Pb_coords[0])
        cross_PbVectors = np.cross(Pb_v1, Pb_v2)
        nVect_PbPlane = cross_PbVectors / linalg.norm(cross_PbVectors)
        return nVect_PbPlane

    def inPlaneAngle(Pb_coords, I_coord):
        nVector_PbPlane = PbPlane(Pb_coords)
        PbI_v = np.subtract(I_coord, Pb_coords[0])
        # Need to project PbI-vector onto the normal to the PbPlane first
        proj_PbI_norm = nVector_PbPlane * np.dot(PbI_v, nVector_PbPlane)
        # Next, need to subtract the above projection from Pb-I vector.
        # This lets us obtain the part of the Pb-I vector that's on the
        # Pb-plane, basically.
        proj_PbI_PbPlane = np.subtract(PbI_v, proj_PbI_norm)
        # Then, we have to get the position of the projection of the
        # iodine atom onto the PbPlane!
        proj_I_PbPlane = np.add(proj_PbI_PbPlane, Pb_coords[0])
        # Finally, let's get the two vectors pointing away from
        # proj_I_PbPlane. Then, find the inPlaneAngle from the dot prod.
        v1 = np.subtract(Pb_coords[0], proj_I_PbPlane)
        v2 = np.subtract(Pb_coords[3], proj_I_PbPlane)
        inPlaneAngle = np.arccos(np.dot(v1, v2) / \
                (linalg.norm(v1) * linalg.norm(v2))) * 180/np.pi
        return inPlaneAngle

    def outPlaneAngle(Pb_coords, I_coord):
        # Get projection of Iodine atom onto the PbPlane norm
        nVector_PbPlane = PbPlane(Pb_coords)
        PbI_v = np.subtract(I_coord, Pb_coords[0])
        proj_PbI_norm = nVector_PbPlane*np.dot(PbI_v, nVector_PbPlane)
        # Now get projection of PbI vector (PbI_v) onto Pb-Pb vector
        Pb_vector = np.subtract(Pb_coords[3], Pb_coords[0])
        proj_PbI_PbPbvector = (np.dot(Pb_vector, PbI_v) / \
                (linalg.norm(Pb_vector)**2)) * Pb_vector
        # Now get position of Iodine that is in the plane that is
        # perpendicular to the Pb-Pb plane
        proj_I_outofPlane = np.add(Pb_coords[0], proj_PbI_PbPbvector, \
                proj_PbI_norm)
        # Now get the angles!
        v1 = np.subtract(Pb_coords[0], proj_I_outofPlane)
        v2 = np.subtract(Pb_coords[3], proj_I_outofPlane)
        outPlaneAngle = np.arccos(np.dot(v1, v2) / (linalg.norm(v1) * \
                linalg.norm(v2))) * 180/np.pi
        return outPlaneAngle
    
    # Calculates bond angles that matches up to VESTA bond angles
    def tiltingAngle_and_correctIodine(octahedral_array):
        #Get all the shared iodines (there should be 4)
        shared_iodines = []
        for j in range(1, len(octahedral_array[0])):
            if octahedral_array[0][j] in octahedral_array[3] and \
                    octahedral_array[0][j] not in shared_iodines:
                shared_iodines.append(octahedral_array[0][j])
        #Get the correct iodine that is within the unit cell!
        for shared_iodine in shared_iodines:
            bond_angle = struct.get_angle(octahedral_array[0][0], \
                shared_iodine, octahedral_array[3][0]) 
            if bond_angle > 100:
                correct_bond_angle = bond_angle
                shared_iodine_index = shared_iodine
                break
        return correct_bond_angle, shared_iodine_index

    Pb_coords = []
    for i in range(0, 4):
        Pb_coords.append(struct[octahedral_array[i][0]].coords)

    tilting_Angle, shared_iodine_index = \
            tiltingAngle_and_correctIodine(octahedral_array)
    I_coord = struct[shared_iodine_index].coords
    print(tilting_Angle)
    tilting_distortion = 180 - tilting_Angle
    inPlane_distortion = 180 - inPlaneAngle(Pb_coords, I_coord)
    outPlane_distortion = 180 - outPlaneAngle(Pb_coords, I_coord)
    print(tilting_distortion, inPlane_distortion, outPlane_distortion)
    return tilting_distortion, inPlane_distortion, outPlane_distortion

#---------------------------------------------------------------------------

struct = vasp.inputs.Poscar.from_file("PBE_CONTCAR").structure
species_arr = struct.species
Pb_indexes, I_indexes, Pb_coords, I_coords = getSpeciesInfo(species_arr)
Pb_coordination = 6
octahedral_array = get_octahedrals(Pb_indexes, I_indexes, Pb_coordination)

tilting_distortion, inPlane_distortion, outPlane_distortion = \
        calc_DistortionAngles(octahedral_array, struct)

# Supercell the structure
struct.make_supercell([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
supercell_lattice_vectors = struct.lattice
supercell_species_arr = struct.species
Pb_indexes, I_indexes, Pb_coords, I_coords = getSpeciesInfo(supercell_species_arr)
octahedral_array = get_octahedrals(Pb_indexes, I_indexes, Pb_coordination)
