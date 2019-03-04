#Utilize the functions to calculate octahedral distortion parameters. 
#Find definitions in octahedral distortion literature search. 

def bond_length_distortion(bond_length_array):
    sum_total = 0
    bond_length_distortion = 0
    for bond_length in bond_length_array:
        sum_total += bond_length
    avg_bond_length = sum_total/(len(bond_length_array))
    print("Average bond length: ", avg_bond_length)
    for bond_length in bond_length_array:
        bond_length_distortion += (((bond_length - avg_bond_length)/avg_bond_length)**2)/(len(bond_length_array))
    return bond_length_distortion

#Distance from center to vertex is h = (a*2^(1/2))/2
# V_octahedron = (2^(1/2)/3)*a^3 , where a = octahedron edge length

def octahedral_elongation(bond_length_array, volume):
    octahedron_edge = (3*volume/(2**(1/2)))**(1/3)
    d_0 = octahedron_edge*(1/(2**(1/2)))
    total = 0
    for bond_length in bond_length_array:
        total += (bond_length/d_0)**2
    octahedral_elongation = total/len(bond_length_array)
    return octahedral_elongation

def average_halide_distance(halide_dist_array):
    total = 0
    for distance in halide_dist_array:
        total += distance
    avg_halide_dist = total/len(halide_dist_array)
    return avg_halide_dist

def halide_distance_deviation(halide_dist_array):
    avg_halide_dist = average_halide_distance(halide_dist_array)
    print("Average halide distance: ", avg_halide_dist)
    result = 0
    for distance in halide_dist_array:
        result += abs(distance - avg_halide_dist)*avg_halide_dist/len(halide_dist_array)
    return result

def volume_discrepancy(halide_dist_array, volume):
    avg_halide_dist = average_halide_distance(halide_dist_array)
    ideal_octahedral = (avg_halide_dist**3)*(2**(1/2))/3
    print("Ideal octahedral volume: ", ideal_octahedral)
    vol_discrepancy = (ideal_octahedral - volume)*100/ideal_octahedral
    return vol_discrepancy

bond_length_array = eval(input("Type in bond lengths in array format (6): "))
volume = eval(input("Type in volume of octahedral: "))
halide_dist_array = eval(input("Type in halide distances in array format (12): "))
print("Bond length distortion: ", bond_length_distortion(bond_length_array))
print("Octahedral elongation: ", octahedral_elongation(bond_length_array, volume))
print("Halide distance deviation: ", halide_distance_deviation(halide_dist_array))
print("Volume discrepancy: ", volume_discrepancy(halide_dist_array, volume))
