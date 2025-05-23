import os
import osmnx as ox
import numpy as np
import json
import argparse
from matplotlib import pyplot as plt
from collections import defaultdict

import multiprocessing
import pandas as pd
import gzip

import lz4.frame
import msgpack

import concurrent
from concurrent.futures import ProcessPoolExecutor, as_completed

from sps_utils import perform_sampling, get_pair_list
from sps_utils import NEW_parallel_func_compute_paths_r_IG, parallel_compute_measures_spi_r

# ================================
# Argument Parsing
# ================================


def parse_arguments():
    """Parse command-line arguments with detailed help descriptions for route generation and analysis."""
    parser = argparse.ArgumentParser(description="Generate routes and evaluate spatial sampling around city centers.")

    # Required Arguments
    parser.add_argument('-c', '--city', required=True,
                        help="City name used to identify and load the corresponding road network.")
    parser.add_argument('-i', '--identifier', type=str, required=True,
                        help="Unique experiment identifier used for result storage and logging.")

    # Optional Arguments
    parser.add_argument('--lat', type=float, required=True,
                        help="Latitude of the city center.")
    parser.add_argument('--lng', type=float, required=True,
                        help="Longitude of the city center.")
    parser.add_argument('--njobs', type=int, default=5,
                        help="Number of parallel jobs to use during computation. Default is 5.")
    
    parser.add_argument('-f', '--rfrom', type=int, default=1, 
                        help="Starting radius (km) for radial sampling.")
    parser.add_argument('-t', '--rto', type=int, default=30, 
                        help="Ending radius (km) for radial sampling.")
    parser.add_argument('-s', '--rstep', type=int, default=1, 
                        help="Radius step (km) for concentric circles in radial sampling.")
    
    parser.add_argument('-d', '--thdist', type=float, default=0.5,
                        help="Minimum distance threshold (in km) between sampled points. Default is 0.5.")
    parser.add_argument('-n', '--ncircles', type=int, default=36, 
                        help="Number of samples per circle for radial sampling.")

    parser.add_argument('-p', '--perturbations', type=str, default=None,
    help="Dash-separated list of displacement boundaries in meters (e.g., '0-100-200'). "
         "Defines the edges of perturbation intervals, e.g., '0-100-200' → [[0, 100], [100, 200]].")


    parser.add_argument('--nsectors', type=int, default=8,
        help="Number of angular sectors used to select perturbed destinations within each displacement ring (default: 8).")

    
    parser.add_argument('-r', '--saveroutes', type=int, default=0,
                        help="Whether to save generated routes (1) or not (0). Default is 0 (do not save).")

    return parser.parse_args()



# PARSE THE ARGUMENTS
args = parse_arguments()


# City and network parameters
city = args.city
city_center = (args.lat, args.lng)
exp_id = args.identifier
network_type = "drive"
save_routes = args.saveroutes

N_jobs = args.njobs

# Radius parameters (in km and m)
min_radius_km = args.rfrom
max_radius_km = args.rto
radius_step_km = args.rstep
radius_city_m = (max_radius_km+1) * 1000  # Buffer for the city's maximum radius in meters

n_samples_circle = args.ncircles
th_distance_km = args.thdist


# Path Stability: parameters
if args.perturbations:
    # Split the string into a list of integers
    list_numbers = list(map(int, args.perturbations.split('-')))
    perturbation_factors = get_pair_list(list_numbers) 
else:
    perturbation_factors = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500]]

n_sectors = args.nsectors


# output folder
output_folder = f"../data/results/{city}_{exp_id}/"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


print("\n" + "="*40)
print("         Experiment Configuration         ")
print("="*40)
print(f"City:                       {args.city}")
print(f"City Center Coordinates:    ({args.lat}, {args.lng})")
print(f"Experiment ID:              {args.identifier}")
print("-" * 40)
print(f"Radii (min → max, step):    {args.rfrom} km → {args.rto} km, step {args.rstep} km")
print(f"Samples per Radius Circle:  {args.ncircles}")
print(f"Threshold Distance [km]:    {args.thdist}")
print(f"Angular Sectors (k):        {args.nsectors}")
print("-" * 40)
perturbs = args.perturbations if args.perturbations else "None"
print(f"Displacement Magnitudes:    {perturbs}")
print("-" * 40)
print(f"Save Routes:                {'Yes' if args.saveroutes else 'No'}")
print(f"Parallel Jobs:              {args.njobs}")
print("="*40)



# 1. Load or Download the road network if not already present
network_file = f"../data/road_networks/{city}_{network_type}_{radius_city_m}.graphml.gz"
uncompressed_network_file = f"../data/road_networks/{city}_{network_type}_{radius_city_m}.graphml"

if os.path.exists(network_file):
    # If the compressed network file exists, decompress and load it
    print(f"Compressed network found: {network_file}")
    print(f"Decompressing and loading network...")
    print("-" * 40)
    with gzip.open(network_file, 'rb') as f_in:
        with open(uncompressed_network_file, 'wb') as f_out:
            f_out.write(f_in.read())
    
    # Load the graph from the decompressed file
    G = ox.load_graphml(uncompressed_network_file)

    # Optionally, remove the uncompressed file after loading
    os.remove(uncompressed_network_file)

else:
    
    print(f"Network file not found: {network_file}")
    print(f"Downloading road network for {city}...")
    print(f"Using center: {city_center} with radius {radius_city_m} meters")
    print("-" * 40)
    G = ox.graph_from_point(city_center, dist=radius_city_m, network_type=network_type, simplify=True)
    # Impute missing edge speeds and calculate edge travel times with the speed module
    G = ox.speed.add_edge_speeds(G)
    G = ox.speed.add_edge_travel_times(G)

    # Save the network as a regular GraphML file
    ox.save_graphml(G, uncompressed_network_file)

    # Compress the saved GraphML file
    with open(uncompressed_network_file, 'rb') as f_in:
        with gzip.open(network_file, 'wb') as f_out:
            f_out.writelines(f_in)

    # Optionally, remove the uncompressed file after compressing
    os.remove(uncompressed_network_file)
    

# 3. Perform the radial sampling  
r_list = np.arange(min_radius_km, max_radius_km+.01, radius_step_km)
sampling_info = perform_sampling(G, r_list, city_center, n_samples_circle, th_distance_km)

print("\n" + "="*40)
print("         Radial Sampling Information         ")
print("="*40)
print(f"Performing Radial Sampling with the following settings:")
print(f"Radii Range (km):          {min_radius_km} → {max_radius_km} (Step: {radius_step_km})")
print(f"Threshold Distance (km):   {th_distance_km}")
print(f"Samples per Circle:        {n_samples_circle}")
print("-" * 40)

print(f"Computed Radii List:       {r_list}")
print("="*40)


# plot the radial sampling
fig, ax = plt.subplots(figsize=(8, 8))
for r in sampling_info["gpd_points"]:
    sampling_info["gpd_points"][r].plot(ax=ax, color="lightblue", markersize=5)
count_outside = 0
for r in sampling_info["points_outside"]:
    for po in sampling_info["points_outside"][r]:
        plt.scatter(po[0], po[1], marker="x", color="red", s=50)
        count_outside+=1
plt.scatter([], [], marker='x', color='red', label=f'# {count_outside}')
plt.legend()
plt.savefig(output_folder+"fig_sampling.pdf", bbox_inches='tight')

# 4. Compute the perturbed shortest/fastest paths (multi-processing)
routing_elements = sampling_info["sampled_nodes"]

print("Launching the experiments...")

manager = multiprocessing.Manager()
spi_dict = manager.dict()
processes = []

for r in routing_elements.keys():

    
    process = multiprocessing.Process(target=NEW_parallel_func_compute_paths_r_IG, 
                                      args=(G, r, routing_elements, perturbation_factors, n_sectors, spi_dict))
    process.start()
    processes.append(process)
    
for p in processes:
    p.join()

for p in processes:
    p.close()

    
dict_spi_routes = dict(spi_dict)


# 4b. Save the generated paths (in case of some errors the slow part is stored)
print("\n" + "="*40)
print("         Route Computation Completed         ")
print("="*40)
print(f"Saving generated routes to file: {output_folder}generated_routes.lz4")
print("-" * 40)
with lz4.frame.open(output_folder+f"generated_routes.lz4", "wb") as f:
    packed_data = msgpack.packb(dict_spi_routes)
    f.write(packed_data)


# 4c. Save the sampling infos, jaccard and the generated paths
_ = sampling_info.pop("gpd_points")

output_file = open(output_folder+"sampling_info.json", "w")
json.dump(sampling_info, output_file)
output_file.close()


# 5. Compute measures on the computed paths (for now Jaccard)
manager = multiprocessing.Manager()
jaccard_dict = manager.dict()
processes = []

for r in routing_elements.keys():
    
    process = multiprocessing.Process(target=parallel_compute_measures_spi_r, 
                                      args=(G, r, routing_elements, dict_spi_routes, jaccard_dict))
    process.start()
    processes.append(process)
    
for p in processes:
    p.join()

for p in processes:
    p.close()

print("\n" + "="*40)
print("         Results Computation Completed         ")
print("="*40)
print(f"Saving results to file: {output_folder}weighted_jaccards.json")
print("-" * 40)

jaccard_dict = dict(jaccard_dict)

# 6 save the dictionary that contains the computed results (jaccard)
output_file = open(output_folder+"weighted_jaccards.json", "w")
json.dump(jaccard_dict, output_file)
output_file.close()



