import numpy as np
from geopy import distance
from shapely.geometry import Point
import geopandas as gpd
from scipy.spatial import cKDTree
import networkx as nx
import osmnx as ox
from collections import defaultdict
from math import atan2, radians, sin, cos, degrees
from matplotlib import pyplot as plt
import igraph as ig

import warnings
warnings.filterwarnings("ignore")



import subprocess

def compute_sps_city(
    city: str,
    lat: float,
    lng: float,
    exp_id: str,
    rfrom: int = 1,
    rto: int = 30,
    rstep: int = 1,
    ncircles: int = 36,
    thdist: float = 0.5,
    perturbations: str = "0-100-200-300-400-500",
    nsectors: int = 8,
    saveroutes: int = 1,
    njobs: int = 5
):
    """
    Compute shortest path stability (SPS) for a given city using compute_sps.py.

    Parameters:
        city (str): City name.
        lat (float): Latitude of city center.
        lng (float): Longitude of city center.
        exp_id (str): Unique identifier for experiment output folder.
        rfrom (int): Start radius for radial sampling in km.
        rto (int): End radius for radial sampling in km.
        rstep (int): Step size for radial sampling in km.
        ncircles (int): Number of samples per circle.
        thdist (float): Minimum distance between samples in km.
        perturbations (str): Dash-separated string of displacement magnitudes (e.g., "25-50-100").
        nsectors (int): Number of angular sectors per perturbation ring.
        saveroutes (int): Whether to save the computed paths (0 or 1).
        njobs (int): Number of parallel jobs.
    """

    cmd = [
        "python", "compute_sps.py",
        "-c", city,
        "--lat", str(lat),
        "--lng", str(lng),
        "-i", exp_id,
        "-f", str(rfrom),
        "-t", str(rto),
        "-s", str(rstep),
        "-n", str(ncircles),
        "-d", str(thdist),
        "-p", perturbations,
        "--nsectors", str(nsectors),
        "-r", str(saveroutes),
        "--njobs", str(njobs)
    ]

    print(f"Running Shortest Path Stability (SPS) computation for {city}...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"Shortest Path Stability (SPS) computation completed successfully for {city}")
    else:
        print(f"Shortest Path Stability (SPS) computation failed for {city}")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)






EARTH_RADIUS = 6371000  # Radius of the Earth in meters


def convert_to_cartesian(lng, lat):
    """
    Convert geographic coordinates (longitude, latitude) to Cartesian coordinates (x, y)
    using an equirectangular projection.

    Parameters:
    -----------
    lng : float
        The longitude in degrees.
    lat : float
        The latitude in degrees.

    Returns:
    --------
    tuple
        A tuple (x, y) representing the Cartesian coordinates in meters, where `x` is the 
        horizontal distance and `y` is the vertical distance from the equator.

    Example:
    --------
    >>> convert_to_cartesian(-73.985428, 40.748817)
    (-8240216.806799497, 4480358.041798708)

    Notes:
    ------
    - This conversion assumes the Earth is a perfect sphere and uses an equirectangular 
      projection, which works best for small distances or areas.
    - The result represents meters relative to the Earth's radius.
    """
    
    lat_rad = np.radians(lat)
    lng_rad = np.radians(lng)
    
    x = EARTH_RADIUS * np.cos(lat_rad) * lng_rad
    y = EARTH_RADIUS * lat_rad
    
    return x, y



def get_node_coords(G, node_id):

    node_lng = G.nodes[node_id]["x"]
    node_lat = G.nodes[node_id]["y"]

    return (node_lng, node_lat)


def fixed_radius_sampling(center, radius, nb_samples):
    """
    Generate a fixed number of equidistant sample points on the circumference of a circle 
    defined by a given center and radius.

    Parameters:
    -----------
    center : tuple
        A tuple containing the coordinates (latitude, longitude) of the center point.
    radius : float
        The radius of the circle in kilometers on which the sample points will be placed.
    nb_samples : int
        The number of sample points to generate on the circle's circumference.

    Returns:
    --------
    list of tuple
        A list of tuples where each tuple represents the coordinates (longitude, latitude) 
        of a sample point on the circle's circumference.
    
    Example:
    --------
    >>> fixed_radius_sampling((40.748817, -73.985428), 5, 12)
    [(-73.985428, 40.803654), (-73.930601, 40.800425), ...]

    Notes:
    ------
    - The function uses geographic coordinates and assumes the Earth is a perfect sphere, 
      so it provides approximate results over long distances.
    - The output list contains `nb_samples` points distributed evenly in a circular pattern 
      around the center point.
    """
    
    res = []

    for theta in range(nb_samples):
        point = distance.distance(kilometers=radius).destination(center, theta * (360/nb_samples))
        res = res + [(point[1], point[0])]

    return res


def perform_sampling(G, r_list, city_center, n_samples_circle, th_distance):#, kd_tree):
    
    
    # GeoPandas representing the sampled points
    gpd_points = {}

    # (lng, lat) representing the sampled points
    sampled_points = {}

    # closest graph node associated with the sampled points
    sampled_nodes = {}

    # sampled nodes coordinates
    sampled_nodes_coordinates = {}

    # list of points that don't have a graph node within th_distance km
    points_outside = {}

    # distance point to closest node
    sampled_node_distance = {}

    # Create a KD-tree
    list_node_ids_kdt = list(G.nodes())
    node_coordinates = [(data['x'], data['y']) for node, data in G.nodes(data=True)]
    # Convert all node coordinates to the Cartesian projection
    cartesian_coordinates = [convert_to_cartesian(lon, lat) for lon, lat in node_coordinates]
    # Create the KDTree using the Cartesian coordinates
    kd_tree = cKDTree(cartesian_coordinates)
        

    for r in r_list:

        points_r_km = fixed_radius_sampling(city_center, r, n_samples_circle)
        sampled_points[r] = points_r_km

        edges_r_km, nodes_r_km, points_outside_r, sampled_node_distance_r = [], [], [], []

        for (lng, lat) in points_r_km:

            # retrieve closest Graph node
            query_point_cartesian = convert_to_cartesian(lng, lat)
            dist, nearest_node_index = kd_tree.query(query_point_cartesian)
            node_id_nearest_node = list_node_ids_kdt[nearest_node_index]

            node_data = G.nodes[node_id_nearest_node]
            lng_node = node_data["x"]
            lat_node = node_data["y"]
            
            haversine_distance_km = ox.distance.great_circle(lat, lng, lat_node, lng_node)/1e3

            if haversine_distance_km <= th_distance:
                nodes_r_km.append(node_id_nearest_node)

            else:
                points_outside_r.append([lng, lat])
                nodes_r_km.append({})

            sampled_node_distance_r.append(haversine_distance_km)
            sampled_nodes_coordinates[node_id_nearest_node] = (lng_node, lat_node)
            

        sampled_nodes[r] = nodes_r_km
        points_outside[r] = points_outside_r
        sampled_node_distance[r] = sampled_node_distance_r
        

        geometry = [Point(xy) for xy in points_r_km]
        gdf_points = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326")

        gpd_points[r] = gdf_points
        
        
        sampling_info = {}

        sampling_info["sampling_parameters"] = {"r_list": list(r_list), 
                                                "n_samples_circle": n_samples_circle, 
                                                "th_distance": th_distance}

        sampling_info["sampled_points"] = sampled_points
        sampling_info["gpd_points"] = gpd_points
        
        sampling_info["sampled_nodes"] = sampled_nodes
        sampling_info["points_outside"] = points_outside
        sampling_info["sampled_node_distance"] = sampled_node_distance
        sampling_info["sampled_nodes_coordinates"] = sampled_nodes_coordinates

        
    return sampling_info




def find_nodes_in_radius(node_id, G, kd_tree, list_node_ids_kdt, radius_m_min=150, radius_m_max=200):
    
    # Get node destination data
    coords_dest = get_node_coords(G, node_id)
    
    # Convert to Cartesian coordinates
    query_point_cartesian = convert_to_cartesian(*coords_dest)
    
    # Query for nodes within the maximum radius
    indices_max = kd_tree.query_ball_point(query_point_cartesian, radius_m_max, return_sorted=True)
    
    # Query for nodes within the minimum radius
    indices_min = kd_tree.query_ball_point(query_point_cartesian, radius_m_min, return_sorted=True)
    
    # Find nodes within the ring (between min and max radii)
    indices_ring = list(set(indices_max) - set(indices_min))
    
    # Get the corresponding node IDs
    closest_nodes_kd = [list_node_ids_kdt[idx] for idx in indices_ring]
    
    # Return the nodes in the ring
    return closest_nodes_kd




def calculate_initial_compass_bearing(pointA, pointB):
    """
    Calculate the bearing between two points with input as (longitude, latitude).
    The formula used is the following:
        θ = atan2(sin(Δλ) * cos(φ2), cos(φ1) * sin(φ2) − sin(φ1) * cos(φ2) * cos(Δλ))    
    :param pointA: Tuple of longitude and latitude for the first point (lng1, lat1)
    :param pointB: Tuple of longitude and latitude for the second point (lng2, lat2)
    :return: The bearing in degrees
    """
    lon1, lat1 = pointA
    lon2, lat2 = pointB
    
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    diff_long = radians(lon2 - lon1)

    x = sin(diff_long) * cos(lat2)
    y = cos(lat1) * sin(lat2) - (sin(lat1) * cos(lat2) * cos(diff_long))

    initial_bearing = atan2(x, y)

    # Convert from radians to degrees and normalize the bearing
    initial_bearing = degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing


'''
def OLD_filter_nodes_by_quadrant(G, central_node, list_node_ids, num_quadrants):

    central_point = get_node_coords(G, central_node)

    # Dictionary to hold nodes per quadrant
    nodes_in_quadrants = {i: [] for i in range(num_quadrants)}

    # Iterate over the node IDs and assign each node to a quadrant/sector
    for node_id in list_node_ids:
        # Get the coordinates of the node from the graph G
        node_data = G.nodes[node_id]
        node_coords = (node_data['x'], node_data['y'])  # Assuming 'x' is longitude and 'y' is latitude
        
        # Convert coordinates to polar coordinates relative to the central point
        angle = calculate_initial_compass_bearing(central_point, node_coords)
        
        # Calculate which quadrant/sector the point belongs to
        angle = (angle + 360) % 360  # Normalize angle to range [0, 360)
        sector = int(angle // (360 / num_quadrants))  # Assign to sector based on angle

        # Add the node to the corresponding sector
        nodes_in_quadrants[sector].append(node_id)

    # Select one node per quadrant
    selected_node_ids = []
    for sector, sector_nodes in nodes_in_quadrants.items():
        if sector_nodes:
            # Select the first node in each sector (or you can apply another selection criteria)
            selected_node_ids.append(sector_nodes[0])

    return selected_node_ids
'''

def get_middle_angle_of_sector(sector, num_quadrants):
    """ Calculate the middle angle of a sector. """
    sector_width = 360 / num_quadrants
    middle_angle = sector * sector_width + (sector_width / 2)
    return middle_angle

def calculate_angle_difference(angle1, angle2):
    """ Calculate the smallest difference between two angles. """
    diff = abs(angle1 - angle2) % 360
    if diff > 180:
        diff = 360 - diff
    return diff



def filter_nodes_by_quadrant(G, central_node, list_node_ids, num_quadrants):
    """
    Filters nodes by dividing the area into quadrants (or sectors) and selecting the most central node in each sector.

    The function computes the angular sector each node falls into based on its position relative to a central point. 
    For each sector, it selects the node closest to the middle angle of that sector, ensuring balanced distribution 
    of nodes across all quadrants/sectors.

    Args:
        central_point (tuple): A tuple of the central point's longitude and latitude (lon, lat).
        list_node_ids (list): A list of node IDs to filter.
        num_quadrants (int): The number of quadrants (or sectors) to divide the area into.
        G (networkx.Graph): The graph containing the node data (including longitude and latitude).

    Returns:
        list: A list of filtered node IDs, one node per quadrant/sector, where each selected node is closest 
        to the middle of its respective sector.
    """
    
    central_point = get_node_coords(G, central_node)

    # Dictionary to hold nodes per quadrant
    nodes_in_quadrants = {i: [] for i in range(num_quadrants)}

    # Dictionary to store the angle of each node
    node_angles = {}

    # Iterate over the node IDs and assign each node to a quadrant/sector
    for node_id in list_node_ids:
        # Get the coordinates of the node from the graph G
        node_data = G.nodes[node_id]
        node_coords = (node_data['x'], node_data['y'])  # Assuming 'x' is longitude and 'y' is latitude
        
        # Convert coordinates to polar coordinates relative to the central point
        angle = calculate_initial_compass_bearing(central_point, node_coords)
        
        # Calculate which quadrant/sector the point belongs to
        angle = (angle + 360) % 360  # Normalize angle to range [0, 360)
        sector = int(angle // (360 / num_quadrants))  # Assign to sector based on angle

        # Add the node to the corresponding sector and store its angle
        nodes_in_quadrants[sector].append(node_id)
        node_angles[node_id] = angle

    # Select one node per quadrant based on proximity to the middle of the sector
    selected_node_ids = []
    for sector, sector_nodes in nodes_in_quadrants.items():
        if sector_nodes:
            # Calculate the middle angle of the sector
            middle_angle = get_middle_angle_of_sector(sector, num_quadrants)
            
            # Find the node closest to the middle angle
            closest_node = min(sector_nodes, key=lambda node_id: calculate_angle_difference(node_angles[node_id], middle_angle))
            
            # Add the closest node to the selected nodes
            selected_node_ids.append(closest_node)

    return selected_node_ids


def weighted_jaccard_similarity(list1, list2, dict_weights):
    
    set1 = set(list1)
    set2 = set(list2)
    
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    
    # Calculate total length of edges in intersection and union
    intersection_length = sum(dict_weights[edge] for edge in intersection)
    union_length = sum(dict_weights[edge] for edge in union)
    
    weighted_jaccard = intersection_length / union_length if union_length > 0 else 0
    
    return weighted_jaccard

def get_pair_list(route):
    consecutive_pairs = [(route[i], route[i+1]) for i in range(len(route)-1)]
    return consecutive_pairs


def initialize_nested_dict(d, keys):
    """Helper function to ensure nested dictionaries are initialized."""
    current = d
    for key in keys:
        if key not in current:
            current[key] = {}
        current = current[key]
    return current




def route_from_ig_to_nx(route, dict_ig_to_nx):
    return [dict_ig_to_nx[edge] for edge in route]
    


def parallel_func_compute_paths_r_IG(G, r, sampled_elements, perturbation_factors, n_sectors, spi_dict):

    # Optimization: uses IGraph to compute shortest paths to a target list of nodes
    # rather than building the entire shortest path tree, significantly reducing computation time.
    
    # Convert the nx network to ig
    G_ig = ig.Graph.from_networkx(G)

    # dict for quick ID lookup
    dict_nx_to_ig = {}
    dict_ig_to_nx = {}
    
    for ind_n, n in enumerate(G_ig.vs()):
        dict_nx_to_ig[n["_nx_name"]] = ind_n
        dict_ig_to_nx[ind_n] = n["_nx_name"]

    # Create a KD-tree
    list_node_ids_kdt = list(G.nodes())
    node_coordinates = [(data['x'], data['y']) for node, data in G.nodes(data=True)]
    # Convert all node coordinates to the Cartesian projection
    cartesian_coordinates = [convert_to_cartesian(lon, lat) for lon, lat in node_coordinates]
    # Create the KDTree using the Cartesian coordinates
    kd_tree = cKDTree(cartesian_coordinates)

    
    n_samples_circle = len(sampled_elements[r])
    dict_tmp = {}

    # Pre-compute the target nodes
    
    for ind_src in range(n_samples_circle):
        
        node_src = sampled_elements[r][ind_src]
        
        if not node_src:  # Skip if node_src is an empty dictionary
            continue

        list_target_nodes = []
        
        # Compute the target nodes (destinations + displacements)
        for ind_dest in range(n_samples_circle):
            if ind_src == ind_dest:
                continue
            
            node_dest = sampled_elements[r][ind_dest]
            
            if not node_dest:  # Skip if node_dest is an empty dictionary
                continue

            list_target_nodes.append(node_dest)

            # perturbations
            for min_r_pert, max_r_pert in perturbation_factors:
                # Create the perturbation of the nodes
                nodes_in_radius = find_nodes_in_radius(node_dest, G, kd_tree, 
                                                       list_node_ids_kdt, 
                                                       radius_m_min=min_r_pert, 
                                                       radius_m_max=max_r_pert)
                
                perturbed_nodes = filter_nodes_by_quadrant(G, node_dest, nodes_in_radius, n_sectors)
    
                list_target_nodes += perturbed_nodes
                
        
        for weight in ["length", "travel_time"]:
            
            # Compute the shortest paths using IGraph!
            list_target_nodes_ig = [dict_nx_to_ig[id_nx] for id_nx in list_target_nodes]
            node_src_ig = dict_nx_to_ig[node_src]
            sp_tree_ig = G_ig.get_shortest_paths(node_src_ig, list_target_nodes_ig, weights=weight, mode="OUT", output="vpath")

            #print("Target nodes:", list_target_nodes)
            
            sp_tree_clean = {tg: route_from_ig_to_nx(sp_tree_ig[ind_tg], dict_ig_to_nx) for ind_tg, tg in enumerate(list_target_nodes)}

            #print("sp_tree_clean", sp_tree_clean.keys())
            
            for ind_dest in range(n_samples_circle):
                if ind_src == ind_dest:
                    continue

                node_dest = sampled_elements[r][ind_dest]
                if not node_dest:  # Skip if node_dest is an empty dictionary
                    continue

                #if node_dest in sp_tree_clean:
                # if the dest is reacheable
                if len(sp_tree_clean[node_dest]) > 0:
                    
                    sp_original = sp_tree_clean[node_dest]

                    for min_r_pert, max_r_pert in perturbation_factors:
                        # Create the perturbation of the nodes
                        nodes_in_radius = find_nodes_in_radius(node_dest, G, kd_tree, 
                                                               list_node_ids_kdt, 
                                                               radius_m_min=min_r_pert, 
                                                               radius_m_max=max_r_pert)
                        perturbed_nodes = filter_nodes_by_quadrant(G, node_dest, nodes_in_radius, n_sectors)

                        # Find the shortest path from the origin to all the perturbed nodes
                        dict_perturbed_node_route = {}
                        for perturbed_node in perturbed_nodes:
                            perturbed_route = sp_tree_clean[perturbed_node]
                            #if len(perturbed_route)>0:
                            dict_perturbed_node_route[perturbed_node] = perturbed_route
  

                        # Prepare keys and initialize nested dictionaries
                        str_od = f"{ind_src}_{ind_dest}"
                        str_r_pert = f"{min_r_pert}_{max_r_pert}"
                        nested_dict = initialize_nested_dict(dict_tmp, [str_od, weight])
                        
                        # NEW
                        nested_dict["optimal"] = sp_original
                        nested_dict["from"] = node_src
                        nested_dict["to"] = node_dest

                        if "perturbations" not in nested_dict:
                            nested_dict["perturbations"] = {}
                        
                        # Store the results
                        nested_dict["perturbations"][str_r_pert] = dict_perturbed_node_route
                        
                        '''
                        # Store the results
                        nested_dict[str_r_pert] = {
                            "from": node_src,
                            "to": node_dest,
                            "optimal": sp_original,
                            "perturbations": dict_perturbed_node_route
                        }'''
        #else:
        #    print(f"Disconnected: r={r}, src={node_src}, dest={node_dest}")

    spi_dict[r] = dict_tmp

def replace_duplicates(lst):
    seen = set()  # Set to track unique values
    result = []   # List to store the final result
    for item in lst:
        # Check if the item is already in the seen set or is a dictionary
        if isinstance(item, dict) or item in seen:
            result.append({})  # Append {} for duplicates
        else:
            seen.add(item)     # Add the unique item to the set
            result.append(item)  # Append the first occurrence
    return result

# this is the one that is executed
def NEW_parallel_func_compute_paths_r_IG(G, r, sampled_elements, perturbation_factors, n_sectors, spi_dict):

    # Convert the nx network to ig
    G_ig = ig.Graph.from_networkx(G)
    
    # dict for quick ID lookup
    dict_nx_to_ig = {}
    dict_ig_to_nx = {}
    
    for ind_n, n in enumerate(G_ig.vs()):
        dict_nx_to_ig[n["_nx_name"]] = ind_n
        dict_ig_to_nx[ind_n] = n["_nx_name"]
    
    # Create a KD-tree
    list_node_ids_kdt = list(G.nodes())
    node_coordinates = [(data['x'], data['y']) for node, data in G.nodes(data=True)]
    # Convert all node coordinates to the Cartesian projection
    cartesian_coordinates = [convert_to_cartesian(lon, lat) for lon, lat in node_coordinates]
    # Create the KDTree using the Cartesian coordinates
    kd_tree = cKDTree(cartesian_coordinates)

    # replace duplicates with -> {} (to fix the bug where origin and destination are the same)
    sampled_elements[r] = replace_duplicates(sampled_elements[r])
    
    n_samples_circle = len(sampled_elements[r])
    dict_tmp = {}

    # compute the perturbations for each node

    dict_dest_to_perturbed = {}
    
    for ind_node in range(n_samples_circle):
    
        list_perturbed_nodes = []
        
        node_name = sampled_elements[r][ind_node]
    
        if node_name == {}:
            continue
            
        for min_r_pert, max_r_pert in perturbation_factors:
            
            # Create the perturbation of the nodes
            nodes_in_radius = find_nodes_in_radius(node_name, G, kd_tree, list_node_ids_kdt, 
                                                   radius_m_min=min_r_pert, radius_m_max=max_r_pert)
            
            perturbed_nodes = filter_nodes_by_quadrant(G, node_name, nodes_in_radius, n_sectors)
    
            list_perturbed_nodes += perturbed_nodes
    
        dict_dest_to_perturbed[node_name] = list_perturbed_nodes

    # FOR DIFFERENT WEIGHTS...
    
    for weight in ["length", "travel_time"]:
    
        dict_sp_trees = {}
        
        for ind_node in range(n_samples_circle):
            
            node_name = sampled_elements[r][ind_node]
        
            if node_name == {}:
                continue
            
            # build all the possible targets, that are the perturbations of itself + all the destinations and their perturbations.
            list_target_nodes = list(k for k in dict_dest_to_perturbed.keys() if k != node_name) + [item for sublist in dict_dest_to_perturbed.values() for item in sublist]
        
            sp_tree_node = get_shortest_path_tree(G_ig, dict_nx_to_ig, dict_ig_to_nx, node_name, list_target_nodes, weight=weight)
            
            dict_sp_trees[node_name] = sp_tree_node
    
        
        # START COMPUTATION for each OD
    
        for ind_src in range(n_samples_circle):
                
            node_src = sampled_elements[r][ind_src]
            
            if not node_src:  # Skip if node_src is an empty dictionary
                continue
        
            for ind_dest in range(n_samples_circle):
                if ind_src == ind_dest:
                    continue
        
                node_dest = sampled_elements[r][ind_dest]
                if not node_dest:  # Skip if node_dest is an empty dictionary
                    continue
        
                # Shortest path three Source (for A->B and A->B's)
                sp_tree_source = dict_sp_trees[node_src]
                    
                if len(sp_tree_source[node_dest]) > 0:
                    
                    # the shortest path from the origin to the destination location
                    sp_original = sp_tree_source[node_dest]
        
                    # Shortest path three Source (B->B's)
                    sp_tree_dest = dict_sp_trees[node_dest]
        
        
                    for min_r_pert, max_r_pert in perturbation_factors:
                        
                        # Create the perturbation of the nodes
                        nodes_in_radius = find_nodes_in_radius(node_dest, G, kd_tree, 
                                                               list_node_ids_kdt, 
                                                               radius_m_min=min_r_pert, 
                                                               radius_m_max=max_r_pert)
                        
                        perturbed_nodes = filter_nodes_by_quadrant(G, node_dest, nodes_in_radius, n_sectors)
        
                        # Find the shortest path from the origin to all the perturbed nodes
                        # A -> B'
                        dict_routes_A_pertB = {}
                        for perturbed_node in perturbed_nodes:
                            perturbed_route = sp_tree_source[perturbed_node]
                            dict_routes_A_pertB[perturbed_node] = perturbed_route
                    
                        # Find the shortest path from the destination to all the perturbed nodes
                        # B -> B'
                        dict_routes_B_pertB = {}
                        for perturbed_node in perturbed_nodes:
                            perturbed_route = sp_tree_dest[perturbed_node]
                            dict_routes_B_pertB[perturbed_node] = perturbed_route
        
        
                        # Prepare keys and initialize nested dictionaries
                        str_od = f"{ind_src}_{ind_dest}"
                        str_r_pert = f"{min_r_pert}_{max_r_pert}"
                        nested_dict = initialize_nested_dict(dict_tmp, [str_od, weight])
                        
                        # NEW
                        nested_dict["optimal"] = sp_original
                        nested_dict["from"] = node_src
                        nested_dict["to"] = node_dest
        
                        if "perturbations" not in nested_dict:
                            nested_dict["perturbations"] = {}
                        
                        # Store the results
                        nested_dict["perturbations"][str_r_pert] = {"routes_A_pertB": dict_routes_A_pertB, 
                                                                    "routes_B_pertB": dict_routes_B_pertB}
                        
                        '''
                        # Store the results
                        nested_dict[str_r_pert] = {
                            "from": node_src,
                            "to": node_dest,
                            "optimal": sp_original,
                            "perturbations": dict_perturbed_node_route
                        }'''
        #else:
        #    print(f"Disconnected: r={r}, src={node_src}, dest={node_dest}")
        
        spi_dict[r] = dict_tmp



def parallel_compute_measures_spi_r(graph, radius, routing_elements, dict_spi_routes, dict_shared_results):
    """
    Compute weighted Jaccard similarity between optimal and perturbed routes for a given radius
    and calculate the total path lengths for different perturbation types.
    
    Parameters:
    graph: NetworkX graph with route information.
    radius: Radius for which the routes are being computed.
    routing_elements: Dictionary of routing elements for each radius.
    dict_spi_routes: Dictionary containing route information for each OD pair and perturbation.
    dict_shared_results: Shared dictionary to store the computed results.
    """

    # Initialize results dictionary
    results_for_radius = {}
    
    # Precompute edge lengths from the graph
    edge_lengths = {(u, v): data["length"] for u, v, data in graph.edges(data=True)}
    
    # Convert radius to float and get the number of OD pairs
    radius = float(radius)
    num_od_pairs = len(routing_elements[radius])
    
    # Iterate through each OD pair (source and destination indices)
    for source_idx in range(num_od_pairs):
        for dest_idx in range(num_od_pairs):
            
            # Create a unique key for the OD pair
            od_pair_key = f"{source_idx}_{dest_idx}"
    
            # Check if routes exist for this OD pair
            if od_pair_key in dict_spi_routes[radius]:
    
                for weight_type, info_routes in dict_spi_routes[radius][od_pair_key].items():
    
                    # Get the optimal route and convert it to an edge list
                    optimal_route = info_routes["optimal"]
                    optimal_edge_list = get_pair_list(optimal_route)
                    optimal_path_length = np.sum([edge_lengths[edge] for edge in get_pair_list(optimal_route)])
    
                    for perturbation_level in info_routes["perturbations"]:
    
                        dict_res_pert_level = {}
                        
                        for type_perturbation_od in ["routes_A_pertB", "routes_B_pertB"]:
                    
                            perturbation_routes = info_routes["perturbations"][perturbation_level][type_perturbation_od]
                        
                            if type_perturbation_od == "routes_A_pertB":
                                
                                # jaccard
                                perturbation_results_jaccard = compute_perturbation_jaccards(optimal_edge_list, perturbation_routes, edge_lengths)
                                dict_res_pert_level = perturbation_results_jaccard
                                
                                # cost
                                for node, route in perturbation_routes.items():
                                    total_path_length = np.sum([edge_lengths[edge] for edge in get_pair_list(route)])
                                    dict_res_pert_level[node]["total_length_A_pertB"] = round(total_path_length, 3)
                    
                            if type_perturbation_od == "routes_B_pertB":
                                # cost
                                for node, route in perturbation_routes.items():
                                    total_path_length = np.sum([edge_lengths[edge] for edge in get_pair_list(route)])
                                    dict_res_pert_level[node]["total_length_B_pertB"] = round(total_path_length, 3)


                        # Store the results in the result dictionary
                        if od_pair_key not in results_for_radius:
                            results_for_radius[od_pair_key] = {}
                        if weight_type not in results_for_radius[od_pair_key]:
                            results_for_radius[od_pair_key][weight_type] = {}
                        
                        results_for_radius[od_pair_key][weight_type][perturbation_level] = dict_res_pert_level

                    results_for_radius[od_pair_key][weight_type]["total_length_A_B"] = round(optimal_path_length, 3)

    # Store results for this radius in the shared dictionary
    dict_shared_results[radius] = results_for_radius


def compute_perturbation_jaccards(optimal_edge_list, perturbations, edge_lengths):
    """
    Compute the weighted Jaccard similarity for all perturbed routes against the optimal route.
    """
    jaccard_results = {}
    for perturbed_node, perturbed_route in perturbations.items():
        perturbed_edge_list = get_pair_list(perturbed_route)
        if len(perturbed_route)>0:
            jaccard_similarity = weighted_jaccard_similarity(optimal_edge_list, perturbed_edge_list, edge_lengths)
        else:
            jaccard_similarity = -1
            
        jaccard_results[perturbed_node] = {"w_jaccard": jaccard_similarity}
    return jaccard_results


def compare_shortest_vs_fastest(dict_spi_routes, radius, od_pair_key, perturbation_level, results_for_radius, edge_lengths):
    """
    Compare the weighted Jaccard similarity between the shortest and fastest routes.
    """
    # Retrieve the shortest and fastest routes
    shortest_route = dict_spi_routes[radius][od_pair_key]["length"]["optimal"]
    fastest_route = dict_spi_routes[radius][od_pair_key]["travel_time"]["optimal"]

    # Convert to edge lists
    shortest_edge_list = get_pair_list(shortest_route)
    fastest_edge_list = get_pair_list(fastest_route)

    # Compute weighted Jaccard similarity
    jaccard_similarity = weighted_jaccard_similarity(shortest_edge_list, fastest_edge_list, edge_lengths)

    # Store the comparison result
    results_for_radius[od_pair_key]["short_vs_fast"] = {"w_jaccard": jaccard_similarity}


def get_shortest_path_tree(G_ig, dict_nx_to_ig, dict_ig_to_nx, node_src, list_target_nodes, weight):
    
    # Convert target nodes and source node to igraph node ids
    list_target_nodes_ig = [dict_nx_to_ig[id_nx] for id_nx in list_target_nodes]
    node_src_ig = dict_nx_to_ig[node_src]
    
    # Get the shortest paths using igraph's get_shortest_paths method
    sp_tree_ig = G_ig.get_shortest_paths(node_src_ig, list_target_nodes_ig, weights=weight, mode="OUT", output="vpath")
    
    # Convert the paths from igraph format back to networkx format
    sp_tree_clean = {tg: route_from_ig_to_nx(sp_tree_ig[ind_tg], dict_ig_to_nx) for ind_tg, tg in enumerate(list_target_nodes)}
    
    return sp_tree_clean




