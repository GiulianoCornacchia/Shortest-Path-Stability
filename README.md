
<div align="center">
  <h1>Shortest Path Stability in Urban Road Networks</h1>
  <h3>A framework to measure how sensitive routes are to small perturbations in destination</h3>
</div>

### Authors:

* Giuliano Cornacchia <sup>1</sup> [<img src="https://img.shields.io/badge/ORCID-0000--0003--2263--7654-brightgreen?logo=orcid&logoColor=white" alt="ORCID" height="16">](https://orcid.org/0000-0003-2263-7654)  
* Mirco Nanni <sup>1</sup> [<img src="https://img.shields.io/badge/ORCID-0000--0003--3534--4332-brightgreen?logo=orcid&logoColor=white" alt="ORCID" height="16">](https://orcid.org/0000-0003-3534-4332)

<sup>1</sup> ISTI-CNR, Pisa, Italy  

____

📄 *Article coming soon*

____

## Built with

![python](https://img.shields.io/badge/Python-3776AB.svg?style=for-the-badge&logo=Python&logoColor=white)
![jupyter](https://img.shields.io/badge/Jupyter-F37626.svg?style=for-the-badge&logo=Jupyter&logoColor=white)
![numpy](https://img.shields.io/badge/NumPy-013243.svg?style=for-the-badge&logo=NumPy&logoColor=white)
![pandas](https://img.shields.io/badge/pandas-150458.svg?style=for-the-badge&logo=pandas&logoColor=white)
![osm](https://img.shields.io/badge/OpenStreetMap-7EBC6F.svg?style=for-the-badge&logo=OpenStreetMap&logoColor=white)

### Requirements

Tested with:

- ![Python](https://img.shields.io/badge/Python-3.9.18-blue)

---

## Overview

This project implements the methodology presented in our article *"The path is the goal: A study on the nature and effects of shortest-path stability under perturbation of destination"* which introduces **Shortest Path Stability** as a measure of how sensitive shortest paths are to small perturbations in destination location.

It allows researchers to:
- Simulate controlled displacements around destinations.
- Measure how often and how drastically routes change.
- Compare the stability of shortest paths across different cities and spatial regions.

---

## What is Shortest Path Stability?

Given a small displacement of your destination (e.g., 50m to the left), does your optimal route change completely? **Shortest Path Stability** quantifies how stable optimal routes are to such minor shifts.

We use:
- **Controlled radial sampling** around a city center.
- **Displacement intervals** and **angular sectors** to perturb destination locations.
- **Weighted Jaccard similarity** to compare the original and perturbed routes.

---

## Outputs

- `generated_routes.lz4`: Compressed dictionary of shortest paths.
- `sampling_info.json`: Metadata on sampled points and OD pairs.
- `weighted_jaccards.json`: Route stability scores.
- `fig_sampling.pdf`: Visualization of radial sampling layout.

---

## Example Usage

```bash
python compute_sps.py \
  --city pisa \
  --lat 43.715906 --lng 10.401866 \
  --identifier sps_pisa \
  --rfrom 1 --rto 10 --rstep 1 \
  --ncircles 6 \
  --thdist 0.5 \
  --perturbations 0-100-200-300 \
  --nsectors 8 \
  --saveroutes 1 \
  --njobs 4
```

---

## Parameters

| Argument           | Description                                                                                 | Required | Default         |
|--------------------|---------------------------------------------------------------------------------------------|----------|-----------------|
| `--city`           | Name of the city                                                                            | ✅       | —               |
| `--lat`, `--lng`   | Latitude and longitude of the city center                                                   | ✅       | —               |
| `--identifier`     | Unique ID for output folder                                                                 | ✅       | —               |
| `--rfrom`, `--rto`, `--rstep` | Radial range and step in km for sampling                                      | ✅       | —               |
| `--ncircles`       | Number of samples per circle                                                                | ✅       | —               |
| `--thdist`         | Minimum distance in km between sampled points                                               | ✅       | —               |
| `--perturbations`  | Dash-separated boundaries for displacement intervals (e.g., `0-100-200`) → `[[0,100],[100,200]]` | ❌       | `0-100-200-300-400-500` |
| `--nsectors`       | Angular sectors to sample within each perturbation ring                                     | ❌       | 8               |
| `--saveroutes`     | Save computed paths (1 = yes, 0 = no)                                                       | ❌       | 0               |
| `--njobs`          | Number of parallel workers                                                                  | ❌       | 5               |

---

## Repository Structure

- `compute_sps.py`: Main script to compute SPS metrics.
- `compute_sps_city.py`: Batch callable function for cities.
- `sps_utils.py`: Utility functions for sampling, path generation, and stability computation.
- `notebooks/`: Jupyter notebooks for running and visualizing experiments.
