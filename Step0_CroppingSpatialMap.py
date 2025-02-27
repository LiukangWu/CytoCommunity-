import numpy as np
import os
import shutil
import pandas as pd
import datetime


# Hyperparameters
InputFolderName = "./Input/"
CellPatchNum = 20000 #细胞分割阈值

# Output folder
ThisStep_OutputFolderName = "./Step0_Output/"
if os.path.exists(ThisStep_OutputFolderName):
    shutil.rmtree(ThisStep_OutputFolderName)
os.makedirs(ThisStep_OutputFolderName)

# Import image name list.
Region_filename = InputFolderName + "ImageNameList.txt"
region_name_list = pd.read_csv(
        Region_filename,
        sep="\t",  # tab-separated
        header=None,  # no heading row
        names=["Image"],  # set our own names for the columns
    )

print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
print("Cropping each image/sample into small patches...")

# Recursive function for splitting regions
def split_region(region_name, coordinates, cell_types, graph_labels, x_range, y_range, prefix="Patch"):
    x_min, x_max = x_range
    y_min, y_max = y_range
    x_mid = (x_min + x_max) / 2
    y_mid = (y_min + y_max) / 2
    
    sub_ranges = [
        ([x_min, x_mid], [y_min, y_mid]),  # Top-left
        ([x_mid, x_max], [y_min, y_mid]),  # Top-right
        ([x_min, x_mid], [y_mid, y_max]),  # Bottom-left
        ([x_mid, x_max], [y_mid, y_max])   # Bottom-right
    ]

    for i, (sub_x_range, sub_y_range) in enumerate(sub_ranges):
        # Filter coordinates within the sub-region
        sub_coordinates = []
        sub_cell_types = []
        for idx, (x, y) in enumerate(coordinates):
            if sub_x_range[0] <= x < sub_x_range[1] and sub_y_range[0] <= y < sub_y_range[1]:
                sub_coordinates.append((x, y))
                sub_cell_types.append(cell_types[idx])
        
        if len(sub_coordinates) == 0:
            continue  # Skip empty regions
        
        # Check if the patch meets the cell count requirement
        if len(sub_coordinates) > CellPatchNum:
            # If too many cells, split further
            new_prefix = f"{prefix}_{i}"
            split_region(region_name, sub_coordinates, sub_cell_types, graph_labels, sub_x_range, sub_y_range, new_prefix)
        else:
            # Save patch data
            patch_name = f"{prefix}_{i}-{region_name}"
            coord_file = f"{ThisStep_OutputFolderName}{patch_name}_Coordinates.txt"
            type_file = f"{ThisStep_OutputFolderName}{patch_name}_CellTypeLabel.txt"
            label_file = f"{ThisStep_OutputFolderName}{patch_name}_GraphLabel.txt"
            
            with open(coord_file, 'w') as cf, open(type_file, 'w') as tf, open(label_file, 'w') as lf:
                for coord in sub_coordinates:
                    cf.write(f"{coord[0]}\t{coord[1]}\n")
                for cell_type in sub_cell_types:
                    tf.write(f"{cell_type}\n")
                for label in graph_labels:
                    lf.write(f"{label}\n")
            
            # Append patch name to the list
            with open(f"{ThisStep_OutputFolderName}ImagePatchNameList.txt", "a") as f0:
                f0.write(f"{patch_name}\n")

# Main loop for processing each image
for graph_index in range(len(region_name_list)):
    print(f"This is image-{graph_index}")
    region_name = region_name_list.Image[graph_index]
    
    # Import target graph x/y coordinates
    GraphCoord_filename = InputFolderName + region_name + "_Coordinates.txt"
    coordinates = []
    with open(GraphCoord_filename, 'r') as file:
        for line in file:
            parts = line.split()
            coordinates.append((float(parts[0]), float(parts[1])))

    # Import cell type info in the target graph
    CellType_filename = InputFolderName + region_name + "_CellTypeLabel.txt"
    cell_types = []
    with open(CellType_filename, 'r') as file:
        for line in file:
            cell_types.append(line.strip())

    # Import target graph label
    GraphLabel_filename = InputFolderName + region_name + "_GraphLabel.txt"
    graph_labels = []
    with open(GraphLabel_filename, 'r') as file:
        for line in file:
            graph_labels.append(line.strip())

    # Get boundary of the entire region
    coordinates_array = np.array(coordinates)
    x_min, y_min = np.min(coordinates_array, axis=0)
    x_max, y_max = np.max(coordinates_array, axis=0)

    # Start recursive splitting
    split_region(region_name, coordinates, cell_types, graph_labels, [x_min, x_max], [y_min, y_max])

print("Step0 done!")
print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
