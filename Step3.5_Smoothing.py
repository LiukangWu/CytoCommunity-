from sklearn.neighbors import NearestNeighbors
from collections import Counter
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sci_palettes
import os
import shutil
import datetime
import matplotlib as mpl
import random

mpl.rcParams['pdf.fonttype'] = 42     # make text in plot editable in AI.
#print(sci_palettes.PALETTES.keys())     # used for checking all color schemes of different journals.
#sci_palettes.PALETTES["d3_category20"]     # see detailed color code
sci_palettes.register_cmap("d3_category20")    # register a specific palette for TCN coloring. "d3_category20" for 20 colors.

# Hyperparameters
KNN_K = 20
Num_TCN = 10
smoothing_range = 50  # 平滑区域大小


InputFolderName = "./Postx-Dx_KNN_Input/"
LastStep_OutputFolderName = "./Step3_Output/"
ThisStep_OutputFolderName = "./Step3.5_Output_Smoothed/"
PatchFolderName = "./Step0_Output/"

if os.path.exists(ThisStep_OutputFolderName):
    shutil.rmtree(ThisStep_OutputFolderName)
os.makedirs(ThisStep_OutputFolderName)

OutputFolderName_1 = ThisStep_OutputFolderName + "TCN_Plot/"
os.mkdir(OutputFolderName_1)
OutputFolderName_2 = ThisStep_OutputFolderName + "Smooth_Plot/"
os.mkdir(OutputFolderName_2)
OutputFolderName_3 = ThisStep_OutputFolderName + "ResultTable_File/"
os.mkdir(OutputFolderName_3)


# Import image name list.
Region_filename = InputFolderName + "ImageNameList.txt"
region_name_list = pd.read_csv(
        Region_filename,
        sep="\t",  # tab-separated
        header=None,  # no heading row
        names=["Image"],  # set our own names for the columns
    )

unique_cell_type_df = pd.read_csv(
        "./Step1_Output/UniqueCellTypeList.txt",
        sep="\t",  # tab-separated
        header=None,  # no heading row
        names=["UniqueCellType"],  # set our own names for the columns
    )
UniqueCellType_vec = unique_cell_type_df['UniqueCellType'].values.tolist()

## Initialize a TCN code list used for matching color palettes across different TCN plots.
UniqueTCN_vec = list(range(1, Num_TCN+1))
UniqueTCN_vec = [str(element) for element in UniqueTCN_vec]


for graph_index in range(len(region_name_list)):

    region_name = region_name_list.Image[graph_index]

    target_graph_map_Merged = pd.DataFrame()
    ## Import target cellular spatial graph x/y coordinates.
    for i in range(2):
        for j in range(2):
            file_prefix = f'Patch_{i}_{j}'
            PatchCoord_filename = PatchFolderName + file_prefix + "-" + region_name + "_Coordinates.txt"
            x_y_coordinates = pd.read_csv(
                PatchCoord_filename,
                sep="\t",  # tab-separated
                header=None,  # no heading row
                names=["x_coordinate", "y_coordinate"],  # set our own names for the columns
            )
            target_graph_map = x_y_coordinates

            ## Import cell type label.
            CellType_filename = PatchFolderName + file_prefix + "-" + region_name + "_CellTypeLabel.txt"
            cell_type_label = pd.read_csv(
                CellType_filename,
                sep="\t",  # tab-separated
                header=None,  # no heading row
                names=["cell_type"],  # set our own names for the columns
            )
            ## Add cell type labels to target graph x/y coordinates.
            target_graph_map["Cell_Type"] = cell_type_label.cell_type

            ## Import the final TCN labels to target graph map.
            MajorityVoting_FileName = LastStep_OutputFolderName + "ImageCollection/" + file_prefix + "-" + region_name + "/TCNLabel_MajorityVoting.csv"
            target_graph_map["TCN_Label"] = np.loadtxt(MajorityVoting_FileName, dtype='int', delimiter=",")

            ## Merge all patches for each original image.
            target_graph_map_Merged = pd.concat([target_graph_map_Merged, target_graph_map],
                                                ignore_index=True)  # "ignore_index=True" will re-generate a new index.

    ## Converting integer list to string list for making color scheme discrete.
    target_graph_map_Merged.TCN_Label = target_graph_map_Merged.TCN_Label.astype(str)
    ## Below is for matching color palettes across different TCN plots, which is quite useful for supervised tasks.
    target_graph_map_Merged["TCN_Label"] = pd.Categorical(target_graph_map_Merged["TCN_Label"], UniqueTCN_vec)
    ## Below is for matching color palettes across different cell type plots, which is quite useful for supervised tasks.
    target_graph_map_Merged["Cell_Type"] = pd.Categorical(target_graph_map_Merged["Cell_Type"], UniqueCellType_vec)

    # -----------------------------------------获取边界处的smooth_cells------------------------------------------------- #
    # 创建target_graph_map_Merged_Smooth，存储平滑后的结果
    target_graph_map_Merged_Smooth = target_graph_map_Merged
    target_graph_map_Merged_Smooth["Smooth_Label"] = "0"  # 标记smooth_cells(初始化为0)，用于可视化

    x_coords = target_graph_map_Merged['x_coordinate'].values
    y_coords = target_graph_map_Merged['y_coordinate'].values
    coordinates_array = np.column_stack((x_coords, y_coords))  # target_graph_map_Merged中的坐标顺序和原始input不同!!!

    # 计算x和y的最值
    x_min, y_min = np.min(coordinates_array, axis=0)
    x_max, y_max = np.max(coordinates_array, axis=0)
    # 计算边界
    x_mid = (x_min + x_max) / 2
    y_mid = (y_min + y_max) / 2

    # 确定smooth cells
    smooth_cells = []

    for index, (x, y) in enumerate(coordinates_array):
        if abs(x - x_mid) <= smoothing_range or abs(y - y_mid) <= smoothing_range:
            smooth_cells.append(index)

    # --------------------------------------更新smooth_cells的TCN_Label---------------------------------------------- #
    for cell in smooth_cells:
        target_graph_map_Merged_Smooth.loc[cell, "Smooth_Label"] = "1"

        # 查找平滑细胞的KNN
        K = KNN_K
        nbrs = NearestNeighbors(n_neighbors=K + 1)  # 加1表示初始包含自身
        nbrs.fit(coordinates_array)
        distances, indices = nbrs.kneighbors(coordinates_array[cell].reshape(1, -1))

        KNN_neighbors = indices[0][1:]

        # for i in KNN_neighbors:
        #     target_graph_map_Merged_Smooth.loc[i, "Smooth_Label"] = "2"
        # print(f"Node {cell} {K} nearest neighbors: {KNN_neighbors}")

        # 获取K个邻居对应的TCN_label
        neighbor_labels = target_graph_map_Merged.loc[KNN_neighbors, 'TCN_Label'].values  # K个邻居Merge后的label
        new_labels = Counter(neighbor_labels).most_common(1)  # (label,count)列表
        target_graph_map_Merged_Smooth.loc[cell, "TCN_Label"] = new_labels[0][0]

    # -----------------------------------------绘图：突出显示边界处的smooth_cells---------------------------------------------- #
    smooth_plot_dict = {"0": "#beaed4", "1": "#7fc97f", "2": "Red"}
    smooth_plot = sns.scatterplot(x="x_coordinate", y="y_coordinate", data=target_graph_map_Merged_Smooth, hue="Smooth_Label",
                               palette=smooth_plot_dict, alpha=1.0, s=2,
                               legend="full")  # "d3_category20" for 20 colors.
    smooth_plot.spines.right.set_visible(False)
    smooth_plot.spines.left.set_visible(False)
    smooth_plot.spines.top.set_visible(False)
    smooth_plot.spines.bottom.set_visible(False)
    smooth_plot.set(xticklabels=[])  # remove the tick label.
    smooth_plot.set(yticklabels=[])
    smooth_plot.set(xlabel=None)  # remove the axis label.
    smooth_plot.set(ylabel=None)
    smooth_plot.tick_params(bottom=False, left=False)  # remove the ticks.

    plt.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0)
    # Save the CURRENT figure.
    smooth_fig_filename1 = OutputFolderName_2 + "SmoothPlot_" + region_name + ".png"
    plt.savefig(smooth_fig_filename1)
    plt.close()

    # # -----------------------------------------Generate plots------------------------------------------------- #
    ## Plot x/y map with "TCN_Label" coloring.
    TCN_plot = sns.scatterplot(x="x_coordinate", y="y_coordinate", data=target_graph_map_Merged_Smooth, hue="TCN_Label",
                               palette="d3_category20", alpha=1.0, s=2,
                               legend="full")  # "d3_category20" for 20 colors.
    # Hide all four spines
    TCN_plot.spines.right.set_visible(False)
    TCN_plot.spines.left.set_visible(False)
    TCN_plot.spines.top.set_visible(False)
    TCN_plot.spines.bottom.set_visible(False)
    TCN_plot.set(xticklabels=[])  # remove the tick label.
    TCN_plot.set(yticklabels=[])
    TCN_plot.set(xlabel=None)  # remove the axis label.
    TCN_plot.set(ylabel=None)
    TCN_plot.tick_params(bottom=False, left=False)  # remove the ticks.

    plt.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0)
    # Save the CURRENT figure.
    TCN_fig_filename1 = OutputFolderName_1 + "TCN_" + region_name + ".pdf"
    plt.savefig(TCN_fig_filename1)
    TCN_fig_filename2 = OutputFolderName_1 + "TCN_" + region_name + ".png"
    plt.savefig(TCN_fig_filename2)
    plt.close()

    ## Export result dataframe: "target_graph_map_Merged".
    TargetGraph_dataframe_filename = OutputFolderName_3 + "ResultTable_" + region_name + ".csv"
    target_graph_map_Merged_Smooth = target_graph_map_Merged_Smooth.drop('Smooth_Label', axis=1)
    target_graph_map_Merged_Smooth.to_csv(TargetGraph_dataframe_filename, na_rep="NULL", index=False)  # remove row index.

