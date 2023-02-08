import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tabulate import tabulate
import pprint
import csv


data = pd.read_csv(
    'GSE64881_segmentation_at_30000bp.passqc.multibam.txt', delimiter="\t")
# Read only the nps and columns, eliminates unwanted data
nps = data.iloc[:, 3:]  # read all rows, start reading from column 3 on
sum_columns = nps.sum()  # produce a series of the sum of columns
sum_rows = nps.sum(axis=1)  # produce a series of the rows of columns
# Isolation of Hist1 region
# Grab only the hist1 region and create a new dataframe
hist1_df = data.iloc[69714:69795].copy()
# read all rows, start reading from column 3 on
nps_hist_1 = hist1_df.iloc[:, 3:]
hist1_sum_columns = nps_hist_1.sum()  # produce a series of the sum of columns
hist1_sum_rows = nps_hist_1.sum(axis=1)  # produce a series of the sum of rows
# Isolate relevant Nps
# creates new dataframe by using boolean indexing(ex. df = df[specified values])
# return columns where the series value is not equal to zero (ex. series_sum[series_sum != 0])
# .index returns the column names to form the new data frame
nps_hist_1 = nps_hist_1[hist1_sum_columns[hist1_sum_columns != 0].index]
# create a new series of sums for the isolated nps
hist1_sum_columns = nps_hist_1.sum()


def num_windows():
    return len(nps.iloc[0:])
# count the nuclear profiles, all start with F


def count_npfs():
    count = 0
    for column_name in data:
        if column_name[0] == 'F':
            count += 1

    return count
# Windows present in a nuclear profile


def average_num_windows_per_np():
    count = sum_columns.sum()  # gets the full number of detections
    return count / num_windows()  # num windows per np


def smallest_num_windows_in_np():
    return sum_columns.min()


def largest_num_windows_in_np():
    return sum_columns.max()


def average_nps_per_window():
    count = sum_rows.sum()
    return count / count_npfs()  # detections of nps per window


def smallest_num_nps_in_window():
    return sum_rows.min()


def largest_num_nps_in_window():
    return sum_rows.max()


def radial_position_data():
    detections = sum_columns
    plt.scatter(detections.index, detections.values)
    plt.xlabel('Nuclear Profiles')
    plt.ylabel('Detections')
    plt.title('Frequency of window detection in each Nuclear Profile')
    plt.show()


def locus_volume_data():
    detections = sum_rows
    plt.scatter(detections.index, detections.values)
    plt.xlabel('Windows')
    plt.ylabel('Detections')
    plt.title('Frequency of nuclear profile detections in each window')
    plt.show()


def radial_position():
    sections = (largest_num_windows_in_np() - smallest_num_windows_in_np()) / 5
    section_1 = round(sections)
    section_2 = round(sections * 2)
    section_3 = round(sections * 3)
    section_4 = round(sections * 4)
    for value in sum_columns:
        if value <= section_1:
            print("Strongly apical")
        elif value <= section_2 and value > section_1:
            print("Somewhat apical")
        elif value <= section_3 and value > section_2:
            print("Neither apical nor equatorial")
        elif value <= section_4 and value > section_3:
            print("Somewhat equatorial")
        elif value > section_4:
            print("Strongly equatorial")


def compaction():
    sections = (largest_num_nps_in_window() -
                smallest_num_nps_in_window()) / 10
    section_1 = round(sections)
    section_2 = round(sections * 2)
    section_3 = round(sections * 3)
    section_4 = round(sections * 4)
    section_5 = round(sections * 5)
    section_6 = round(sections * 6)
    section_7 = round(sections * 7)
    section_8 = round(sections * 8)
    section_9 = round(sections * 9)
    for value in sum_rows:
        if value <= section_1:
            print("10")
        elif (value <= section_2 and value > section_1):
            print("9")
        elif (value <= section_3 and value > section_2):
            print("8")
        elif (value <= section_4 and value > section_3):
            print("7")
        elif (value <= section_5 and value > section_4):
            print("6")
        elif (value <= section_6 and value > section_5):
            print("5")
        elif (value <= section_7 and value > section_6):
            print("4")
        elif (value <= section_8 and value > section_7):
            print("3")
        elif (value < section_9 and value > section_8):
            print("2")
        else:
            print("1")

# Hist1 Functions


def hist1_num_windows():
    return len(nps_hist_1.iloc[0:])


def hist1_average_num_windows_per_np():
    count = hist1_sum_columns.sum()
    return count / hist1_num_windows()


def hist1_smallest_num_windows_in_np():
    return hist1_sum_columns.min()


def hist1_largest_num_windows_in_np():
    return hist1_sum_columns.max()


def hist1_average_nps_per_window():
    count = hist1_sum_rows.sum()
    return count / count_npfs()


def hist1_smallest_num_nps_in_window():
    return hist1_sum_rows.min()


def hist1_largest_num_nps_in_window():
    return hist1_sum_rows.max()


def hist1_radial_position():
    sections = (hist1_largest_num_windows_in_np() -
                hist1_smallest_num_windows_in_np()) / 5
    section_1 = round(sections)
    section_2 = round(sections * 2)
    section_3 = round(sections * 3)
    section_4 = round(sections * 4)
    for index, value in zip(hist1_sum_columns.index, hist1_sum_columns.values):
        if value <= section_1:
            print(index, ":Strongly apical")
        elif value <= section_2 and value > section_1:
            print(index, ":Somewhat apical")
        elif value <= section_3 and value > section_2:
            print(index, ":Neither apical nor equatorial")
        elif value <= section_4 and value > section_3:
            print(index, ":Somewhat equatorial")
        elif value > section_4:
            print(index, ":Strongly equatorial")


def hist1_compaction():
    sections = (largest_num_nps_in_window() -
                smallest_num_nps_in_window()) / 10
    section_1 = round(sections)
    section_2 = round(sections * 2)
    section_3 = round(sections * 3)
    section_4 = round(sections * 4)
    section_5 = round(sections * 5)
    section_6 = round(sections * 6)
    section_7 = round(sections * 7)
    section_8 = round(sections * 8)
    section_9 = round(sections * 9)
    temp_num = 0

    # Use only values from hist1 region, and compare to whole genome
    for index, value in zip(hist1_sum_rows.index, hist1_sum_rows.values):
        if value <= section_1:
            print(index, "10")
            temp_num += 10
        elif (value <= section_2 and value > section_1):
            print(index, "9")
            temp_num += 9
        elif (value <= section_3 and value > section_2):
            print(index, "8")
            temp_num += 8
        elif (value <= section_4 and value > section_3):
            print(index, "7")
            temp_num += 7
        elif (value <= section_5 and value > section_4):
            print(index, "6")
            temp_num += 6
        elif (value <= section_6 and value > section_5):
            print(index, "5")
            temp_num += 5
        elif (value <= section_7 and value > section_6):
            print(index, "4")
            temp_num += 4
        elif (value <= section_8 and value > section_7):
            print(index, "3")
            temp_num += 3
        elif (value < section_9 and value > section_8):
            print(index, "2")
            temp_num += 2
        else:
            print(index, "1")
            temp_num += 1

    average_compaction = temp_num / hist1_num_windows()

    print("Average compaction: ", average_compaction)

# Hist1 Jaccard index


def store_value():
    # Create a dictionary of arrays with a key value
    nps_arrays = {}
    for column in nps_hist_1.columns:
        # convert each column to its own array
        nps_array = nps_hist_1[column].to_numpy()
        # push back the newly formed array into the dictionary
        nps_arrays[column] = nps_array

    return nps_arrays


def hist1_jaccard_index(np1: np.ndarray, np2: np.ndarray):
    m11 = 0
    m01 = 0
    m10 = 0
    for index in range(len(np1)):
        if np1[index] == 1 and np2[index] == 1:
            m11 += 1
        elif np1[index] == 0 and np2[index] == 1:
            m01 += 1
        elif np1[index] == 1 and np2[index] == 0:
            m10 += 1

    jaccard_index = m11 / (m11 + m10 + m01)
    return jaccard_index


def store_jaccard_index(nps_arrays: dict):
    jacard_matrix = pd.DataFrame(
        index=nps_arrays.keys(), columns=nps_arrays.keys())
    for x in jacard_matrix.columns:
        for y in jacard_matrix.index:
            jacard_matrix.at[x, y] = hist1_jaccard_index(
                nps_arrays[x], nps_arrays[y])

    jacard_matrix = pd.read_csv('jaccard_matrix.csv')
    print(jacard_matrix)

    sns.heatmap(jacard_matrix)
    plt.show()


def hist1_distance_index(np1: np.ndarray, np2: np.ndarray):
    m11 = 0
    m01 = 0
    m10 = 0
    for index in range(len(np1)):
        if np1[index] == 1 and np2[index] == 1:
            m11 += 1
        elif np1[index] == 0 and np2[index] == 1:
            m01 += 1
        elif np1[index] == 1 and np2[index] == 0:
            m10 += 1

    jaccard_index = m11 / (m11 + m10 + m01)
    return 1 - jaccard_index


def store_distance_matrix(nps_arrays: dict):
    distance_matrix = pd.DataFrame(
        index=nps_arrays.keys(), columns=nps_arrays.keys())
    for x in distance_matrix.columns:
        for y in distance_matrix.index:
            distance_matrix.at[x, y] = hist1_distance_index(
                nps_arrays[x], nps_arrays[y])

    distance_matrix = pd.read_csv('distance_matrix.csv')
    print(distance_matrix)
    sns.heatmap(distance_matrix)
    plt.show()


def normalize_jaccard_index(np1: np.ndarray, np2: np.ndarray):
    m11 = 0
    m01 = 0
    m10 = 0
    for index in range(len(np1)):
        if np1[index] == 1 and np2[index] == 1:
            m11 += 1
        elif np1[index] == 0 and np2[index] == 1:
            m01 += 1
        elif np1[index] == 1 and np2[index] == 0:
            m10 += 1

    jaccard_index = m11 / min(m10 + m11, m01 + m11)
    return jaccard_index


def store_normalize_jaccard_matrix(nps_arrays: dict):
    normalize_jaccard_matrix = pd.DataFrame(
        index=nps_arrays.keys(), columns=nps_arrays.keys())
    for x in normalize_jaccard_matrix.columns:
        for y in normalize_jaccard_matrix.index:
            normalize_jaccard_matrix.at[x, y] = normalize_jaccard_index(
                nps_arrays[x], nps_arrays[y])

    normalize_jaccard_matrix = pd.read_csv('normalize_jaccard_matrix.csv')
    print(normalize_jaccard_matrix)
    sns.heatmap(normalize_jaccard_matrix)
    plt.show()


# For Heat Maps
# Stores each NP detections into a dictionary (Ex. [0,1,1,0,0,1,0], [0,1,0,0,0,0,1]) with the key being the specific NP
nps_arrays = store_value()


def select_k_value():
    normalize_df = pd.read_csv('normalize_jaccard_matrix.csv')
    normalize_df = normalize_df.sample(n=3, replace=True)
    print(normalize_df)
    c1 = []
    c2 = []
    c3 = []
    for x in normalize_df.columns:
        column = normalize_df[x]
        min_value = column.min()
        # If there are duplicates in the colum
        duplicates = column[column == min_value]
        if len(duplicates) > 1:
            min_index = np.random.choice(duplicates.index)
        else:
            min_index = column.idxmin()

        if min_index == normalize_df.index[0]:
            c1.append(x)

        if min_index == normalize_df.index[1]:
            c2.append(x)

        if min_index == normalize_df.index[2]:
            c3.append(x)

    np_relation = {
        normalize_df.index[0]: c1, normalize_df.index[1]: c2, normalize_df.index[2]: c3}
    pprint.pprint(np_relation)


"""
# Print the data
print("-----------------------------------------------------------------------")
print(data)
print("DATA for Activity 1")
print("-----------------------------------------------------------------------")
print("1.) Number of Genetic Windows:", num_windows())
print("2.) Number of Nuclear Profiles:", count_npfs())
print("3.) Average number of windows per nuclear profile:",
      average_num_windows_per_np())
print("4.) Smallest number of windows present in a NP:",
      smallest_num_windows_in_np())
print("    Largest number of windows present in a NP:",
      largest_num_windows_in_np())
print("5.) Average Number of nuclear profiles per window: ",
      average_nps_per_window())
print("    Smallest number nps in a window:", smallest_num_nps_in_window())
print("    Largest number nps in a window: ", largest_num_nps_in_window())
print("-----------------------------------------------------------------------")
print("-----------------------------------------------------------------------")
print("Data for Activity 2")
radial_position_data()
locus_volume_data()
print("-----------------------------------------------------------------------")
print("-----------------------------------------------------------------------")
print("Data for Activity 3")
print(hist1_df)
print("1.) Number of Genetic Windows:", hist1_num_windows())
print("2.) Number of Nuclear Profiles:", count_npfs())
print("3.) Average number of windows per nuclear profile:",
      hist1_average_num_windows_per_np())
print("4.) Smallest number of windows present in a NP:",
      hist1_smallest_num_windows_in_np())
print("    Largest number of windows present in a NP:",
      hist1_largest_num_windows_in_np())
print("5.) Average Number of nuclear profiles per window: ",
      hist1_average_nps_per_window())
print("    Smallest number nps in a window:",
      hist1_smallest_num_nps_in_window())
print("    Largest number nps in a window: ",
      hist1_largest_num_nps_in_window())
print("------------------------------------------------------------------------")
print("Radial position of Nuclear profiles in the Hist1 region:")
hist1_radial_position()
print("------------------------------------------------------------------------")
print("Compaction of Nuclear windows in the hist1 region:")
hist1_compaction()
print(" Data for Activity 4")
print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
nps_arrays = store_value()
print("Jaccard Matrix")
print("________________________________________________________________________")
store_jaccard_index(nps_arrays)
print("------------------------------------------------------------")
print("Difference Matrix")
store_distance_matrix(nps_arrays)

print("-----------------------------------------------------------------------_")
print("Activity 5: Clustering")
print("------------------------------------------------------------------------")
print("Normalized Jaccard Distance Matrix")
store_normalize_jaccard_matrix(nps_arrays)
"""
select_k_value()
quit()


# Is the HLB in different state? (Not formed, forming, formed)
# To determine this we must find difference and similarity between samples. We need to find three different states of these NPS
