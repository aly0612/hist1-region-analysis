import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv(
    'GSE64881_segmentation_at_30000bp.passqc.multibam.txt', delimiter="\t")
nps = data.iloc[:, 3:] #Read only the nps and columns, eliminates unwanted data
sum_columns = nps.sum()
sum_rows = nps.sum(axis=1)

#Isolation of Hist1 region
hist1_df = data.iloc[69714:69795].copy()
nps_hist_1 = hist1_df.iloc[:,3:]
hist1_sum_columns = nps_hist_1.sum()
hist1_sum_rows = nps_hist_1.sum(axis=1)

def num_windows():
    return len(nps.iloc[0:])
#count the nuclear profiles, all start with F
def count_npfs():
    count = 0
    for column_name in data:
        if column_name[0] == 'F':
            count += 1

    return count
#Windows present in a nuclear profile
def average_num_windows_per_np():
    count = sum_columns.sum()
    return count / num_windows()

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
    sections = (largest_num_nps_in_window() - smallest_num_nps_in_window()) / 10
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
        if value <= section_1 :
            print ("10")
        elif (value <= section_2 and value > section_1):
            print ("9")
        elif (value <= section_3 and value > section_2):
            print ("8")
        elif(value <= section_4 and value > section_3):
            print ("7")
        elif (value <= section_5 and value > section_4):
            print("6")
        elif (value <= section_6 and value > section_5):
            print("5")
        elif(value <= section_7 and value > section_6):
            print("4")
        elif(value <= section_8 and value > section_7):
            print("3")
        elif(value < section_9 and value > section_8):
            print("2")
        else:
            print("1")

#Hist1 Functions

def hist1_num_windows():
    return len(nps_hist_1.iloc[0:])

def hist1_average_num_windows_per_np():
    count = hist1_sum_columns.sum()
    return count / hist1_num_windows()

def hist1_smallest_num_windows_in_np():
    return hist1_sum_columns.min()

def hist1_largest_num_windows_in_np():
    return hist1_sum_columns.max()

def  hist1_average_nps_per_window():
    count = hist1_sum_rows.sum()
    return count / count_npfs()

def hist1_smallest_num_nps_in_window():
    return hist1_sum_rows.min()

def hist1_largest_num_nps_in_window():
    return hist1_sum_rows.max()

def hist1_radial_position():
    sections = (hist1_largest_num_windows_in_np() - hist1_smallest_num_windows_in_np()) / 5
    section_1 = round(sections)
    section_2 = round(sections * 2)
    section_3 = round(sections * 3)
    section_4 = round(sections * 4)
    for index, value in zip(hist1_sum_columns.index, hist1_sum_columns.values):
        if value <= section_1:
            print(index,":Strongly apical")
        elif value <= section_2 and value > section_1:
            print(index,":Somewhat apical")
        elif value <= section_3 and value > section_2:
            print(index,":Neither apical nor equatorial")
        elif value <= section_4 and value > section_3:
            print(index,":Somewhat equatorial")
        elif value > section_4:
            print(index, ":Strongly equatorial")

def hist1_compaction():
    sections = (largest_num_nps_in_window() - smallest_num_nps_in_window()) / 10
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

    for index, value in zip(hist1_sum_rows.index, hist1_sum_rows.values): #Use only values from hist1 region, and compare to whole genome
        if value <= section_1 :
            print (index, "10")
            temp_num += 10
        elif (value <= section_2 and value > section_1):
            print (index, "9")
            temp_num += 9
        elif (value <= section_3 and value > section_2):
            print (index, "8")
            temp_num += 8
        elif(value <= section_4 and value > section_3):
            print (index, "7")
            temp_num += 7
        elif (value <= section_5 and value > section_4):
            print(index, "6")
            temp_num += 6
        elif (value <= section_6 and value > section_5):
            print(index, "5")
            temp_num += 5
        elif(value <= section_7 and value > section_6):
            print(index, "4")
            temp_num += 4
        elif(value <= section_8 and value > section_7):
            print(index, "3")
            temp_num += 3
        elif(value < section_9 and value > section_8):
            print(index, "2")
            temp_num += 2
        else:
            print(index, "1")
            temp_num += 1
    
    average_compaction = temp_num / hist1_num_windows()

    print("Average compaction: ", average_compaction)

    print(hist1_sum_rows)


print(data)
print("DATA for Activity 1")
print("-----------------------------------------------------------------------")
print("1.) Number of Genetic Windows:", num_windows())
print("2.) Number of Nuclear Profiles:", count_npfs())
print("3.) Average number of windows per nuclear profile:", average_num_windows_per_np())
print("4.) Smallest number of windows present in a NP:", smallest_num_windows_in_np())
print("    Largest number of windows present in a NP:", largest_num_windows_in_np())
print("5.) Average Number of nuclear profiles per window: ", average_nps_per_window())
print("    Smallest number nps in a window:", smallest_num_nps_in_window())
print("    Largest number nps in a window: ", largest_num_nps_in_window())
print("-----------------------------------------------------------------------")
print("-----------------------------------------------------------------------")
print("Data for Activity 2")
radial_position_data()
locus_volume_data()
#radial_position()
#compaction()
print("-----------------------------------------------------------------------")
print("-----------------------------------------------------------------------")
print("Data for Activity 3")
print(hist1_df)
print("1.) Number of Genetic Windows:", hist1_num_windows())
print("2.) Number of Nuclear Profiles:", count_npfs())
print("3.) Average number of windows per nuclear profile:", hist1_average_num_windows_per_np())
print("4.) Smallest number of windows present in a NP:", hist1_smallest_num_windows_in_np())
print("    Largest number of windows present in a NP:", hist1_largest_num_windows_in_np())
print("5.) Average Number of nuclear profiles per window: ", hist1_average_nps_per_window())
print("    Smallest number nps in a window:", hist1_smallest_num_nps_in_window())
print("    Largest number nps in a window: ", hist1_largest_num_nps_in_window())
print("------------------------------------------------------------------------")
print("Radial position of Nuclear profiles in the Hist1 region:")
hist1_radial_position()
print("------------------------------------------------------------------------")
print("Compaction of Nuclear windows in the hist1 region:")
hist1_compaction()
