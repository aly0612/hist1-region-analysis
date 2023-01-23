import numpy as np
import pandas as pd

data = pd.read_csv('GSE64881_segmentation_at_30000bp.passqc.multibam.txt', delimiter="\t")
nps = data.iloc[:,4:]
windows = data.T.iloc[3:411]

def count_npfs():
  count = 0
  for column_name in data:
    if column_name[0] == 'F':
      count += 1
  
  return count

def average_windows():
  count = 0
  for column in data: #access each column
    if column[0] == 'F': #if nuclear profile
      for i in data[column].values: #loop the data within the column
        if i == 1: #if NP is present
          count += 1 #increment counter

  return count_npfs() / count     #divide total NPs by number of detections
                                  #nuclear profiles / detection in window

def smallest_windows():
 sum_obj_nps = nps.sum()
 return sum_obj_nps.min()

def largest_windows():
  sum_obj_nps = nps.sum()
  return sum_obj_nps.max()
  
def average_nps():
  sum_obj_windows = windows.sum()
  count = sum_obj_windows.sum()
  
  return count / 90876 #detections of nps per window

def smallest_nps():
  sum_obj_windows = windows.sum()
  return sum_obj_windows.min()

def largest_nps():
  sum_obj_windows = windows.sum()
  return sum_obj_windows.max()

print("DATA")
print(data)
print("-----------------------------------------------------------------------")
print("Smallest number of windows present:", smallest_windows())
print("1.) Number of Genetic Windows:", 90876) 
print("2.) Number of Nuclear Profiles:",count_npfs())
print("3.) Average number of windows on a Nuclear Profile", average_nps())
print("4.) Smallest number of windows present in a NP:", smallest_windows())
print("    Largest number of windows present in a NP: ", largest_windows())
print("5.) Average Number of Nps where a window is detected: ", average_windows())
print("    Smallest number nps in a window: ", smallest_nps())
print("    Largest number nps in a window:  ", largest_nps())
print("-----------------------------------------------------------------------")









 




  
  
  
  

  




