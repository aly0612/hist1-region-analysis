import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import csv
from tabulate import tabulate
from scipy.spatial.distance import pdist, squareform, cdist
from copy import deepcopy

original_original_data = pd.read_csv(
    'GSE64881_segmentation_at_30000bp.passqc.multibam.txt', delimiter="\t")
original_data = pd.read_csv('jaccard_matrix.csv')
normalize_data = pd.read_csv('normalize_jaccard_matrix.csv')

def select_initial_clusters(normalize_df: pd.DataFrame):
    # select 3 random samples in the dataframe, if the same samples replace
    normalize_df = normalize_df.sample(n=3, replace=True)
    c1 = []
    c2 = []
    c3 = []
    # iterate over columns in dataframe
    for x in normalize_df.columns:
        # retrieves the column with the name x
        column = normalize_df[x]
        # get the maximum value in the column x
        max_value = column.max()
        # If there are duplicates in the column, store the duplicates in a series
        duplicates = column[column == max_value]
        if len(duplicates) > 1:  # If there are duplicates
            # select a random index from the duplciates
            max_index = np.random.choice(duplicates.index)
            # assign the max index to a random
        else:
            # get the index of the max value
            max_index = column.idxmax()
        # assign the column to the cluster with the corresponding index label
        if max_index == normalize_df.index[0]:
            c1.append(x)

        if max_index == normalize_df.index[1]:
            c2.append(x)

        if max_index == normalize_df.index[2]:
            c3.append(x)
    # create the dictionary mapping each np to a specific cluster
    np_relation = {
        normalize_df.index[0]: c1, normalize_df.index[1]: c2, normalize_df.index[2]: c3}
    # create the table using tabulate
    table = tabulate(np_relation.items(), headers=['Cluster', 'NPs'])

# write the table to a file
    with open('clusters.csv', 'w', newline='') as csvfile:
      writer = csv.writer(csvfile)
      for line in table.split('\n'):
        writer.writerow(line.split('\t')) 
    
    return np_relation, normalize_df


clusters, cluster_df = select_initial_clusters(normalize_data)

print("Initial Clusters DF")
print("*******************")
print(cluster_df)
print("*******************")


# Find the medoid of each cluster
def find_medoid(cluster_df: pd.DataFrame, clusters: dict):
    medoids = []
    for np_key, np_values in clusters.items():
        cluster_np_values = cluster_df[np_values]
        distance_matrix = squareform(pdist(cluster_np_values, metric='jaccard'))
        # calculate the average dissimilarity between each NP and all other NPs in the cluster
        avg_dissimilarity = np.mean(distance_matrix, axis=1)
        # find the index of the NP with the minimal average dissimilarity
        medoid_index = np.argmin(avg_dissimilarity)
        # add the medoid to the list of medoids
        medoids.append(np_values[medoid_index])
    return medoids
medoids = find_medoid(cluster_df, clusters)
print("1st Medoids", medoids)
print("*******************")

def get_new_df(medoids: list, normalize_data: pd.DataFrame):
    nps = normalize_data.index.tolist()
   # new_nps = [np for np in nps if np not in medoids]
    new_df = normalize_data.loc[medoids, nps]
    return new_df

def assign_to_cluster(medoids: list, normalize_data: pd.DataFrame):
    new_df = get_new_df(medoids, normalize_data)
    c1 = []
    c2 = []
    c3 = []
    for x in new_df.columns:
        column = new_df[x]
        max_value = column.max()
        duplicates = column[column == max_value]
        if len(duplicates) > 1:
            max_index = np.random.choice(duplicates.index)
        else:
            max_index = column.idxmax()
        if max_index == new_df.index[0]:
            c1.append(x)

        if max_index == new_df.index[1]:
            c2.append(x)

        if max_index == new_df.index[2]:
            c3.append(x)
    np_relation = {
        new_df.index[0]: c1, new_df.index[1]: c2, new_df.index[2]: c3}
    table = tabulate(np_relation.items(), headers=['Cluster', 'NPs'])
    with open('cluster2.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for line in table.split('\n'):
            writer.writerow(line.split('\t'))
    return np_relation

medoids = find_medoid(get_new_df(medoids, normalize_data), assign_to_cluster(medoids, normalize_data))

def calculate_variance(medoids: list, normalize_data: pd.DataFrame):
    new_df = get_new_df(medoids, normalize_data)
    distance_matrix = squareform(pdist(new_df, metric='jaccard'))
    variance = np.sum(distance_matrix) / (len(medoids) * len(medoids))
    return variance



threshold = 0.005 # set the threshold for change in variance

prev_variance = 0  # initialize the previous variance
for i in range(50):
    medoids = find_medoid(get_new_df(medoids, normalize_data), assign_to_cluster(medoids, normalize_data))
    variance = calculate_variance(medoids, normalize_data)
    print("Variance", variance)
    print("Medoids", medoids)
    print("*******************")

    if abs(variance - prev_variance) < threshold:  # check if change in variance is less than threshold
        print("Variance change is less than threshold. Stopping iterations.")
        break

    prev_variance = variance  # update previous variance with current variance






def kmedoids_clustering(normalize_data, k=3, num_repeats=15):
    quality_metrics = []  # list to store the quality metric for each repetition
    best_medoids = None
    min_variance = float('inf')
    
    for i in range(num_repeats):
        # perform k-medoids clustering with random initial clusters
        clusters, cluster_df = select_initial_clusters(normalize_data)
        medoids = find_medoid(cluster_df, clusters)
        while True:
            new_medoids = find_medoid(cluster_df, clusters)
            if new_medoids == medoids:
                break
            else:
                medoids = new_medoids
                clusters = assign_to_cluster(medoids, normalize_data)
        
        # calculate the quality metric for this repetition
        variance = calculate_variance(medoids, normalize_data)
        quality_metrics.append(variance)
        
        # check if this set of clusters has the lowest variance so far
        if variance < min_variance and len(medoids) == k:
            min_variance = variance
            best_medoids = medoids
    
    return best_medoids, min_variance, quality_metrics


  
best_medoids, min_variance, quality_metrics = kmedoids_clustering(normalize_data)
print('Best medoids:', best_medoids)
print('Minimum variance:', min_variance)

medoid1 = best_medoids[0]
medoid2 = best_medoids[1]
medoid3 = best_medoids[2]
medoid1_df = original_original_data.loc[:,medoid1]
medoid2_df = original_original_data.loc[:,medoid2]
medoid3_df = original_original_data.loc[:,medoid3]
medoid1_df = medoid1_df.to_frame()
medoid2_df = medoid2_df.to_frame()
medoid3_df = medoid3_df.to_frame()
fig, axs = plt.subplots(1, 3, figsize=(12, 4))
sns.heatmap(medoid1_df, cmap='Blues', ax=axs[0])
sns.heatmap(medoid2_df, cmap='Blues', ax=axs[1])
sns.heatmap(medoid3_df, cmap='Blues', ax=axs[2])
plt.show()





