import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import csv
from tabulate import tabulate
from scipy.spatial.distance import pdist, squareform, cdist
from math import pi
original_original_data = pd.read_csv(
    'GSE64881_segmentation_at_30000bp.passqc.multibam.txt', delimiter="\t")
original_data = pd.read_csv('jaccard_matrix.csv')
normalize_data = pd.read_csv('normalize_jaccard_matrix.csv')

nps = original_original_data.iloc[:, 3:]
sum_columns = nps.sum()
sum_rows = nps.sum(axis=1)


def smallest_num_windows_in_np():
    return sum_columns.min()


def largest_num_windows_in_np():
    return sum_columns.max()


def radial_position():
    sections = (largest_num_windows_in_np() - smallest_num_windows_in_np()) / 5
    section_1 = round(sections)
    section_2 = round(sections * 2)
    section_3 = round(sections * 3)
    section_4 = round(sections * 4)

    # initialize an empty dictionary to store the radial positions of NPs
    radial_positions = {}

    for np_name, value in sum_columns.items():
        if value <= section_1:
            radial_positions[np_name] = "StrA"
        elif value <= section_2 and value > section_1:
            radial_positions[np_name] = "SomeA"
        elif value <= section_3 and value > section_2:
            radial_positions[np_name] = "Neither"
        elif value <= section_4 and value > section_3:
            radial_positions[np_name] = "SomeE"
        elif value > section_4:
            radial_positions[np_name] = "StrE"

    return radial_positions


radial_position_data = radial_position()


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

# Find the medoid of each cluster


def find_medoid(cluster_df: pd.DataFrame, clusters: dict):
    medoids = []
    for np_key, np_values in clusters.items():
        cluster_np_values = cluster_df[np_values]
        distance_matrix = squareform(
            pdist(cluster_np_values, metric='jaccard'))
        # calculate the average dissimilarity between each NP and all other NPs in the cluster
        avg_dissimilarity = np.mean(distance_matrix, axis=1)
        # find the index of the NP with the minimal average dissimilarity
        medoid_index = np.argmin(avg_dissimilarity)
        # add the medoid to the list of medoids
        medoids.append(np_values[medoid_index])
    return medoids


medoids = find_medoid(cluster_df, clusters)


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


medoids = find_medoid(get_new_df(medoids, normalize_data),
                      assign_to_cluster(medoids, normalize_data))


def calculate_variance(medoids: list, normalize_data: pd.DataFrame):
    new_df = get_new_df(medoids, normalize_data)
    distance_matrix = squareform(pdist(new_df, metric='jaccard'))
    variance = np.sum(distance_matrix) / (len(medoids) * len(medoids))
    return variance


threshold = 0.005  # set the threshold for change in variance
prev_variance = 0  # initialize the previous variance

for i in range(50):
    medoids = find_medoid(get_new_df(medoids, normalize_data),
                          assign_to_cluster(medoids, normalize_data))
    variance = calculate_variance(medoids, normalize_data)

    # check if change in variance is less than threshold
    if abs(variance - prev_variance) < threshold:
        break

    prev_variance = variance  # update previous variance with current variance


def kmedoids_clustering(normalize_data, k=3, num_repeats=50):
    quality_metrics = []  # list to store the quality metric for each repetition
    best_medoids = None
    min_variance = float('inf')
    best_clusters = None
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
            best_clusters = assign_to_cluster(best_medoids, normalize_data)

    return best_medoids, min_variance, quality_metrics, best_clusters


best_medoids, min_variance, quality_metrics, best_clusters = kmedoids_clustering(
    normalize_data)
print("Feature selection - Activity1")
print("*******************")
print("Optimizing Clusters")
print("*******************")
print("Best Medoids", best_medoids)
print("Min Variance", min_variance)

c1 = []
c2 = []
c3 = []
#creating a list of the best clusters, this will allow us to create a new dataframe with the best clusters
#the best clusters are the ones that have the best medoids
for key in best_clusters:
    if key == best_medoids[0]:
        c1 = best_clusters[key]
    if key == best_medoids[1]:
        c2 = best_clusters[key]
    if key == best_medoids[2]:
        c3 = best_clusters[key]
#creating dataframes for each cluster based on the best medoids
#A 1 will indicate that the NP is detected in the window within that specific cluster
cluster1_data = original_original_data.loc[:, c1].T
cluster2_data = original_original_data.loc[:, c2].T
cluster3_data = original_original_data.loc[:, c3].T

print("*******************")
print("Creating Heatmap")

# create grid of plots
fig, axs = plt.subplots(1, 3, figsize=(12, 4))

# create heatmap for each cluster and add to grid of plots
sns.heatmap(cluster1_data,   ax=axs[0])
sns.heatmap(cluster2_data,  ax=axs[1])
sns.heatmap(cluster3_data,   ax=axs[2])

# set titles for each heatmap

axs[0].set_title('Cluster 1')
axs[1].set_title('Cluster 2')
axs[2].set_title('Cluster 3')

# adjust spacing between plots
# this will allows all plots to be merged into one image
plt.subplots_adjust(wspace=0.3) 

# display the plot
plt.show()

print("*******************")
print("Feature Selection - Activity 2")

feature_df = pd.read_csv('features.csv')
print("*******************")
print("All Features DataFrame")
print(feature_df)
print("*******************")

def isolate_features(feature_df: pd.DataFrame):
    isolated_features = feature_df.loc[:, ['Hist1', 'LAD']]
    return isolated_features

print("Isolated Features DataFrame - Hist1 and LAD")
print(isolate_features(feature_df))


def isolate_all_features(feature_df: pd.DataFrame):
    isolated_features = feature_df.loc[:, ['Hist1', 'Vmn', 'LAD', 'RNAPII-S2P', 'RNAPII-S5P', 'RNAPII-S7P',
                                           'Enhancer', 'H3K9me3', 'H3K20me3', 'h3k27me3', 'H3K36me3', 'NANOG', 'pou5f1', 'sox2', 'CTCF-7BWU']]
    return isolated_features


def isolate_hist1_region(cluster_df: pd.DataFrame):
    cluster_df = cluster_df.iloc[69716:69797].copy()
    return cluster_df
print("*******************")
print("Isolated Hist1 Region DataFrame - Cluster 1")
print(isolate_hist1_region(cluster1_data.T))
print("*******************")
print("Isolated Hist1 Region DataFrame - Cluster 2")
print(isolate_hist1_region(cluster2_data.T))
print("*******************")
print("Isolated Hist1 Region DataFrame - Cluster 3")
print(isolate_hist1_region(cluster3_data.T))
print("*******************")


isolate_hist1_region(cluster1_data.T).to_csv('cluster1.csv', index=False)
isolate_hist1_region(cluster2_data.T).to_csv('cluster2.csv', index=False)
isolate_hist1_region(cluster3_data.T).to_csv('cluster3.csv', index=False)


print(isolate_all_features(feature_df))


def calculate_np_feature_percentage(cluster_df, feature_df):
    cluster_df = isolate_hist1_region(cluster_df)
    cluster_df.index = range(len(cluster_df.index))
    feature_df = isolate_features(feature_df)
    np_feature_counts = {}
    total_windows = len(cluster_df)

    for np in cluster_df.columns:
        np_feature_counts[np] = {}
        for feature in feature_df.columns:
            count = 0
            for index, row in cluster_df.iterrows():
                if row[np] == 1 and feature_df.loc[index, feature] >= 1:
                    count += 1

            np_feature_counts[np][feature] = count / total_windows

    percentages_df = pd.DataFrame(np_feature_counts)
    return percentages_df


def calculate_np_feature_percentage__all(cluster_df, feature_df):
    cluster_df = isolate_hist1_region(cluster_df)
    cluster_df.index = range(len(cluster_df.index))
    feature_df = isolate_all_features(feature_df)
    np_feature_counts = {}
    total_windows = len(cluster_df)

    for np in cluster_df.columns:
        np_feature_counts[np] = {}
        for feature in feature_df.columns:
            count = 0
            for index, row in cluster_df.iterrows():
                if row[np] == 1 and feature_df.loc[index, feature] >= 1:
                    count += 1

            np_feature_counts[np][feature] = count / total_windows

    percentages_df = pd.DataFrame(np_feature_counts)
    return percentages_df


cluster1_feature_percent = calculate_np_feature_percentage(
    cluster1_data.T, feature_df)

cluster1_feature_percent_all = calculate_np_feature_percentage__all(
    cluster1_data.T, feature_df)
cluster2_feature_percent_all = calculate_np_feature_percentage__all(
    cluster2_data.T, feature_df)
cluster3_feature_percent_all = calculate_np_feature_percentage__all(
    cluster3_data.T, feature_df)


def create_his1_box_plot():
    cluster1_percentages_df = calculate_np_feature_percentage(
        cluster1_data.T, feature_df)
    cluster2_percentages_df = calculate_np_feature_percentage(
        cluster2_data.T, feature_df)
    cluster3_percentages_df = calculate_np_feature_percentage(
        cluster3_data.T, feature_df)

    # combine the percentage data for each cluster into a single dataframe
    all_cluster_percentages = pd.concat([cluster1_percentages_df.loc['Hist1'],
                                         cluster2_percentages_df.loc['Hist1'],
                                         cluster3_percentages_df.loc['Hist1']],
                                        axis=1)
    all_cluster_percentages.columns = ['Cluster1', 'Cluster2', 'Cluster3']

    # create the box plot using catplot
    sns.set(style='whitegrid')
    sns.catplot(data=all_cluster_percentages, kind='box')
    plt.xlabel('Cluster')
    plt.ylabel('Percentage of windows with Hist1 feature')
    plt.title('Hist1')
    plt.show()


def create_LAD_box_plot():
    cluster1_percentages_df = calculate_np_feature_percentage(
        cluster1_data.T, feature_df)
    cluster2_percentages_df = calculate_np_feature_percentage(
        cluster2_data.T, feature_df)
    cluster3_percentages_df = calculate_np_feature_percentage(
        cluster3_data.T, feature_df)

    # combine the percentage data for each cluster into a single dataframe
    all_cluster_percentages = pd.concat([cluster1_percentages_df.loc['LAD'],
                                         cluster2_percentages_df.loc['LAD'],
                                         cluster3_percentages_df.loc['LAD']],
                                        axis=1)
    all_cluster_percentages.columns = ['Cluster1', 'Cluster2', 'Cluster3']

    # create the box plot using catplot
    sns.set(style='whitegrid')
    sns.catplot(data=all_cluster_percentages, kind='box')
    plt.xlabel('Cluster')
    plt.ylabel('Percentage of windows with LAD feature')
    plt.title('LAD')
    plt.show()

print("*******************")
print("Creating box plots for Hist1 and LAD features")
create_his1_box_plot()
create_LAD_box_plot()
print("*******************")

cluster1_Nps = cluster1_data.T.columns.to_list()
cluster2_Nps = cluster2_data.T.columns.to_list()
cluster3_Nps = cluster3_data.T.columns.to_list()


def calculate_radial_position_percentages(cluster, radial_positions: dict):
    # Initialize a dictionary to count the number of NPs that fall into each radial position category
    radial_position_counts = {"StrA": 0, "SomeA": 0,
                              "Neither": 0, "SomeE": 0, "StrE": 0}

    # Iterate through the NP names in the cluster and count the number of NPs that fall into each radial position category
    for np_name in cluster:
        radial_position = radial_positions.get(np_name)
        if radial_position is not None:  # Make sure the NP has a radial position assigned
            radial_position_counts[radial_position] += 1

    # Calculate the percentage of NPs that fall into each radial position category
    total_nps = sum(radial_position_counts.values())
    radial_position_percentages = {
        k: v/total_nps for k, v in radial_position_counts.items()}

    return radial_position_percentages


c1_radial_percentages = calculate_radial_position_percentages(
    cluster1_Nps, radial_position_data)
c2_radial_percentages = calculate_radial_position_percentages(
    cluster2_Nps, radial_position_data)
c3_radial_percentages = calculate_radial_position_percentages(
    cluster3_Nps, radial_position_data)


def generate_bargraph(radial_percentges: dict):
    # Create a bar graph of the radial position percentages
    plt.bar(radial_percentges.keys(), radial_percentges.values())
    plt.ylabel('Percentage of NPs')
    plt.show()

print("*******************")
print("Feature Selection - Activity 3")
print("Creating bar graphs for radial position percentages")
generate_bargraph(c1_radial_percentages)
generate_bargraph(c2_radial_percentages)
generate_bargraph(c3_radial_percentages)


def calculate_average_percentage(cluster1_feature_percentages: pd.DataFrame):
    average_percentages = {}
    for features in cluster1_feature_percentages.index:
        average_percentages[features] = cluster1_feature_percentages.loc[features].mean(
        ) * 100
    average_percentages_df = pd.DataFrame(average_percentages, index=[0])
    return average_percentages_df

cluster1_average_percentages = calculate_average_percentage(
    cluster1_feature_percent_all)
cluster2_average_percentages = calculate_average_percentage(
    cluster2_feature_percent_all)
cluster3_average_percentages = calculate_average_percentage(
    cluster3_feature_percent_all)

all_cluster_feature_percentages = pd.concat(
    [cluster1_average_percentages, cluster2_average_percentages, cluster3_average_percentages], axis=0)
all_cluster_feature_percentages.index = ['Cluster1', 'Cluster2', 'Cluster3']

print("*******************")
print("Feature Selection - Activity 4")
print("Average percentage of windows with each feature for each cluster")
print(all_cluster_feature_percentages)


def create_radar_chart(all_cluster_feature_percentages: pd.DataFrame):
    # number of variable
    categories = list(all_cluster_feature_percentages)
    N = len(categories)
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    # loop through the rows of the data frame and plot each one as a separate data point
    for i, row in all_cluster_feature_percentages.iterrows():
        values = row.values.flatten().tolist()
        values += values[:1]
        angles = [n / float(N) * 2 * pi for n in range(N)]
        angles += angles[:1]
        ax.plot(angles, values, linewidth=1,
                linestyle='solid', label=f'Cluster {i}')
        ax.fill(angles, values, alpha=0.1)

    plt.xticks(angles[:-1], categories, color='grey', size=8)

    ax.set_rlabel_position(0)
    plt.yticks([5, 10, 15, 20], ["5", "10", "15", "20"], color="grey", size=7)
    plt.ylim(0, 20)

    ax.plot(angles, values, linewidth=1, linestyle='solid')

    ax.fill(angles, values, 'b', alpha=0.1)

    plt.show()

print("*******************")
print("Radar chart of average percentage of windows with each feature for each cluster")
create_radar_chart(all_cluster_feature_percentages)
