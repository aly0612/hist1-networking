import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def read_data():
    linkage_matrix = pd.read_csv(
        "normalized_linkage.csv", delimiter=',')

    # convert labels to start at 0
    linkage_matrix.rename(columns=lambda x: "window " + str(int(x) - 69716),
                          index=lambda y: "window " + str(int(y) - 69716),
                          inplace=True)
    print("Linkage Matrix")
    print("-----------------------------")
    print(linkage_matrix)
    print("-----------------------------")

    return linkage_matrix


def read_feature_data():
    feature_df = pd.read_csv("features.csv", delimiter=',')
    isolate_data = feature_df[['name', 'Hist1', 'LAD']].copy()
    isolate_data.loc[:, 'name'] = [
        'Window {}'.format(i) for i in range(len(isolate_data))]
    print("Feature Data")
    print("-----------------------------")
    print(isolate_data)
    print("-----------------------------")
    return isolate_data


def calculate_centrality(linkage_matrix: pd.DataFrame):
    num_windows = len(linkage_matrix)
    degree_centrality = {}
    edge_matrix = contains_edge(linkage_matrix)

    for window in range(num_windows):
        count = edge_matrix.iloc[window].sum()
        degree_centrality[f"Window {window}"] = count / (num_windows - 1)

    return degree_centrality


def contains_edge(linkage_matrix: pd.DataFrame):
    edge_matrix = pd.DataFrame(
        index=linkage_matrix.index, columns=linkage_matrix.columns, dtype=int)
    average = calculate_l_average(linkage_matrix)

    for i in range(1, len(linkage_matrix)):
        for j in range(0, i):
            linkage_value = linkage_matrix.iloc[i, j]
            if linkage_value > average:
                edge_matrix.iloc[i, j] = 1
                edge_matrix.iloc[j, i] = 1
            else:
                edge_matrix.iloc[i, j] = 0
                edge_matrix.iloc[j, i] = 0

    # Set diagonal elements to 0
    for i in range(len(linkage_matrix)):
        edge_matrix.iloc[i, i] = 0

    return edge_matrix


def calculate_l_average(linkage_matrix: pd.DataFrame):
    linkage_values = []

    for i in range(0, len(linkage_matrix)):
        for j in range(0, i):

            linkage_values.insert(i, linkage_matrix.iloc[i, j])

    sum_values = sum(linkage_values)
    return sum_values / len(linkage_values)


def calculate_stats(degree_centrality: dict):

    # Calculate the average degree centrality
    sum_values = sum(degree_centrality.values())
    average = sum_values / len(degree_centrality)
    print("-----------------------------")
    print(f"Average degree centrality: {average}")
    print("-----------------------------")

    # caluclate max degree centrality
    max_value = max(degree_centrality.values())
    print(f"Max degree centrality: {max_value}")
    print("-----------------------------")

    # calculate min degree centrality
    min_value = min(degree_centrality.values())
    print(f"Min degree centrality: {min_value}")
    print("-----------------------------")

    # Sorted list of degree centrality
    sorted_degree_centrality = sorted(
        degree_centrality.items(), key=lambda x: x[1], reverse=False)

    degree_centrality_df = pd.DataFrame(sorted_degree_centrality, columns=[
                                        "Windows", "Degree Centrality"])

    print("Sorted list of degree centrality")
    print("-----------------------------")
    print(degree_centrality_df)
    print("-----------------------------")

    return degree_centrality_df


def find_communities(centrality_df: pd.DataFrame, edge_matrix: pd.DataFrame):
    # Get the top 5 hubs in descending order
    hubs = centrality_df.tail(5).iloc[::-1]
    print("Hubs")
    print("-----------------------------")
    print(hubs)
    print("-----------------------------")
    network = []

    for index, row in hubs.iterrows():
        hub_window = row["Windows"]
        hub_index = int(hub_window.split(" ")[1])
        # Extract the window number as an integer
        connected_windows = []

        for window, value in edge_matrix.iloc[hub_index].items():
            if value == 1:
                connected_windows.append(window)

        hub_size = len(connected_windows)
        community = {"Hub": hub_window, "Size": hub_size,
                     "Neighbors": connected_windows}

        network.append(community)
    print("Network")
    print("-----------------------------")
    for hubs in network:
        print(hubs)
        print("-----------------------------")

    return network


def feature_percentage(feature_df: pd.DataFrame, hubs: list):

    hub_window = hubs['Hub']
    count_Hist1 = 0
    count_LAD = 0

    for neighbor in hubs['Neighbors']:
        neighbor_index = int(neighbor.split(" ")[1])
        neighbor_features = feature_df.iloc[neighbor_index]
        if neighbor_features['Hist1'] == 1:
            count_Hist1 += 1
        if neighbor_features['LAD'] == 1:
            count_LAD += 1
    Hist1_Avg = count_Hist1 / hubs['Size']
    LAD_Avg = count_LAD / hubs['Size']

    community_features = {"Hub": hub_window,
                          'Hist1': Hist1_Avg, 'LAD': LAD_Avg}

    return community_features


def show_heatmap(linkage_matrix: pd.DataFrame, neighbors: list, position, title):
    # Create an empty matrix of zeros with the shape (81, 81)
    heatmap_data = pd.DataFrame(
        0, index=linkage_matrix.index, columns=linkage_matrix.columns)

    # Update the heatmap_data with the linkage values only for the community nodes
    for node1 in neighbors:
        for node2 in neighbors:
            if node1 != node2:
                heatmap_data.at[node1, node2] = linkage_matrix.at[node1, node2]
                # set negative values to 0
                if heatmap_data.at[node1, node2] < 0:
                    heatmap_data.at[node1, node2] = 0
    # Plot the heatmap
        # Plot the heatmap
    ax = plt.subplot(3, 2, position)
    sns.heatmap(heatmap_data, cmap="coolwarm", square=True, ax=ax)
    ax.set_title(title)


def main():
    linkage_matrix = read_data()
    calculate_centrality(linkage_matrix)
    centrality_df = calculate_stats(calculate_centrality(linkage_matrix))
    edge_matrix = contains_edge(linkage_matrix)
    my_network = find_communities(centrality_df, edge_matrix)
    feature_df = read_feature_data()
    print("Feature Percentages for Hist1 and LAD")
    print("-----------------------------")
    for hubs in my_network:
        print(feature_percentage(feature_df, hubs))
        print("-------------------------------")

    plt.figure(figsize=(20, 20))
    for idx, hub in enumerate(my_network):
        show_heatmap(linkage_matrix, hub['Neighbors'], idx + 1, hub['Hub'])

    plt.tight_layout()
    plt.show()

    sys.exit(0)


if __name__ == "__main__":
    main()
