import sys
import numpy as np
import pandas as pd


def read_data():
    linkage_matrix = pd.read_csv(
        "normalized_linkage.csv", delimiter=',')
    print(linkage_matrix)
    return linkage_matrix

def calculate_centrality(linkage_matrix: pd.DataFrame):
    vertices = contains_vertex(linkage_matrix)
    degree_centrality = {}
    num_windows = len(linkage_matrix)
    for window in range(num_windows):
        count = 0
        for key, value in vertices.items():
            if value == 1 and key.startswith(f"Window {window}"):
                    count += 1
        
        degree_centrality[f"Window {window}"] = count / (num_windows - 1)
    
    

    calculate_stats(degree_centrality)
    return degree_centrality

def contains_vertex(linkage_matrix: pd.DataFrame):
    vertices = {}
    average = calculate_l_average(linkage_matrix)
    for i in range(1, len(linkage_matrix)):
        for j in range(0, i):
            linkage_value = linkage_matrix.iloc[i, j]
            key = f"Window {i} and Window {j}"
            if linkage_value > average:
                vertices[key] = 1

            else:
                vertices[key] = 0
      


    return vertices

def calculate_l_average(linkage_matrix: pd.DataFrame):
    linkage_values = []

    for i in range(1, len(linkage_matrix)):
        for j in range(0, i):

            linkage_values.insert(i, linkage_matrix.iloc[i, j])

    sum_values = sum(linkage_values)
    return sum_values / len(linkage_values)

def calculate_stats(degree_centrality: dict):
    # Calculate the average degree centrality
    sum_values = sum(degree_centrality.values())
    average = sum_values / len(degree_centrality)
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
    
    
    #Show full dataframe
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)
   
    degree_centrality_df = pd.DataFrame(sorted_degree_centrality, columns=["Windows", "Degree Centrality"])
    
    print("Sorted list of degree centrality")
    print("-----------------------------")
    print(degree_centrality_df)

def main():
    linkage_matrix = read_data()
    calculate_centrality(linkage_matrix)
    sys.exit(0)


if __name__ == "__main__":
    main()
