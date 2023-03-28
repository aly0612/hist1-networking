import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def read_data():
    data = pd.read_csv(
        "GSE64881_segmentation_at_30000bp.passqc.multibam.txt", delimiter='\t')
    # Isolate Hist1 region
    hist1_data = data.iloc[69716:69797].copy()
    hist1_data.insert(0, 'Windows', hist1_data.apply(lambda row: ' '.join(
        [str(row['chrom']), str(row['start']), str(row['stop'])]), axis=1))
    hist1_data.drop(['chrom', 'start', 'stop'], axis=1, inplace=True)
    # Remove columns with all 0s
    for column in hist1_data:
        if hist1_data[column].sum() == 0:
            hist1_data.drop(column, axis=1, inplace=True)
    return hist1_data
# calculates Fa for each window

def calculate_frequency(hist1_data: pd.DataFrame):
    Np_freq_per_window = {}
    # get num NPs
    total_nps = hist1_data.iloc[:, 1:].columns.size
    count = 0
    # Store frequency of NPs in a dictionary for each window
    for index, row in hist1_data.iterrows():
        # This line calculates the sum of the detection of NPs in each window
        row_sum = row.iloc[1:].sum()
        frequency = row_sum / total_nps
        Np_freq_per_window[row.iloc[0]] = frequency

        count += 1

    return Np_freq_per_window
# Calculates Fab for each window

def calculate_coseg(hist1_data: pd.DataFrame):
    total_nps = hist1_data.iloc[:, 1:].columns.size
    coseg_values = {}
    # Store coseg values in a dictionary for each window
    for i in range(len(hist1_data)):
        for j in range(i + 1, len(hist1_data)):
            # Get the two windows that are being compared
            window1 = hist1_data.iloc[i]
            window2 = hist1_data.iloc[j]

            # This line calculates the sum of the common NPs in both windows that are being compared
            common_nps = (window1.iloc[1:] & window2.iloc[1:]).sum()
            coseg_value = common_nps / total_nps
            
            # Make the key for the dictionary
            window_pair_name = f"{window1.iloc[0]} and {window2.iloc[0]}"
            # Add the key and value to the dictionary
            coseg_values[window_pair_name] = coseg_value
    return coseg_values

def calculate_normalized_linkage(hist1_data: pd.DataFrame):
    # Calculate Fa and Fab
    fa_dict = calculate_frequency(hist1_data)
    fab_dict = calculate_coseg(hist1_data)
    normalized_linkage_df = pd.DataFrame(
        index=hist1_data.index, columns=hist1_data.index)
    
    for window_pair, fab in fab_dict.items():
        # Split the window pair name into two separate window names
        window1, window2 = window_pair.split(" and ")
        # Get the index of the window in the dataframe
        window1_index = hist1_data[hist1_data['Windows'] == window1].index[0]
        window2_index = hist1_data[hist1_data['Windows'] == window2].index[0]
        # Get the frequency of the window
        fa = fa_dict[window1]
        fb = fa_dict[window2]
        # Calculate the normalized linkage
        d_value = fab - (fa * fb)

        if d_value < 0:
            Dmax = min(fa * fb, (1 - fa) * (1 - fb))
        else:
            Dmax = min((1 - fa) * fb, fa * (1 - fb))

        if Dmax == 0:
            normalized_D = 0
        else:
            normalized_D = d_value / Dmax
        # Add the normalized linkage to the dataframe
        normalized_linkage_df.at[window1_index, window2_index] = normalized_D
        normalized_linkage_df.at[window2_index, window1_index] = normalized_D
        # Sets diagonal element to 1
        np.fill_diagonal(normalized_linkage_df.values, 1)
        # Converts the dataframe to float
        normalized_linkage_df = normalized_linkage_df.astype(float)

    print(normalized_linkage_df)
    return normalized_linkage_df

def generate_heatmap(normalized_linkage_df: pd.DataFrame):
    plt.figure(figsize=(20, 20))
    sns.heatmap(normalized_linkage_df, cmap="coolwarm", vmin=-1, vmax=1)
    plt.title("Normalized Linkage Heatmap", fontsize=20)
    plt.show()

def main():
    hist1_data = read_data()
    calculate_normalized_linkage(hist1_data)
    generate_heatmap(calculate_normalized_linkage(hist1_data))

    sys.exit(0)

if __name__ == "__main__":
    main()
