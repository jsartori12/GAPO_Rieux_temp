#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 12:55:25 2024

@author: joao.sartori
"""

import pandas as pd
import matplotlib.pyplot as plt
import Funcs_output_ana
import modeling_functions_pyrosetta
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO


path_outputs = "Results/Naturals"

#apt_esm_esm.txt  apt_esm_random.txt  apt_rosetta_esm.txt  apt_rosetta_random.txt
#####
apt_esm_esm= Funcs_output_ana.Read_data(f"{path_outputs}/apt_esm_esm.txt", 1)
unique_apt_esm_esm= apt_esm_esm.drop_duplicates(subset=['Sequence'])
#####
apt_esm_random = Funcs_output_ana.Read_data(f"{path_outputs}/apt_esm_random.txt", 2)
unique_apt_esm_random = apt_esm_random.drop_duplicates(subset=['Sequence'])
#####

apt_rosetta_esm= Funcs_output_ana.Read_data(f"{path_outputs}/apt_rosetta_esm.txt", 1)
unique_apt_rosetta_esm= apt_rosetta_esm.drop_duplicates(subset=['Sequence'])

apt_rosetta_random= Funcs_output_ana.Read_data(f"{path_outputs}/apt_rosetta_random.txt", 2)
unique_apt_rosetta_random= apt_rosetta_random.drop_duplicates(subset=['Sequence'])

min_dG_row = unique_apt_rosetta_random.loc[unique_apt_rosetta_random['dG'].idxmin()]

# Get the sequence associated with the lowest dG value
# sequence_with_lowest_dG = min_dG_row

# print(sequence_with_lowest_dG)
apt_esm_esm_best, apt_esm_random_best = Funcs_output_ana.List_min_dgs([apt_esm_esm, apt_esm_random],["apt_esm_esm", "apt_esm_random"], 
                                                                      apt_function = "dG", min_or_max=max, out = "esm")

apt_rosetta_esm_best, apt_rosetta_random_best = Funcs_output_ana.List_min_dgs([apt_rosetta_esm, apt_rosetta_random],["apt_rosetta_esm", "apt_rosetta_random"], 
                                                                              apt_function = "dG", out = "rosetta")



unique_apt_esm_esm_algn = Funcs_output_ana.calculate_identity_matrix_distance2(df=unique_apt_esm_esm, out = "esm_esm")
unique_apt_esm_esm_algn = pd.DataFrame(unique_apt_esm_esm_algn)

unique_apt_esm_random_algn = Funcs_output_ana.calculate_identity_matrix_distance2(df=unique_apt_esm_random, out = "esm_random")
unique_apt_esm_random_algn = pd.DataFrame(unique_apt_esm_random_algn)


unique_apt_rosetta_esm_algn = Funcs_output_ana.calculate_identity_matrix_distance2(df=unique_apt_rosetta_esm, out = "rosetta_esm")
unique_apt_rosetta_esm_algn = pd.DataFrame(unique_apt_rosetta_esm_algn)

unique_apt_rosetta_random_algn = Funcs_output_ana.calculate_identity_matrix_distance2(df=unique_apt_rosetta_random, out = "rosetta_random")
unique_apt_rosetta_random_algn = pd.DataFrame(unique_apt_rosetta_random_algn)

# unique_apt_esm_esm_algn = pd.read_csv("esm_esm_distance_matrix.csv")
# unique_apt_esm_random_algn = pd.read_csv("esm_random_distance_matrix.csv")

# unique_apt_rosetta_esm_algn = pd.read_csv("rosetta_esm_distance_matrix.csv")
# unique_apt_rosetta_random_algn = pd.read_csv("rosetta_random_distance_matrix.csv")

# unique_apt_esm_esm_algn = Funcs_output_ana.Create_align_fasta(df = unique_apt_esm_esm, out = "esm_esm")
# unique_apt_esm_random_algn = Funcs_output_ana.Create_align_fasta(df = unique_apt_esm_random, out = "esm_random")

# unique_apt_rosetta_esm_algn = Funcs_output_ana.Create_align_fasta(df = unique_apt_rosetta_esm, out = "rosetta_esm")
# unique_apt_rosetta_random_algn = Funcs_output_ana.Create_align_fasta(df = unique_apt_rosetta_random, out = "rosetta_random")


# teste_dist = Funcs_output_ana.Diss_dist(unique_apt_esm_esm_algn)

# # Define the integer value
# int_value = -11

# # Create a Pandas Series with values ranging from 0 to 100
# series = pd.Series(range(101))

# # Replace the values with the given integer
# series = pd.Series([int_value] * len(series))


#testegif = Funcs_output_ana.Plot_bests(df1= rep1, df2=rep2, df3=rep3, natural=series, out="teste")

# all_sequences_esm = pd.concat([apt_esm_esm,apt_esm_random])
# all_sequences_rosetta = pd.concat([apt_rosetta_esm,apt_rosetta_random])

# len(apt_esm_esm["Sequence"][0])

# CDRs = list(range(0,len(apt_esm_esm["Sequence"][0])))

# esm_esm_vhse = Funcs_output_ana.vhse_zcales_converter(df=apt_esm_esm, CDR_positions = CDRs, descriptor = "vhse")
# esm_random_vhse = Funcs_output_ana.vhse_zcales_converter(df=apt_esm_random, CDR_positions = CDRs, descriptor = "vhse")

# rosetta_esm_vhse = Funcs_output_ana.vhse_zcales_converter(df=apt_rosetta_esm, CDR_positions = CDRs, descriptor = "vhse")
# rosetta_random_vhse = Funcs_output_ana.vhse_zcales_converter(df=apt_rosetta_random, CDR_positions = CDRs, descriptor = "vhse")


# all_vhse_esm = pd.concat([esm_esm_vhse, esm_random_vhse])
# all_vhse_rosetta = pd.concat([rosetta_esm_vhse, rosetta_random_vhse])


# all_conds_cycleesm = pd.concat([unique_apt_esm_esm["Cycle"], unique_apt_esm_random["Cycle"]])
# all_conds_repesm = pd.concat([unique_apt_esm_esm["Condition"], unique_apt_esm_random["Condition"]])

# all_conds_cyclerose = pd.concat([unique_apt_rosetta_esm["Cycle"], unique_apt_rosetta_random["Cycle"]])
# all_conds_reprose = pd.concat([unique_apt_rosetta_esm["Condition"], unique_apt_rosetta_random["Condition"]])


# all_conds_cycle_dfesm = pd.DataFrame(all_conds_cycleesm)
# all_conds_rep_dfesm = pd.DataFrame(all_conds_repesm)

# all_conds_cycle_dfrosetta = pd.DataFrame(all_conds_cyclerose)
# all_conds_rep_dfrosetta = pd.DataFrame(all_conds_reprose)


esm_esm_dist = pd.read_csv("esm_esm_distance_matrix.csv")
esm_random_dist = pd.read_csv("esm_esm_distance_matrix.csv")

rosetta_esm_dist = pd.read_csv("esm_esm_distance_matrix.csv")
rosetta_random_dist = pd.read_csv("esm_esm_distance_matrix.csv")

##descomentar
tsne_esm_esm = Funcs_output_ana.red_tsne(data = esm_esm_dist, dimension=2, prp=100, label_list = unique_apt_esm_esm["Cycle"].tolist(), out = "esm_esm")
tsne_esm_random = Funcs_output_ana.red_tsne(data = esm_esm_dist, dimension=2, prp=100, label_list = unique_apt_esm_esm["Cycle"].tolist(), out = "esm_random")

tsne_rosetta_esm = Funcs_output_ana.red_tsne(data = esm_esm_dist, dimension=2, prp=100, label_list = unique_apt_esm_esm["Cycle"].tolist(), out = "rosetta_esm")
tsne_rosetta_random = Funcs_output_ana.red_tsne(data = esm_esm_dist, dimension=2, prp=100, label_list = unique_apt_esm_esm["Cycle"].tolist(), out = "rosetta_random")

#tsne_esm_esm = pd.read_csv("esm_esm_tsne.csv")
#tsne_esm_random = pd.read_csv("esm_random_tsne.csv")

#tsne_rosetta_esm =  pd.read_csv("rosetta_esm_tsne.csv")
#tsne_rosetta_random = pd.read_csv("rosetta_random_tsne.csv")



def treat_tsne_data(tsne_csv, df, out,dG_min=0.80, dG_max=0.92):
    # Read the t-SNE CSV file
    tsne = pd.read_csv(tsne_csv)
    
    # Drop the first column (assuming it's an index column)
    tsne = tsne.iloc[:, 1:]
    
    # Rename the column 'Conditions' to 'Population'
    tsne = tsne.rename(columns={'Conditions': 'Population'})
    
    # Ensure both dataframes have matching rows
    if len(df) >= len(tsne):
        # If df is larger or equal, trim df to match tsne's number of rows
        df_trimmed = df.iloc[:len(tsne), :]
    else:
        # If df is smaller, trim tsne to match df's number of rows
        tsne = tsne.iloc[:len(df), :]
        df_trimmed = df

    # Assign the corresponding dG values from df_trimmed to the tsne dataframe
    tsne["dG"] = df_trimmed["dG"].tolist()
    
    # Call the plot_panel_3d function to generate the plot
    print(tsne)
    Funcs_output_ana.plot_panel_3d(tsne, dG_min=dG_min, dG_max=dG_max, out = out) 

treat_tsne_data("esm_esm_tsne.csv", unique_apt_esm_esm, dG_min= 0.70, dG_max= 0.95,out = "esm_esm")
treat_tsne_data("esm_random_tsne.csv", unique_apt_esm_random, dG_min= 0.70, dG_max= 0.95,out = "esm_random")
treat_tsne_data("rosetta_esm_tsne.csv", unique_apt_rosetta_esm, dG_min= -450, dG_max=-200, out = "rosetta_esm")
treat_tsne_data("rosetta_random_tsne.csv", unique_apt_rosetta_random, dG_min= -450, dG_max=-200, out = "rosetta_random")


#unique_apt_rosetta_esm["dG"].max()




















# tsne_all = tsne_all.rename(columns={'Conditions': 'Population'})

# tsne_all["rep"] = all_conds_rep_df["Condition"].tolist()

# tsne_all_esm.to_csv("tSNE_all_reps_esm.csv", index = False)
# tsne_all_rosetta.to_csv("tSNE_all_reps_rosetta.csv", index = False)

# tsne_all = pd.read_csv("tSNE_all_reps.csv")

# Funcs_output_ana.Ploting_tSNE_panel(tsne_all_esm, out = "esm", column_name = "rep")
# Funcs_output_ana.Ploting_tSNE_panel(tsne_all_rosetta, out = "rosetta", column_name = "rep")


# Funcs_output_ana.Ploting_tSNE(tsne_all, out = "tsne_dG", color_method="dG")
# Funcs_output_ana.Ploting_tSNE(tsne_all, out = "tsne_population", color_method="population")

# Funcs_output_ana.Ploting_tSNE_panel(tsne_all, out = "tsne_population", column_name = "rep")
# Funcs_output_ana.Ploting_tSNE_panel(tsne_all, out = "teste_bugado", column_name = "rep", plot_lowest_dG=True)
# Funcs_output_ana.Plotting_tSNE_panel_lines(tsne_all, out = "teste_bugado", column_name = "rep")

# tsne_all["Population"]

# tsne_all[tsne_all["rep"] == 1]




# Funcs_output_ana.Ploting_tSNE(tsne_all[tsne_all["rep"] == 1], out = "dg_rep1", color_method="dG")
# Funcs_output_ana.Ploting_tSNE(tsne_all[tsne_all["rep"] == 2], out = "dg_rep2", color_method="dG")
# Funcs_output_ana.Ploting_tSNE(tsne_all[tsne_all["rep"] == 3], out = "dg_rep3", color_method="dG")

# Funcs_output_ana.Ploting_tSNE(tsne_all[tsne_all["rep"] == 1], out = "pop_rep1", color_method="population")
# Funcs_output_ana.Ploting_tSNE(tsne_all[tsne_all["rep"] == 2], out = "pop_rep2", color_method="population")
# Funcs_output_ana.Ploting_tSNE(tsne_all[tsne_all["rep"] == 3], out = "pop_rep3", color_method="population")


# def histogram(df):
#     # Criar o histograma
#     plt.hist(df["dG"], bins=10, edgecolor='black')

#     # Títulos e rótulos
#     plt.title('Histograma Exemplo')
#     plt.xlabel('Valores')
#     plt.ylabel('Frequência')

#     # Mostrar o gráfico
#     plt.show()


# histogram(rep1)
# histogram(rep2)
# histogram(rep3)



# # Standardize the data
# scaler = StandardScaler()
# X_scaled = scaler.fit_transform(tsne_all.iloc[:, :2])

# # Perform PCA
# pca = PCA(n_components=2)  # We are interested in the first 2 components
# X_pca = pca.fit_transform(X_scaled)

# x=0
# explained_variance = [(x+pca_features.explained_variance_[i]/np.sum(pca_features_all.explained_variance_))*100 for i in range(dimensions)]

# # Print explained variance for the first 2 components
# explained_variance = pca.explained_variance_ratio_
# print(f'Explained variance for the first 2 components: {explained_variance}')

# 15150/3

# 5050*2



# def pca_dataset(dataset, dimensions):

#     pca_features= PCA(n_components=dimensions)
#     principalComponents_features = pca_features.fit_transform(dataset)

#     pca_features_all = PCA(n_components= min(np.shape(dataset)))
#     principalComponents_features_all = pca_features_all.fit_transform(dataset)

#     x=0
#     explained_variance = [(x+pca_features.explained_variance_[i]/np.sum(pca_features_all.explained_variance_))*100 for i in range(dimensions)]
    
#     return principalComponents_features, explained_variance

# teste, variance = pca_dataset(dataset = all_vhse, dimensions = 2)







# tsne = Funcs_output_ana.red_tsne(data = all_vhse, dimension = 2, prp = 20, label_list = all_sequences["Condition"].tolist())

# tsne_all = tsne_all.rename(columns={'Conditions': 'Population'})
# tsne_all["rep"] = 0
# tsne_all.iloc[0:2800:,3] = 1
# tsne_all.iloc[2801:-1:,3] = 2
# Funcs_output_ana.Ploting_tSNE(tsne_rep1)
# Funcs_output_ana.Ploting_tSNE(tsne_rep2)
# Funcs_output_ana.Ploting_tSNE(tsne_all)






# last_cycles = all_Sequences[all_Sequences["Cycle"] == str(150)]

# id_matrix_rep1 = Funcs_output_ana.calculate_identity_matrix(rep1)

# tsne = Funcs_output_ana.red_tsne(data = id_matrix, dimension = 2, prp = 20, col_cond = all_best["Condition"].tolist())


# Funcs_output_ana.Ploting_tSNE(tsne)

# import umap

# reducer = umap.UMAP()



# myata = pd.read_csv("miyata_table.csv", index_col=None)
# myata.rename(columns={'Unnamed: 0': 'Aminoacid'}, inplace=True)
# myata.set_index('Aminoacid', inplace = True)

# bla1, bla2, bla3 = Funcs_output_ana.List_min_dgs(condition_1, condition_2, condition_3)


# dist_matrix_cond3 = Funcs_output_ana.Distance_matrix_init_final(df = condition_3, out = "teste")



# min(condition_2[condition_2["dG"]])

# min(condition_2["dG"])

# teste1, teste2 = Funcs_output_ana.algn_pairwise(seq1 = min_c1["Sequence"], seq2 = min_c2["Sequence"])

# id_teste = Funcs_output_ana.calculate_indentity(seq1 = min_c1["Sequence"], seq2 = min_c2["Sequence"])


# len(teste[1])

# teste1 = Funcs_output_ana.compute_identical_matches(align1 = teste[0], align2 = teste[1])

# min_c1 = Funcs_output_ana.Get_best_ind(condition_1)

# min_c2 = Funcs_output_ana.Get_best_ind(condition_2)


# min_c3 = Funcs_output_ana.Get_best_ind(condition_3)

# dada = pd.concat([min_c1, min_c2, min_c3])

# dada["Sequence"][0]

# modeling_functions_pyrosetta.exc





# min_c3["Sequence"].name










# min_c3.index()





# lowests_cond3 = [min(condition_3[condition_3["Cycle"] == str(i)]["dG"]) for i in range(len(set(condition_3["Cycle"])))]

# min(condition_3[condition_3["Cycle"] == str(0)]["dG"])
    
    
# condition_3['dG']

# min(condition_3[condition_3["Cycle"] == str(0)]float(["dG"]))

# plot1, = plt.plot(bla1, linestyle='-', label='Natural Sequences')
# plot1, = plt.plot(bla2, linestyle='-', label='Natural Sequences')
# plot1, = plt.plot(bla3, linestyle='-', label='Natural Sequences')

# # Your list of values
# data1 = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]
# data2 = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]

# # Plot the first set of values with a label
# plot1, = plt.plot(data1, marker='o', linestyle='-', label='Data Points 1')

# # Plot the second set of values with a label
# plot2, = plt.plot(data2, marker='s', linestyle='--', label='Data Points 2')

# # Add labels and title
# plt.xlabel('X-axis label')
# plt.ylabel('Y-axis label')
# plt.title('Plot of Values')

# # Manually create handles for the legend
# legend_handles = [plot1, plot2]

# # Add a legend with specified handles
# plt.legend(handles=legend_handles)

# # Show the plot
# plt.show()



# df["Column3"]

# len(set(df["Column3"]))

# lowests_cond1 = [float(min(df[df["Column3"] == str(i)]["Column2"])) for i in range(len(set(df["Column3"])))]


# min(df[df["Column3"] == '1']["Column2"])

# # Plot the values
# plt.plot(lowests_cond1, marker='o', linestyle='-')

# # Add labels and title
# plt.xlabel('Cycle')
# plt.ylabel('Best_dG')
# plt.title('lowest_dG_per_cycle (Condition 1)')

# # Show the plot
# plt.show()


min_dG_row = unique_apt_rosetta_random.loc[unique_apt_rosetta_random['dG'].idxmin()]




