from collections import namedtuple
import argparse
import os
import math
import sys
sys.path.append("/storage/mustoe/software/StructureAnalysisTools-Dev")
from ReactivityProfile import ReactivityProfile
import numpy as np
import matplotlib
#Avoid using windows backend as this was written for the TACO computational cluster
#at BCM
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import scipy.stats

#Custom named record for the purpose of housing reactivities
normalized_tuple = namedtuple('normalized_profile_values', 'A C U G G_N7')


#Extract and sort profile names from the user submitted profile folder.
#Always keeps central profile first in the sorted order.
def get_profiles(prof_1, folder):
    profile_strings = os.listdir(folder) 
    profile_strings = list(filter(lambda name: name[-4:] == ".txt", profile_strings))

    profile_strings.remove(prof_1)
    profile_strings.sort(key = lambda name: name.lower())
    profile_strings = [prof_1] + profile_strings
         
    return profile_strings
    
#Extract and return reactivity tuple for each profile 
def get_reactivities(profile_strings, folder):
    
    reactivity_tuples = []
    for profile in profile_strings:
        p = ReactivityProfile(folder + profile)
        p_ga = ReactivityProfile(folder + profile + "ga")

        
        reactivities = []
        for nt in "ACUG":
            reactivities.append(list(p.profile('norm')[p.sequence == nt]))

        reactivities.append(list(p_ga.profile('norm')[p_ga.sequence == "G"]))
        reactivity_tuples.append(normalized_tuple(A=reactivities[0], C=reactivities[1], U=reactivities[2], G=reactivities[3], G_N7=reactivities[4]))


        
    return reactivity_tuples

#Calculate linear regression, pearson correlation of two lists of reactivities.
#Plot information to a given subplot.
def plot_reactivity(ax, original_x, y, x_title, y_title = None):
    
    x = [elem for elem in original_x]
    #print("x_title: {}".format(x_title))
    if x_title == "norm G_N7":
   #     print("-----------------------")
   #     print("x: ", x)
   #     print("y: ", y)
   #     print("Pop")
        for i in range(len(x)):
            if math.isinf(x[i]):
                x.pop(i)
                y.pop(i)
        for i in range(len(y)):
            if math.isinf(y[i]):
                x.pop(i)
                y.pop(i)
    #    print("x: ", x)
    #    print("y: ", y)

    #calculate linear regression 
    regress = scipy.stats.linregress(x, y)

    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)

    gmin = min(xmin, ymin)
    gmax = max(xmax, ymax)

    #Plot points
    ax.scatter(x, y)

    #Plot the diagonal
    ax.plot([gmin, gmax], [gmin, gmax], 'k--', label='diagonal')

    #Plot the fit
    x_fit = np.linspace(xmin, xmax, num=100)
    y_fit = x_fit * regress.slope + regress.intercept
    ax.plot(x_fit, y_fit, 'r', label='fit')

    #Add labels
    if y_title != None:
        ax.set_ylabel(y_title, fontsize = 12)
    ax.set_title(x_title +  ": R: {0:.3} M: {1:.3}".format(regress.rvalue, regress.slope))


#Generate and plot correlation, slope, and scatter plots for the cartesian product of the main
#profile and all other profiles.
def plot_all_reactivities(profile_strings, normalized_tuple_list, out_prefix):
    main_file = profile_strings[0]
    secondary_files = profile_strings[1:]

    main_reactivity = normalized_tuple_list[0]
    secondary_reactivity = normalized_tuple_list[1:]



    
    x = 5
    y = len(secondary_files)

    fig, axes = plt.subplots(y, x, figsize=(25, 5 * y))

    for x_iter in range(x):
        for y_iter in range(y):
            ax = axes[y_iter, x_iter]
            x_title = "norm " + normalized_tuple._fields[x_iter] 
            y_title = secondary_files[y_iter] 

            y_reac = secondary_reactivity[y_iter]


            if x_iter == 0:
                plot_reactivity(ax, main_reactivity[x_iter], y_reac[x_iter], x_title, y_title)

            else: 
                plot_reactivity(ax, main_reactivity[x_iter], y_reac[x_iter], x_title)

    
    axes[y - 1, 2].set_xlabel(main_file, fontsize = 12)
    plt.tight_layout()
    plt.savefig("{}.pdf".format(args.out_prefix))


#Remove any NT in the normalized tuple list that contains a NaN value as this
#will interfere with downstream correlation calculations.
def remove_NaN(normalized_tuple_list, profile_strings):
    nan_columns = {"A":[], "C":[], "U":[], "G":[], "G_N7":[]}


    for i, base in enumerate(normalized_tuple._fields):

        for nt in range(len(normalized_tuple_list[0][i])):
            nt_values = [normalized_tuple_list[j][i][nt] for j in range(len(normalized_tuple_list))]
            if any(math.isnan(elem) for elem in nt_values):
                nan_columns[base].append(nt)

    filtered_normalized_tuple_list = []
    for norm_tuple in normalized_tuple_list:
        filtered_normalized_tuple_list.append(normalized_tuple(A =  [ nt for index, nt in enumerate(norm_tuple[0]) if index not in nan_columns["A"] ],
                                                               C =  [ nt for index, nt in enumerate(norm_tuple[1]) if index not in nan_columns["C"] ],
                                                               U =  [ nt for index, nt in enumerate(norm_tuple[2]) if index not in nan_columns["U"] ],
                                                               G =  [ nt for index, nt in enumerate(norm_tuple[3]) if index not in nan_columns["G"] ],
                                                               G_N7 = [ nt for index, nt in enumerate(norm_tuple[4]) if index not in nan_columns["G_N7"] ]))
        
    return filtered_normalized_tuple_list 



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prof_1", help = "Name of central profile file for comparison (Including .txt and without path.)")
    parser.add_argument("--folder", help = "Full path to folder housing all profile.txt and profile.txtga files.")
    parser.add_argument("--out_prefix", help = "Prefix of output file (not including .pdf)")
    args = parser.parse_args()

    if args.folder[-1] != "/":
        args.folder += "/"


    profile_strings = get_profiles(args.prof_1, args.folder)
    normalized_tuple_list = get_reactivities(profile_strings, args.folder)
    normalized_tuple_list = remove_NaN(normalized_tuple_list, profile_strings)
    plot_all_reactivities(profile_strings, normalized_tuple_list, args.out_prefix)
