import argparse
import matplotlib.pyplot as plt
import ReactivityProfile as rp
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
import numpy as np



def load_profiles(profile_files):
    profile_list = []
    #ga_profile_list = []
    for f in profile_files:
        profile_list.append(rp.ReactivityProfile(f))
    #for f in ga_profile_files:
    #    ga_profile_list.append(rp.ReactivityProfile(f))

    #return profile_list, ga_profile_list
    return profile_list

""" Both of these functions are key
    functions for downstream map 
    calls.
"""
def average_key_function(x):
    if len(x) == 0:
        return np.nan
    else:
        return sum(x) / len(x)


def std_dev_key_function(x):
    #Note, this calculates the population std dev. Not the sample std dev.
    if len(x) == 0:
        return np.nan
    else:
        return float(np.std(x))

#def calculate_average_nt_reactivity(profile_list, ga_profile_list):
def calculate_average_nt_reactivity(profile_list):
    avg_reac = []
    #avg_N7 = []
    #std_dev = []
    #std_dev_N7 = []

    for NT in profile_list[0].sequence:
        avg_reac.append([])
        #avg_N7.append([])
        #std_dev.append([])
        #std_dev_N7.append([])
        

    for profile in profile_list:
        for NT in range(len(profile.normprofile)):
            if np.isnan(profile.normprofile[NT]):
                pass
            else:
                avg_reac[NT].append(profile.normprofile[NT])
               # std_dev[NT].append(profile.normprofile[NT])

    #for profile in ga_profile_list:
    #    for NT in range(len(profile.normprofile)):
    #        if np.isnan(profile.normprofile[NT]):
    #            pass
    #        else:
    #            avg_N7[NT].append(profile.normprofile[NT])
    #            std_dev_N7[NT].append(profile.normprofile[NT])
            


    #avg_N13 = list(map(average_key_function, avg_N13))
    #avg_N7 = list(map(average_key_function, avg_N7))
#    std_dev = list(map(std_dev_key_function, std_dev))
#    std_dev_N7 = list(map(std_dev_key_function, std_dev_N7))

    #return avg_N13, avg_N7, std_dev, std_dev_N7
    #return list(map(average_key_function, avg_reac))
    return [average_key_function(reac) for reac in avg_reac]

def calculate_nt_stdev(profile_list):
    reacs = []

    for NT in profile_list[0].sequence:
        reacs.append([])

    for profile in profile_list:
        for NT in range(len(profile.normprofile)):
            if np.isnan(profile.normprofile[NT]):
                pass
            else:
                reacs[NT].append(profile.normprofile[NT])

    return np.array([std_dev_key_function(reac) for reac in reacs])

def average_profile(profile_list):
    avg_profile = rp.ReactivityProfile()
    profile_list = load_profiles(profile_list)
    avg_profile.normprofile = np.array(calculate_average_nt_reactivity(profile_list))
    avg_profile.sequence = profile_list[0].sequence
    return avg_profile

def calc_stdev(profile_list):
    profile_list = load_profiles(profile_list)
    return calculate_nt_stdev(profile_list)

#def color_reactivity(scores, threshold, color):
#    """Returns a list of colors corresponding to scores position in the threshold. 
#        
#        scores: rate(s) to be colored
#        threshold: Threshold of colors in descending order.
#    """
#
#    colors = []
#    for score in scores:
#        if score > threshold[0]:
#            colors.append(color[0])
#        elif score < threshold[0] and score > threshold[1]:
#            colors.append(color[1])
#        else:
#            colors.append(color[2])
#
#    return colors
    

#def plot_subplot(avg_N13, avg_N7, std_dev, std_dev_N7):
#    """ Plots subplots.
#    """
#
#    filtered_N13 = [elem for elem in avg_N13 if elem]
#    filtered_std_dev = [elem for elem in std_dev if elem]
#    N13_indexes = [i + 1 for i in range(len(avg_N13)) if avg_N13[i]]
#    N13_colors = color_reactivity(filtered_N13, [.4, .2], ['red', 'orange', 'black'])
#
#
#    for i in range(len(avg_N7)):
#        if not avg_N7[i]:
#            continue
#        elif avg_N7[i] > 10:
#            avg_N7[i] = 10
#
#    filtered_N7 = [elem for elem in avg_N7 if elem]
#
#    N7_colors = color_reactivity(filtered_N7, [5,3], ["purple", "pink", "black"])
#
#    filtered_N7 = [ -1 * elem for elem in filtered_N7]
#
#    filtered_std_dev_N7 = [elem for elem in std_dev_N7 if elem]
#    N7_indexes = [i + 1 for i in range(len(avg_N7)) if avg_N7[i]]
#
#    fig, axs = plt.subplots(2,1, sharex=False, figsize=(.09 * len(avg_N13), 2.8), gridspec_kw={'hspace': 0})
#    axs_top = axs[0]
#    axs_bottom = axs[1]
#
#    axs_top.bar(N13_indexes, filtered_N13, yerr = filtered_std_dev, width=1, color = N13_colors)
#
#
#    axs_top.yaxis.set_ticks([0.0, .2, .4])
#    axs_top.yaxis.set_ticks([.6, .8, 1.0, 1.2, 1.4], minor = True, )
#    axs_top.yaxis.set_ticklabels([.6,.8,1.0,1.2,1.4], minor = True)
#    axs_top.set_ylim(bottom=0, top=1.5) 
#    axs_top.grid(True, axis="y", which = "major")
#    axs_top.tick_params(which='both', width=2)
#    axs_top.set_ylabel("N7 --- N1/3")
#    axs_top.get_yaxis().set_label_coords(-0.02, 0)
#
#    axs_top.set_title("Average Profile Comparison", size=14)
#
#    
#    axs_bottom.set_xlabel("NT(s)")






#    axs_bottom.bar(N7_indexes, filtered_N7, yerr = filtered_std_dev_N7, width=1, color = N7_colors)
#    axs_bottom.set_ylim(bottom=-10, top=0) 
#    axs_bottom.set_yticks([-3, -5, -10])
#    axs_bottom.tick_params(which='both', width=2)
#    axs_bottom.grid(True, axis="y")
#
#    axs_bottom.add_patch(Rectangle((-.5, -3), 1, 3, facecolor = 'black', edgecolor = 'none'))
#    axs_bottom.add_patch(Rectangle((-.5, -5), 1, 2, facecolor = 'pink', edgecolor = 'none'))
#    axs_bottom.add_patch(Rectangle((-.5, -10), 1, 5, facecolor = 'purple', edgecolor = 'none'))
#    axs_top.set_xlim(left=-.5)
#    axs_bottom.set_xlim(left=-.5)
#
#    axs_top.add_patch(Rectangle((-.5, 0), 1, .2, facecolor = 'black', edgecolor = 'none'))
#    axs_top.add_patch(Rectangle((-.5, 0.2), 1, .2, facecolor = 'orange', edgecolor = 'none'))
#    axs_top.add_patch(Rectangle((-.5, 0.4), 1, 1.2, facecolor = 'red', edgecolor = 'none'))
#
#    plt.tight_layout()


#def plot_avg_values(avg_N13, avg_N7, std_dev, std_dev_N7, output_prefix):
#    """ Constructs subplots and formats main plot.
#    """
#    plot_subplot(avg_N13, avg_N7, std_dev, std_dev_N7) 
#
#    plt.savefig(output_prefix + ".pdf", format = "pdf")


#if __name__ == "__main__":
#
#    print("Plotting")
#
#    parser = argparse.ArgumentParser()
#    h =  "Two or more profile.txt files to plot average norm. reactivity / stnd dev for\n"
#    h += "in a nt specific manner."
#    parser.add_argument("--profiles", nargs = '+', help=h)
#    h =  "Two or more profile.txtga files to plot average norm. reactivity / stnd dev for\n"
#    h += "in a nt specific manner."
#    parser.add_argument("--ga_profiles", nargs = '+', help=h)
#    parser.add_argument("--output_prefix", help="Title of output (excluding .pdf)")
#    
#    args = parser.parse_args()
#    profile_list, ga_profile_list = load_profiles(args.profiles, args.ga_profiles)
#    avg_N13, avg_N7, std_dev, std_dev_N7 = average_profiles(profile_list, ga_profile_list)
#    plot_avg_values(avg_N13, avg_N7, std_dev, std_dev_N7, args.output_prefix)
