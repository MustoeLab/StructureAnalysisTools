#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import numpy as np
import scipy.stats
import argparse 
 
import sys, os
from ReactivityProfile import ReactivityProfile


def filter_profiles(profile1, profile2, nts, name, exclhighbg):
    # lets check to make sure the sequences are the same
    if not np.array_equal(profile1.sequence, profile2.sequence):
        raise ValueError('Sequences are not the same!')
        
    
    # initiliaze my mask (all false)
    mask = np.zeros(len(profile1.sequence), dtype=bool)
    
    if nts is None:
        mask[:] = True #reassign all values to True
    
    else:
        for n in nts:
            mask = mask | (profile1.sequence == n)
    
    if exclhighbg is not None:
        with np.errstate(invalid='ignore'):
            mask = mask & (profile1.backprofile < exclhighbg) & (profile2.backprofile < exclhighbg)
    
    
    # get the desired profiles
    r1 = profile1.profile(name)     
    r2 = profile2.profile(name)
    
    mask = mask & np.isfinite(r1) & np.isfinite(r2)

    return r1[mask], r2[mask]

def filterProfilesByNt(profile1, profile2, profile1_N7=None, profile2_N7=None, nts=None, name=None, exclhighbg=None):
    """Return matched reactivities from ReactivityProfile objects prof1 and prof2,
    filtering out any NaNs
    If nts argument is None, return all nts. Otherwise, return only these nts   
    name argument is passed to ReactivityProfile to get desired profile    
    """

    # If extracting N7_G just use use old functionality
    if nts == "N7_G":
        profile1 = profile1_N7
        profile2 = profile2_N7
        nts = "G"
    
    filtered_r1, filtered_r2 = filter_profiles(profile1, profile2, nts, name, exclhighbg)

    # If it's a list and N7_G is to be included or user is running with default argument (all)
    if (isinstance(nts, list) and "N7_G" in nts) or (nts is None and (profile1_N7 or profile2_N7)):
        if not (profile1_N7 and profile2_N7):
            raise ValueError("Must pass profile1_N7 and profile2_N7 to run filter_profiles with N7_G positions")

        filtered_r1_N7, filtered_r2_N7 = filter_profiles(profile1_N7, profile2_N7, "G", name, exclhighbg)

        filtered_r1 = np.append(filtered_r1, filtered_r1_N7)
        filtered_r2 = np.append(filtered_r2, filtered_r2_N7)
    
    return filtered_r1, filtered_r2
    
    

def plotCorrelation(var1, var2, ax, title='', xlabel='', ylabel=''):
    if len(var1) == 0 and len(var2) == 0:
        return
    elif sum(var1) == 0 and sum(var2) == 0:
        print("Norm profile not present. Skipping.")
        return
    
    regress = scipy.stats.linregress(var1, var2)
    
    # get min and max values for plotting
    xmin, ymin = min(var1), min(var2)
    gmin = min(xmin, ymin)
    xmax, ymax = max(var1), max(var2)
    gmax = max(xmax, ymax)
    
    ax.scatter(var1,var2)
    
    # plot the diagonal
    ax.plot([gmin, gmax], [gmin, gmax], 'k--', label='diagonal')

    # plot the fit
    x = np.linspace(xmin, xmax, num=100)
    y = x*regress.slope + regress.intercept
    ax.plot(x, y, 'r', label='fit')
    
    # add labels and stuff
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('{} : R={:.3f} ; m={:.1f}'.format(title, regress.rvalue, regress.slope))
    
    
#############################################################################

    


if __name__=='__main__':
    
    parser=argparse.ArgumentParser()
    parser.add_argument('profile1', help='Path to first profile file')
    parser.add_argument('profile2', help='Path to second profile file')
    parser.add_argument('output_prefix', help='Path to output file (Will save as pdf)')
    parser.add_argument("--N7_profile1", help = "Path to first N7 profile file. Should correspond to profile1. (Will add additional N7_G column to analysis.)")
    parser.add_argument("--N7_profile2", help = "Path to second N7 profile file. Should correspond to profile2.")
    parser.add_argument('--exclhighbg', type=float, help='Exclude nts with bg values greater than this cutoff')
    parser.add_argument('--comparison', default='all', help='Type of comparison to perform. Options are raw/sub/back/norm/all (default=all)')
    parser.add_argument("--mindepth", type=int, default = 100, help="Any nucleotide less than this depth is not plotted in any graph. (default = 100)")
    parser.add_argument("--include_in_all", nargs="+", help = "What nucleotides to include in the 'All' column. (Default: A C U G [G_N7 - if N7 included] )")
    
    args=parser.parse_args()

    # Input validation
    if args.comparison not in ('raw','sub', 'back', 'norm','all'):
        exit('Invalid comparison selected. Options are raw/sub/back/norm/all')
    
    if bool(args.N7_profile1) != bool(args.N7_profile2):
        exit("Error! must pass both N7_profile* flags, or none at all.")

    if args.N7_profile1:
        for N1_prof, N7_prof in zip((args.profile1, args.profile2), (args.N7_profile1, args.N7_profile2)):
            N1_prof = os.path.basename(N1_prof)
            N7_prof = os.path.basename(N7_prof)
            if N1_prof.replace(".txt", ".txtga") != N7_prof:
                print(f"Warning! N1 profile {N1_prof} and N7 profile {N7_prof} seem to have come from separate shapemapper runs. Proceed with caution.")

    if args.include_in_all:
        for nt in args.include_in_all:
            if nt not in {"A", "C", "U", "G", "N7_G"}:
                exit(f"Error! {nt} is not a valid nucleotide! Valid nucleotides include: A C U G or G_N7. These arguments must be space delimited. EG --include_in_all A C U")
            elif nt == "N7_G" and args.N7_profile1 is None:
                exit("Error! Including N7_G in all column plots, but not passing N7 profiles!")

    # Setting what we include in all by default
    if args.include_in_all is None:
        if not (args.N7_profile1 is None): # If N7 included
            args.include_in_all = ["A", "C", "U", "G", "N7_G"]
        else:
            args.include_in_all = ["A", "C", "U", "G"]

    # read in the desired files
    p1 = ReactivityProfile(args.profile1, depthcut = args.mindepth)
    p2 = ReactivityProfile(args.profile2, depthcut = args.mindepth)


    # Removing renormalization currently as we should be plotting out the values SM outputs

    #for p in (p1, p2):
    #    p.normprofile[np.isnan(p.subprofile)] = np.nan
    #if args.comparison in ('all','norm') and not args.GA:
    #    if not np.all(np.isnan(p1.subprofile)):
    #        p1.normalize(eDMS=True)
    #    if not np.all(np.isnan(p2.subprofile)):
    #        p2.normalize(eDMS=True)


    # create the matplotlib figure object
    if args.N7_profile1 is None:
        if args.comparison == 'all':
            fig, ax = plot.subplots(4, 5, figsize=(20,16))

            for i,name in enumerate(('raw','sub','back','norm')):
       
                # this loop will compute correlations for all nts combined, and then individually
                for j, n in enumerate( ('All', 'A', 'C', 'U', 'G')):

                       if n=='All':
                           x,y = filterProfilesByNt(p1, p2, nts=args.include_in_all, name=name, exclhighbg=args.exclhighbg)
                       else:
                           x,y = filterProfilesByNt(p1, p2, nts=n, name=name, exclhighbg=args.exclhighbg)
                       
               
                       plotCorrelation(x,y, ax[i,j], title=name+' '+n)
             
            # label the axes
            ax[3,2].set_xlabel(args.profile1.split('/')[-1])
            ax[1,0].set_ylabel(args.profile2.split('/')[-1])

        else:
        
            fig, ax = plot.subplots(1, 5, figsize=(20,4))

            # this loop will compute correlations for all nts combined, and then individually
            for i, n in enumerate( ('All', 'A', 'C', 'U', 'G')):

                if n=='All':
                   x,y = filterProfilesByNt(p1, p2, nts=args.include_in_all, name=args.comparison, exclhighbg=args.exclhighbg)
                else:
                   x,y = filterProfilesByNt(p1, p2, nts=n, name=args.comparison, exclhighbg=args.exclhighbg)
            
                plotCorrelation(x,y, ax[i], title=n)

            # label the axes
            ax[2].set_xlabel(args.profile1.split('/')[-1])
            ax[0].set_ylabel(args.profile2.split('/')[-1])



    # N7 plotting functionality
    else:
        # Initialize N7 reactivity profile
        p1_N7 = ReactivityProfile(args.N7_profile1, depthcut = args.mindepth)
        p2_N7 = ReactivityProfile(args.N7_profile2, depthcut = args.mindepth)

        if args.comparison == 'all':
            fig, ax = plot.subplots(4, 6, figsize=(24,16))

            for i,name in enumerate(('raw','sub','back','norm')):
       
                # this loop will compute correlations for all nts combined, and then individually
                for j, n in enumerate( ('All', 'A', 'C', 'U', 'G', 'N7_G')):

                    if n=='All':
                        x,y = filterProfilesByNt(p1, p2, p1_N7, p2_N7, nts=args.include_in_all, name=name, exclhighbg=args.exclhighbg)
                    else:
                        x,y = filterProfilesByNt(p1, p2, p1_N7, p2_N7, nts=n, name=name, exclhighbg=args.exclhighbg)
                    
           
                    plotCorrelation(x,y, ax[i,j], title=name+' '+n)
             
            # Label axes
            ax[3,2].set_xlabel(args.profile1.split('/')[-1] + "(ga)")
            ax[1,0].set_ylabel(args.profile2.split('/')[-1] + "(ga)")

        else:

            fig, ax = plot.subplots(1, 6, figsize=(24,4))

            # this loop will compute correlations for all nts combined, and then individually
            for i, n in enumerate( ('All', 'A', 'C', 'U', 'G', 'N7_G')):

                if n=='All':
                   x,y = filterProfilesByNt(p1, p2, p1_N7, p2_N7, nts=args.include_in_all, name=args.comparison, exclhighbg=args.exclhighbg)
                else:
                   x,y = filterProfilesByNt(p1, p2, p1_N7, p2_N7, nts=n, name=args.comparison, exclhighbg=args.exclhighbg)
            
                plotCorrelation(x,y, ax[i], title=n)

            # label the axes
            ax[2].set_xlabel(args.profile1.split('/')[-1] + "(ga)")
            ax[0].set_ylabel(args.profile2.split('/')[-1] + "(ga)")

    #make a bit prettier
    fig.tight_layout()
    # save the data in pdf format
    out_file = args.output_prefix
    if not out_file.endswith(".pdf"): out_file += ".pdf"
    fig.savefig(out_file)
