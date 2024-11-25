#!/usr/bin/python

# ---------------------------------------------------------------------------------------
# Flexible Arc Plotting code
# Anthony Mustoe
# Weeks lab, UNC
# 2016, 2017, 2018
# 
# Modifications
# Version 2.1 with plot profile option
#
# ---------------------------------------------------------------------------------------


import sys, os, math, argparse
import pandas as pd
import RNAStructureObjects as RNAtools
import numpy as np
import matplotlib 
matplotlib.use('Agg')
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['xtick.major.width'] = 2.5
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['xtick.minor.size'] = 4 
matplotlib.rcParams['xtick.minor.width'] = 1
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Arial'
import matplotlib.pyplot as plot
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from ReactivityProfile import ReactivityProfile
from mean_reactivity_stdev import average_profile, calc_stdev

#import warnings

#warnings.filterwarnings("ignore")


class ArcLegend(object):
    """Container for Legend for arc plot"""
    
    def __init__(self, title=None, colors=[], labels=[], msg=None):
        
        assert len(colors)==len(labels), "colors and labels must be the same size" 

        self.title=title
        self.colors = colors
        self.labels = labels
        self.msg = msg
    
    
    def append(self, t, colors, labels):
        
        if self.title is not None and t is not None:
            self.title+=' & '+t
        elif self.title is None:
            self.title = t

        self.colors.extend(colors)
        self.labels.extend(labels)


    def drawLegend(self, ax, xbound, ybound):
        """Draw legend on axes, computing location based on axes bounds
        ax     = axes object to draw on
        xbound = xboundary of axis
        ybound = yboundary of axis
        """
        
        spacing = 3
        
        # determine bounding box
        nlines = len(self.colors) + int(self.title is not None) + 1
        height = spacing*nlines       
        width = 15
        if self.title is not None and len(self.title) > width:
            width = len(self.title)
        
        xloc = xbound
        yloc = ybound
        

        # write the title
        if self.title is not None:
            ax.text(xloc, yloc, self.title, horizontalalignment='left',
                    size="11",weight="medium", color="black")
            yloc -= spacing


        for i, c in enumerate(self.colors):

            # self.colors is a list of dict patch specifications, so c is a dict
            ax.add_patch(patches.Rectangle( (xloc, yloc), 1.5, 1.5, color=c, clip_on=False ) )

            # now add label
            ax.text(xloc+2.5, yloc, self.labels[i], horizontalalignment='left',
                    size="8",weight="normal", color="black")
            yloc -= spacing
    
        if self.msg:
            ax.text(xloc, yloc, self.msg, horizontalalignment='left',
                    size="8",weight="normal", color="black")



class ArcPlot(object):

    def __init__(self, title='', fasta=None):
        
        self.title = title

        self.seq = ''
        if fasta:
            self.addFasta(fasta)


        self.seqcolors = ['black']
        self.drawseq = True
        self.grid = False

        self.handleLength = 4*(math.sqrt(2)-1)/3
        self.adjust = 1.4

        self.topPatches = []
        self.botPatches = []
        self.height = [0, 0] # the figure height

        self.intdistance = None

        self.reactprofile = None
        self.reactprofileType = 'SHAPE'

        self.N7profile = None
        self.lower_N7_Plot = None

        self.toplegend = None
        self.botlegend = None

        self.upper_error = None
        self.lower_error = None
        
        self.annotation = [] # list of annotation information


    def addFasta(self, fastafile):
        
        read = False
        with open(fastafile) as inp:
            for line in inp:
   
                if line[0]=='>':
                    read = True
                else:
                    line = line.rstrip()
                    if read and len(line) > 0:
                        self.seq += line
                        
        self.length = len(self.seq)

    def addArcPath(self, outerPair, innerPair=None, panel=1, color = 'black', alpha=0.5, window=1):
        """ add the arPath object for a given set of parameters
        If only outerPair is passed, innerPair is computed from window
        If both outerPair and innerPair are passed, window is not used
        """
        
        if innerPair is None:
            innerPair = [outerPair[0]+0.5 + window-1, outerPair[1]-0.5]
            outerPair = [outerPair[0]-0.5, outerPair[1]+0.5 + window-1]
        else:
            innerPair = [innerPair[0]+0.5, innerPair[1]-0.5]
            outerPair = [outerPair[0]-0.5, outerPair[1]+0.5]
            
        innerRadius = (innerPair[1] - innerPair[0])/2.0
        outerRadius = (outerPair[1] - outerPair[0])/2.0
        
        verts = []
        codes = []

        # outer left
        verts.append( (outerPair[0], 0) )
        codes.append( Path.MOVETO )

        # outer left control 1
        verts.append( (outerPair[0], panel*self.handleLength*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer left control 2
        verts.append( (outerPair[0]+outerRadius*(1-self.handleLength), panel*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer center
        verts.append( (outerPair[0]+outerRadius, panel*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer right control 1
        verts.append( (outerPair[0]+outerRadius*(1+self.handleLength), panel*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer right control 2
        verts.append( (outerPair[1], panel*self.handleLength*outerRadius) )
        codes.append( Path.LINETO )
                       
        # outer right
        verts.append( (outerPair[1], 0) )
        codes.append( Path.LINETO )
                        
        # inner right
        verts.append( (innerPair[1], 0) )
        codes.append( Path.LINETO )
            
        # inner right control 1
        verts.append( (innerPair[1], panel*self.handleLength*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner right control 2
        verts.append( (innerPair[0]+innerRadius*(1+self.handleLength), panel*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner center
        verts.append( (innerPair[0]+innerRadius, panel*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner left control 1 
        verts.append( (innerPair[0]+innerRadius*(1-self.handleLength), panel*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner left control 2
        verts.append( (innerPair[0], panel*self.handleLength*innerRadius) )
        codes.append( Path.LINETO )
 
        # inner left
        verts.append( (innerPair[0], 0) )
        codes.append( Path.LINETO )
 
        # outer left duplicate 
        verts.append( (outerPair[0], 0) )
        codes.append( Path.CLOSEPOLY )
 
        # rescale verts

        if panel == 1:
            adval = self.adjust
        else:
            adval = -(self.adjust-1)
        
        # move arcs away from x-axis
        verts = [ (x,y+adval) for x,y in verts ]
         

        indpath = Path(verts, codes)
        patch = patches.PathPatch(indpath, facecolor=color, alpha=alpha,
                                  linewidth=0, edgecolor='none')
        
        if panel == 1:
            self.topPatches.append(patch)
            if outerRadius > self.height[1]:
                self.height[1] = outerRadius
        
        else:
            self.botPatches.append(patch)
            if outerRadius > self.height[0]:
                self.height[0] = outerRadius



    def plotProfile(self, ax, bounds=None, colthresh = (-10, 0.4, 0.85, 3), heightscale=None, N7=False, bar_height_scale = 1.0):
        """Add a reactivity profile to the axes (self.reactprofile needs to be set)
        
        ax          = axes object to add plot. Expected to be top axis unless N7.
        colthresh   = thresholds to bin reactivity data. 
                      First element indicates "No-data lower bound"
                      Last element indicates "Max reactivity" -- values are trucated above this
        heightscale = scaling to control height of bars
        """
        
        if self.reactprofile is None:
            return
        
        # check to see if DMS and colthresh defaults not overridden
        if self.reactprofileType == 'DMS' and colthresh == (-10,0.4,0.85,3):
            colthresh = (-10, 0.2, 0.4, 1)
            

        xvals = [ [] for i in range(4) ]
        yvals = [ [] for i in range(4) ]
        
        if heightscale is None:
            heightscale = max(4, min(10, len(self.reactprofile)/50.))
            heightscale = min(max(self.height)/4., heightscale)
            
            if self.reactprofileType == 'DMS': # adjust for compressed ploting range
                heightscale *= 2

        heightscale *= bar_height_scale

        if N7: 
            if self.reactprofileType == 'DMS': # Scale plot to better fit with DMS
                heightscale *= 0.30303030303030304

            spltSeq = list(self.seq)
            for x,y in enumerate(self.N7profile):
                if bounds is not None and x<bounds[0] or x>bounds[1]:
                    continue

                if y is None or y != y  or y<colthresh[0]:
                    if spltSeq[x] in ['g', 'G']:
                        if y != -3.3: # Visualize small values
                            xvals[1].append(x+1)
                            yvals[1].append(.2 * heightscale)
                        else:
                            xvals[0].append(x+1)
                            yvals[0].append(-1)
                elif y < colthresh[1]:
                    xvals[1].append(x+1)
                    yvals[1].append(y * heightscale)
                elif y < colthresh[2]:
                    xvals[2].append(x+1)
                    yvals[2].append(y * heightscale)
                else:
                    xvals[3].append(x+1)
                    if y > colthresh[3]:
                        yvals[3].append(colthresh[3]*heightscale)
                    else:
                        yvals[3].append(y * heightscale)
        else:
            for x,y in enumerate(self.reactprofile):

                if bounds is not None and x<bounds[0] or x>bounds[1]:
                    continue

                if y is None or y != y or y<colthresh[0]:
                    xvals[0].append(x+1)
                    yvals[0].append(0.6*heightscale+self.adjust)
                elif y < colthresh[1]:
                    xvals[1].append(x+1)
                    yvals[1].append(y*heightscale)
                elif y < colthresh[2]:
                    xvals[2].append(x+1)
                    yvals[2].append(y*heightscale)
                else:
                    xvals[3].append(x+1)
                    if y > colthresh[3]:
                        yvals[3].append(colthresh[3]*heightscale)
                    else:
                        yvals[3].append(y*heightscale)

        
        if N7:
        # If the values are N7 values, make them negative
            yvals[0] = [elem * -1 for elem in yvals[0]] 
            yvals[1] = [elem * -1 for elem in yvals[1]] 
            yvals[2] = [elem * -1 for elem in yvals[2]] 
            yvals[3] = [elem * -1 for elem in yvals[3]] 
            ax.bar(xvals[0], -.1, alpha=0.7, linewidth=0, color=(179./255, 171./255, 148./255), align='center', clip_on=False)
            ax.bar(xvals[1], yvals[1], alpha=0.7, linewidth=0, color='black', align='center', clip_on=False)
            ax.bar(xvals[2], yvals[2], alpha=0.7, linewidth=0, color='hotpink', align='center', clip_on=False, width = 1.3)
            ax.bar(xvals[3], yvals[3], alpha=0.7, linewidth=0, color='darkviolet', align='center', clip_on=False, width = 1.3)
            colthresh = [elem * -1 for elem in colthresh]

            ax.axes.get_yaxis().set_visible(True)
            ax.tick_params(axis='y', direction='out', labelsize=6, left=True, right=False)
            ax.set_yticks( np.array(colthresh) * heightscale)

            labels = [str(x) for x in colthresh]
            labels[-1] = '>'+labels[-1]
            ax.set_yticklabels( labels )

            ax.set_frame_on(True)
            for l in ('right','top','bottom'):
                ax.spines[l].set_visible(False)
               
            ax.spines['left'].set_bounds(self.adjust, colthresh[3] * heightscale)

        else:
            ax.bar(xvals[0], yvals[0], alpha=0.7, linewidth=0, color=(179./255, 171./255, 148./255),
                   align='center', clip_on=False, bottom=-0.3*heightscale)
        
            ax.bar(xvals[1], yvals[1], alpha=0.7, linewidth=0, color='black',
                   align='center', clip_on=False, bottom=self.adjust)
            ax.bar(xvals[2], yvals[2], alpha=0.7, linewidth=0, color='orange',
                   align='center', clip_on=False, bottom=self.adjust)
            ax.bar(xvals[3], yvals[3], alpha=0.7, linewidth=0, color='red',
                   align='center', clip_on=False, bottom=self.adjust)

            ax.axes.get_yaxis().set_visible(True)
            ax.tick_params(axis='y', direction='out', labelsize=6, left=True, right=False)
            ax.set_yticks( np.array(colthresh[1:])*heightscale+self.adjust )
            
            labels = [str(x) for x in colthresh[1:]]
            labels[-1] = '>'+labels[-1]
            ax.set_yticklabels( labels )
            ax.set_ylabel('Norm. React.', position=(0,0), ha='left', size=6)
               
            ax.set_frame_on(True)
            for l in ('right','top','bottom'):
                ax.spines[l].set_visible(False)
               
            ax.spines['left'].set_bounds(self.adjust, colthresh[3]*heightscale+self.adjust)




    
    def highlightLowReactivity(self, profile, cutoff=0.0001, window=1):
        """Profile is a np array (0-indexed) of values"""
        
        for i,r in enumerate(profile):
           
            if r<cutoff:
                patch = patches.Rectangle((i+0.5,0),window,4,linewidth=0,fill=True,alpha=0.2)
                self.topPatches.append(patch)     
                patch = patches.Rectangle((i+0.5,-4),window,4,linewidth=0,fill=True,alpha=0.2)
                self.botPatches.append(patch)
                
                
    def addANNO(self, ANNO, length, panel, colors=None):
        
        from math import sqrt
        """Add annotations as bars"""
        
        if colors:
            colors = colors.split(',')
        else:
            colors = ['darkblue', 'darkgreen', 'darkred', 'darkorange', 'darkmagenta']
            
        groups = ANNO['group'].unique()
        colordict = {}
        for group, col in zip(groups, colors):
            colordict[group] = col
            
        ANNO['color'] = ANNO['group'].apply(lambda g: colordict[g])
        
        anno_count = 2
        plot_ratio_used = 0.1
        while plot_ratio_used < 0.9:
            if len(ANNO['label'].unique()) <= anno_count:
                break
            else:
                plot_ratio_used += 0.1
                anno_count += 2
        else:
            plot_ratio_used = 0.9
         
        plot_area_height = max(self.height)*plot_ratio_used
        height_per_anno = plot_area_height/len(ANNO['label'].unique())
        annotation_half_width = height_per_anno*0.2
        annotation_text_offset = height_per_anno * 0.25
        
        label = ''
        index = -1
        
        for idx, row in ANNO.iterrows():
            if row['label'] != label:
                index += 1
                label = row['label']

            start = row['tstart']
            end = row['tend']
            coordinate = [(start, (-index-1)*height_per_anno + annotation_half_width), \
                          (start, (-index-1)*height_per_anno - annotation_half_width), \
                          (end, (-index-1)*height_per_anno - annotation_half_width), \
                          (end, (-index-1)*height_per_anno + annotation_half_width), \
                          (0, 0)]
            codes = [Path.MOVETO, \
                     Path.LINETO, \
                     Path.LINETO, \
                     Path.LINETO, \
                     Path.CLOSEPOLY]

            # move arcs away from x-axis
            if panel == 1:
                adval = max(self.height)*0.025
            else:
                adval = -max(self.height)*0.025
            coordinate = [(x,y+adval) for x,y in coordinate]
            
            anno_path = Path(coordinate, codes)
            anno_patch = patches.PathPatch(anno_path, facecolor=row['color'], \
                                           linewidth = 0, label = row['group'], alpha = 0.8)
            
            annotation_x_coor = (start + end)/2
            annotation_y_coor = (-index-1)*height_per_anno + annotation_text_offset + adval
            annotation = (label, annotation_x_coor, annotation_y_coor, row['color'], height_per_anno)
            self.annotation.append(annotation)

            if panel == 1:
                self.topPatches.append(anno_patch)
                self.height[1] = length*sqrt(len(ANNO))*0.02
            else:
                self.botPatches.append(anno_patch)
                self.height[0] = length*sqrt(len(ANNO))*0.02
                
        # Add the legend
        t = 'Group'
        c = list(colordict.values())
        l = list(colordict.keys())
        
        if panel>0:
            self.toplegend = ArcLegend(title=t, colors=c, labels=l)
        else:
            self.botlegend = ArcLegend(title=t, colors=c, labels=l) 
             




    def writePlot(self, outPath="arcs.pdf", bounds=None, write=True,
                  msg=None, msg_pos=(0,1), msg_kwargs={}, bar_height_scale=1.0, **args):
        
        cutbounds = True
        if bounds is None:
            bounds = (0,len(self.seq))
            cutbounds = False
        else:
            bounds = [bounds[0]-0.5, bounds[1]+0.5] # make new copy

        
        doubleplot = len(self.botPatches)>0
        scaleFactor = 0.05
        
        figWidth = (bounds[1]-bounds[0])*scaleFactor
        figHeight = 2*max(self.height)*scaleFactor
        

        fig = plot.figure( figsize=(figWidth, figHeight)) # 500*scaleFactor
        fig.subplots_adjust(hspace=0.0)

        if self.lower_N7_Plot == True and not doubleplot:
            axT = fig.add_subplot(211)
            axB = None
            doubleplot = True
            self.addArcPath((1,2), panel=-1, alpha = 0)

    
        else:
            axT = fig.add_subplot(211)
            axB = None

        for pat in self.topPatches:
            axT.add_patch(pat)
            pat.set_clip_on(cutbounds)

        axT.set_xlim(bounds)
        axT.set_ylim(0, max(self.height))

        axT.set_frame_on(False)
        axT.axes.get_yaxis().set_visible(False)
        
        if self.toplegend:
            self.toplegend.drawLegend(axT, bounds[0], max(20, self.height[1]*.5))

        if doubleplot:
            axB = fig.add_subplot(212, sharex=axT)
            for pat in self.botPatches:
                axB.add_patch(pat)
                pat.set_clip_on(cutbounds)

            axB.set_xlim(bounds)
            axB.set_ylim(-max(self.height), 0)
            axB.set_frame_on(False)
            axB.axes.get_yaxis().set_visible(False)
            axB.tick_params(axis='both', which='both', top=False, bottom=False, labelbottom=False)
            
            if self.botlegend:
                self.botlegend.drawLegend(axB, bounds[0], -self.height[1]*.5) 



        # format sequence ticks
        axT.get_xaxis().tick_bottom()  
        
        majorLocator = MultipleLocator(100)
        minorLocator = MultipleLocator(20)
        tickinterval = 5

        axT.xaxis.set_major_locator(majorLocator)
        axT.xaxis.set_major_formatter( FormatStrFormatter('%i') )
        axT.xaxis.set_minor_locator(minorLocator)
        axT.xaxis.set_minor_formatter( FormatStrFormatter('%i') )

        majlabels = axT.axes.get_xaxis().get_majorticklabels()
        minlabels = axT.axes.get_xaxis().get_minorticklabels()
        majticks = axT.axes.get_xaxis().get_major_ticks()
        minticks = axT.axes.get_xaxis().get_minor_ticks()
        

        for label in majlabels:
            label.set_size(10)
        
        if bounds[0] == 0: 
            majlabels[1].set_visible(False)
            minlabels[1].set_visible(False)
            majticks[1].set_visible(False)
            minticks[1].set_visible(False)
        
        beginint, r = divmod(bounds[0], 20)
        if r == 0:
            beginint+=1

        for i, label in enumerate(minlabels):
            label.set_size(7)
            # Below fix likely only necessary when running in py2 removes duplicate ticks
            if sys.version_info[0] < 3:
                if i % tickinterval == beginint: 
                    label.set_visible(False)


        # plot gridlines
        if self.grid:
            axT.grid(b=True, which="minor", color='black',linestyle='--', alpha=0.5)
            axT.grid(b=True, which="major", color='black',linestyle='-', lw=2, alpha=0.8)

            if doubleplot:
                axB.grid(b=True, which="minor", color='black',linestyle='--', alpha=0.5)
                axB.grid(b=True, which="major", color='black',linestyle='-', lw=2, alpha=0.8)



        # write title (structure names)
        if self.title:
            axT.text(0.0, self.height[1]-2, self.title, horizontalalignment='left',
                      size="24",weight="bold", color="blue")


        # put nuc sequence on axis 
        if self.drawseq:
            fontProp = matplotlib.font_manager.FontProperties(family = "sans-serif", 
                                              style="normal",
                                              weight="extra bold",
                                              size="4")
            
            for i, nuc in enumerate(self.seq):
                
                if i<bounds[0]-1:
                    continue
                elif i>bounds[1]-1:
                    break

                if nuc == "T":
                    nuc = "U"

                try:
                    col = self.seqcolors[i]
                except IndexError:
                    col = "black"
                
                axT.annotate(nuc, xy=(i+0.5, 0),fontproperties=fontProp,color=col,
                             annotation_clip=False,verticalalignment="baseline")

        # add annotation for eCLIP data
        # each element is a tuple = (label, x-coor, y-coor, color, fontsize)
        for annotation in self.annotation:
            axB.text(annotation[1], annotation[2], annotation[0], \
                     color = annotation[3], fontsize = annotation[4], \
                     horizontalalignment='center', verticalalignment='bottom')                

        if self.reactprofile is not None:
            self.plotProfile(axT, bounds, bar_height_scale = bar_height_scale, **args)
    
        if self.lower_N7_Plot == True:
            self.plotProfile(axB, bounds, colthresh = (0, 1.6, 2.3, 3.3), N7=True, **args)

        # Plot st dev for NT averages across multiple profiles
        if self.upper_error is not None:
            self.plot_error_bars(axT, self.upper_error, self.reactprofile, bar_height_scale = bar_height_scale)
        if self.lower_error is not None:
            self.plot_error_bars(axB, self.lower_error, self.N7profile, N7=True)
        
        if msg is not None:
            axT.text(msg_pos[0], msg_pos[1], msg, transform=axT.transAxes, **msg_kwargs)


        # save the figure
        if write and outPath != 'show':
            fig.savefig(outPath, dpi=100, bbox_inches="tight", transparent=True)
            plot.close(fig)
        
        elif write and outPath == 'show':
            plot.show()

        else:
            return fig, axT, axB 



    def addCT(self, ctObj, color='black', alpha=0.7, panel=1):
    
        seqlen = len(ctObj.ct)

        if self.seq == '' or self.seq[0] == ' ':
            self.seq = ctObj.seq
        elif len(self.seq) != seqlen:
            print("Warning:: CT length = {0}; expecting length = {1}".format(seqlen, len(self.seq)))
        

        i = 0
        while i < seqlen:
            
            if ctObj.ct[i] > i:

                outerPair = [ i+1, ctObj.ct[i] ]

                # find the right side of helix
                lastPairedNuc = ctObj.ct[i]
                i += 1
                while i<seqlen and (lastPairedNuc-ctObj.ct[i]) == 1:
                    lastPairedNuc = ctObj.ct[i]
                    i += 1
                i -= 1


                innerPair = [ i+1, ctObj.ct[i] ]

                self.addArcPath(outerPair, innerPair, panel=panel, color=color, alpha=alpha)

            i += 1


    
    def addRings(self, ringfile, N7=False, N17=False, metric='alpha', panel=1, bins=None, contactfilter=(None,None),
                 filterneg=False, skip_legend=False):
        """Add arcs from ringmapper file
        metric  = z, sig or alpha
        default: alpha
        """
        
        if contactfilter[0]:
            print("Performing contact filtering with distance={}".format(contactfilter[0]))
        
        if filterneg:
            print("Filtering out negative correlations")

        colors = [(44,123,182), (44,123,182), (171,217,233), (255,255,255), (253,174,97), (215,25,28), (215,25,28)]


        if N7 or N17:
            colors = [(91, 218, 190), (91, 218, 190), (0, 95, 154), (255, 255, 255), (218, 91, 172), (83, 48, 95), (83, 48, 95)]

        if metric == "alpha":
            colors = colors[3:]

        colors = [tuple([y/255.0 for y in x]) for x in colors]
        
        if metric == "alpha":
            alpha = [0.0, 1.0, 1.0]
            gradiant = [False, True, False]

        else:
            alpha = [1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
            gradiant = [False, True, False, False, True, False]
        
        if metric=='z':
            if bins is None:
                bins = [-50, -5, -1, 0, 1, 5, 50]
            else:
                bins = [-1e10, -bins[1], -bins[0], 0, bins[0], bins[1], 1e10]
        elif metric=='alpha':
            if bins is None:
                bins = [0, 2, 5, 1e10]
            else:
                bins = [0, bins[0], bins[1], 1e10]

        else:
            if bins is None:
                bins = [-1e10, -100, -20, 0, 20, 100, 1e10]
            else:
                bins = [-1e-10, -bins[1], -bins[0], 0, bins[0], bins[1], 1e10]

        allcorrs = []
         
        with open(ringfile) as inp:
            header = inp.readline().split()
            window = int(header[1].split('=')[1])
            mname = header[2].split('=')[1]
            inp.readline()
                
            # define molecule length if not set from something else
            try:
                len(self.seq)/len(self.seq)
            except:
                self.seq = ' '*int(header[0])
                print("self.seq: ", self.seq)
 

            for line in inp:
                spl = line.split()
                i = int(spl[0])
                j = int(spl[1])

                if contactfilter[0] and contactfilter[1].contactDistance(i,j) <= contactfilter[0]:
                    continue
                
                if filterneg and int(spl[3])<0:
                    continue

                if metric == 'alpha' and float(spl[2]) < 20:
                    continue

                if self.seq[i - 1] not in ['A','G','C','U','T'] or self.seq[j - 1] not in ['A','G','C','U','T']:
                    continue


                if metric == 'z':
                    val = float(spl[4])
                elif metric =='alpha':
                    val = float(spl[9])
                else:
                    val = float(spl[2])

                if val > bins[0] and metric == "alpha":
                    allcorrs.append( ( int(spl[0]), int(spl[1]), float(spl[3])*val ) )

                elif val > bins[3]:
                    allcorrs.append( ( int(spl[0]), int(spl[1]), float(spl[3])*val ) )
        

        if len(allcorrs)==0:
            print('WARNING: No RINGs passing filter in {}'.format(ringfile))

        for i in range(len(bins)-1):
            
            corrs = [c for c in allcorrs if bins[i]<c[2]<=bins[i+1]]
            corrs.sort(key=lambda x:x[2]) 
           
            for c in corrs:
                if gradiant[i]:
                    col = self._colorGrad(c[2], colors[i], colors[i+1], bins[i], bins[i+1])
                else:
                    col = colors[i]
                
                self.addArcPath( c[:2], panel=panel, color=col, alpha=alpha[i], window=window)

        # Add the legend
        t = 'RING {0} win={1} {2}'.format(mname, window, metric.upper())

        if metric == "alpha":
            c = (colors[1], colors[2])
            l = ('>{}'.format(bins[1]), '>{}'.format(bins[2]))

        else:
            if filterneg:
                c = (colors[4], colors[5])
                l = ('>{}'.format(bins[4]), '>{}'.format(bins[5]))
            else:
                c = (colors[1], colors[2], colors[4], colors[5])
                l = ('<{}'.format(bins[1]), '<{}'.format(bins[2]), 
                     '>{}'.format(bins[4]), '>{}'.format(bins[5]))
        

        if not skip_legend:
            if panel>0:
                if self.toplegend is not None:
                    self.toplegend.append(t, c, l)
                else:
                    self.toplegend = ArcLegend(title=t, colors=c, labels=l)
            else:
                if self.botlegend is not None:
                    self.botlegend.append(t, c, l)
                else:
                    self.botlegend = ArcLegend(title=t, colors=c, labels=l)



    def addPairProb(self, dpObj, panel=1, bins=None):
        """ add pairing probability arcs from a dotplot object"""
 

        if self.seq == '':
            self.seq = ' '*dpObj.length
        elif len(self.seq) != dpObj.length:
            print("Warning: dp file sequence length = {0}; CT/FASTA length = {1}".format(dpObj.length, len(self.seq)))


        refcolors = [ (150,150,150), (255,204,0),  (72,143,205) ,(81, 184, 72) ]
        refalpha = [0.55, 0.65, 0.75, 0.85]

        gradiant = [False, False, False, False]


        if bins is None:
            bins = [0.03, 0.1, 0.3, 0.8, 1.0]
            colors = refcolors
            alpha = refalpha
        else:
            if len(bins) > 5:
                raise IndexError('{0} PairProb levels specified -- the maximum allowed is 4'.format(len(bins)-1))

            colors = refcolors[-(len(bins)-1):]
            alpha = refalpha[-(len(bins)-1):]


        # convert colors to percents
        colors = [tuple([y/255.0 for y in x]) for x in colors]
        

        for i in range(len(bins)-1):
            
            cutdp = dpObj.requireProb(bins[i], bins[i+1])
            cutpairs = cutdp.pairList()
            
            if len(cutpairs) == 0:
                continue
            
            # sort each of the pairs by its strength
            strength = []
            for p in cutpairs:
                elem = cutdp.partfun[p[0]-1]
                pos = ( elem['pair'] == p[1]-1) # make mask
                strength.append(elem['log10'][pos][0]) # get value

            strengthzip = list(zip(cutpairs, strength))
            strengthzip.sort(key=lambda x: x[1], reverse=True)

            cutpairs, strength = list(zip(*strengthzip))
            
            # iterate through sorted pairs and draw arcs
            for pindex, pair in enumerate(cutpairs):
 
                if gradiant[i]:
                    col = self._colorGrad(strength[pindex], colors[i], colors[i+1], 
                                          bins[i], bins[i+1], log=True)
                else:
                    col = colors[i]

                self.addArcPath(pair, panel=panel, color=col, alpha=alpha[i])


        # Add the legend
        t = 'Pairing Prob.'
        c = colors
        l = ['>{}'.format(x) for x in bins[:-1]]
        
        if panel>0:
            self.toplegend = ArcLegend(title=t, colors=c, labels=l)
        else:
            self.botlegend = ArcLegend(title=t, colors=c, labels=l) 
       




    def addPairMap(self, pmobj, panel=1, plotall=False):
        """plot pairs output by pairmapper
        plotall = True will plot all complementary correlations (not just 1&2)
        """

        if len(pmobj.primary)==0 and len(pmobj.secondary)==0:
            print('WARNING: No PAIR-MaP correlations in {}'.format(pmobj.sourcename))


        def getZ(corrlist):
            return [(x[0], x[1], x[3]) for x in corrlist]
        
        colors = [(100, 100, 100), (30,194,255), (0,0,243)]
        colors = [tuple([x/255.0 for x in c]) for c in colors]
 
        
        if len(self.seq)==0:
            self.seq = ' '*pmobj.molsize


        if plotall:
            self.plotAlphaGradient( getZ(pmobj.remainder), colors[0], (0.0,0.8),
                                    1, 6, window=pmobj.window, panel=panel)
        
        self.plotAlphaGradient( getZ(pmobj.secondary), colors[1], (0.2,0.6),
                                2, 6, window=pmobj.window, panel=panel)
        
        self.plotAlphaGradient( getZ(pmobj.primary), colors[2], (0.5,0.9),
                                2, 6, window=pmobj.window, panel=panel)
        
        
        # make legend
        c = colors[::-1]  # reverse colors so primary on top
        l = ['Principal','Minor','Not passing']
        if not plotall:
            c = c[:-1]
            l = l[:-1]
        
        if panel>0:
            if self.toplegend is not None:
                self.toplegend.append('PairMap', c, l)
            else:
                self.toplegend = ArcLegend(title='PairMap', colors=c, labels=l)
        else:
            if self.botlegend is not None:
                self.botlegend.append('PairMap', c, l)
            else:
                self.botlegend = ArcLegend(title='PairMap', colors=c, labels=l)



    def _colorGrad(self, value, colorMin, colorMax, minValue, maxValue, log=False):
        """ return an interpolated rgb color value """
        
        if log:
            value = 10**-value
            minValue = 10**-minValue
            maxValue = 10**-maxValue

        if value > maxValue:
            return colorMax
        elif value < minValue:
            return colorMin
        else:
            v = value - min(maxValue, minValue)
            v /= abs(maxValue-minValue)
        
            col = []
            for i in range(3):
                col.append( v*(colorMax[i] - colorMin[i]) + colorMin[i] )
            
            return col


    def plotAlphaGradient(self, pairlist, color, alpha_range, lb, ub, window=1, panel=1):
        """plot the list of pairs with the same color but different alpha"""

        if max(color) > 1:
            color = [x/255. for x in color]

        dynamicrange = ub-lb

        for i,j,val in sorted(pairlist, key=lambda x:x[2]):
            
            if val<lb:
                continue
            elif val>ub:
                scale = 1
            else:
                scale = (val-lb)/dynamicrange

            a = (alpha_range[1]-alpha_range[0])*scale+alpha_range[0]
            
            self.addArcPath( (i,j), window=window, color=color, panel=panel, alpha=a)



    def compareCTs(self, refCT, compCT, panel=1):
        
        if len(refCT.ct) != len(compCT.ct):
            raise IndexError('CT objects are different sizes')

        share = refCT.copy()
        refonly = refCT.copy()
        componly = refCT.copy()
        
        # ct files and set to blank (ie 0)
        for c in (share, refonly, componly):
            for i in range(len(c.ct)):
                c.ct[i] = 0
    

        for i in range(len(refCT.ct)):
            
            if refCT.ct[i] == compCT.ct[i]:
                share.ct[i] = refCT.ct[i]
            else:
                refonly.ct[i] = refCT.ct[i]
                componly.ct[i] = compCT.ct[i]
            
        sharedcolor = (150/255., 150/255., 150/255.)
        refcolor = (38/255., 202/255., 145/255.)
        compcolor = (153/255., 0.0, 1.0)
        self.addCT(share, color=sharedcolor, panel=panel)
        self.addCT(refonly, color=refcolor, panel=panel)
        self.addCT(componly, color=compcolor, panel=panel)
 

        sens,ppv,nums = refCT.computePPVSens(compCT, False)
        msg = 'Sens={0:.2f} PPV={1:.2f}'.format(sens, ppv)
        print(msg)

        if panel>0:
            self.toplegend = ArcLegend(colors=[sharedcolor,refcolor,compcolor], 
                                       labels=['Shared', refCT.name, compCT.name],
                                       msg=msg)                  
        else:
            self.botlegend = ArcLegend(colors=[sharedcolor,refcolor,compcolor], 
                                       labels=['Shared', refCT.name, compCT.name],
                                       msg=msg)                  
       
        return sens,ppv
 

    def assignSeqColors(self, shapearr):

        if len(shapearr) != len(self.seq):
            raise IndexError("Shape Array does not match length of sequence")

        col = []       
        for x in shapearr:

            if x < -4: 
                col.append( (160,160,160 ) )  # grey
            elif x > 0.85: 
                col.append( (255,0,0) )  # red
            elif 0.85 >= x >0.4: 
                col.append( (255,164,26) ) # orange
            else: 
                col.append( (0,0,0) ) # black

        self.seqcolors = [tuple([y/255.0 for y in x]) for x in col]


    
    def colorSeqByMAP(self, mapfile):
        
        shape, seq = RNAtools.readSHAPE(mapfile)
        
        if self.seq == '' or self.seq.count(' ')==len(self.seq):
            self.seq = seq
            
        elif len(shape) != len(self.seq):
            raise IndexError("Mapfile does not match length of sequence!")
        
        self.assignSeqColors(shape)
 
    
    def readSHAPE(self, mapfile):

        shape, seq = RNAtools.readSHAPE(mapfile)
        
        if self.seq == '' or self.seq.count(' ') == len(self.seq):
            self.seq = seq
        
        elif len(shape) != len(self.seq):
            raise IndexError("Mapfile does not match length of sequence!")
        
        self.reactprofile = shape
        
        # need to give the plot some height if no other things added
        if max(self.height)==0:
            self.height[1] = max(5, min(10, len(self.reactprofile)/50.))



    def readProfile(self, profilefile, dms=False):

        if isinstance(profilefile, ReactivityProfile):
            self.reactprofile = profilefile.normprofile
            if not isinstance(self.seq, np.ndarray) and (self.seq == '' or self.seq.count(' ') == len(self.seq)):
                self.seq = profilefile.sequence

        else:
            with open(profilefile) as inp:
                line = inp.readline().split()
                if len(line) == 2 or len(line)==4:
                    ftype=1
                else:
                    ftype=2
            

            if ftype==1:
                self.readSHAPE(profilefile)
            else:
                profile = ReactivityProfile(profilefile)
                self.reactprofile = profile.normprofile
                
                if self.seq == '' or self.seq.count(' ') == len(self.seq):
                    self.seq = profile.sequence

        if dms:
            self.reactprofileType = 'DMS'
        
        # need to give the plot some height if no other things added
        if max(self.height)==0:
            #self.height[1] = max(5, min(10, len(self.reactprofile)/50.))
            self.height[1] = max(20, min(40, len(self.reactprofile)/12.5))


    def readN7Profile(self, N7File, panel=-1):

        if isinstance(N7File, ReactivityProfile):
            profile = N7File   
        else:
            profile = ReactivityProfile(N7File)

        react = profile.normprofile

        for ite in range(len(react)):
            if (np.isnan(react[ite])):
                react[ite] = -3.3
        

        if (panel == -1):
            self.lower_N7_Plot=True

        self.N7profile = react

    def splitCatRing(self, ringFile):
        """Splits up concatenated ring file into three arrays by type of ring (N1/3, N7, N1/3-7)"""
        oFile = open(ringFile, "r")
        allLines = [line for line in oFile]
        oFile.close()

        #Removes header
        arcLines = allLines[2:]

        #Appends header to start of growing line lists
        N13Lines = allLines[0:2]
        N7Lines = allLines[0:2]
        N17Lines = allLines[0:2]
        

        #Calculates midpoint given this is a concatenated file
        midpt = len(self.seq)
        

        #Iterates through and sort lines by type
        for index in range(len(arcLines)):
            splLine = arcLines[index].split()

            i = int(splLine[0])
            j = int(splLine[1])

            if ((i > midpt and j <= midpt ) or (i <= midpt and j > midpt)):
                splLine[1] = str(int(splLine[1]) - midpt )

                if int(splLine[0]) > int(splLine[1]):
                    splLine[0], splLine[1] = splLine[1], splLine[0]

                adjustLine = "\t".join(splLine) + "\n"
                N17Lines.append(adjustLine)
            

            elif (i > midpt and j > midpt):
                splLine[0] = str(int(splLine[0]) - midpt )
                splLine[1] = str(int(splLine[1]) - midpt )
                adjustLine = "\t".join(splLine) + "\n"
                N7Lines.append(adjustLine)

            elif (i <= midpt and j <= midpt):
                N13Lines.append(arcLines[index])

        return N13Lines, N7Lines, N17Lines


    def addCatRings(self, ringfile, N7=False, N17=False, metric='alpha', panel=1, bins=None, contactfilter=(None,None),filterneg=False):
        """Writes concatenated ring files to plot. Splits files into N13 on top and N7 / N17 below"""

        
        #Sanity Checking
        if(not list(self.seq)):
            print("Error: No sequence specified. Must add profile or fasta before adding concatenated ring file.")
        else:
            
            #Splits file into an array of lines based on type
            N13Arr, N7Arr, N17Arr = self.splitCatRing(ringfile)
            
            

            #Produces temp files
            N13File = open(".N13_file_temp_ring.txt", "w")
            for line in N13Arr:
                N13File.write(line)
            N13File.close()

            N7File = open(".N7_file_temp_ring.txt", "w")
            for line in N7Arr:
                N7File.write(line)
            N7File.close()

            N17File = open(".N17_file_temp_ring.txt", "w")
            for line in N17Arr:
                N17File.write(line)
            N17File.close()


            #Write to plot
            self.addRings(".N13_file_temp_ring.txt", N7=False, N17=False, metric=metric, panel=1, bins=bins, contactfilter=contactfilter, filterneg=filterneg)

            if N7:
                self.addRings(".N7_file_temp_ring.txt", N7=N7, metric=metric, panel=-1, bins=bins, contactfilter=contactfilter, filterneg=filterneg)

            if N17:
                if N7:
                    self.addRings(".N17_file_temp_ring.txt", N17=N17, metric=metric, panel=-1, bins=bins, contactfilter=contactfilter, filterneg=filterneg, skip_legend=True)
                else:
                    self.addRings(".N17_file_temp_ring.txt", N17=N17, metric=metric, panel=-1, bins=bins, contactfilter=contactfilter, filterneg=filterneg)

            #Removes temp files
            os.remove(".N17_file_temp_ring.txt")
            os.remove(".N7_file_temp_ring.txt")
            os.remove(".N13_file_temp_ring.txt")

    # Record std dev in error bar attribute for average profile arcplot runs
    def add_error_bars(self, nt_stdev, lower = False):
        if lower:
            self.lower_error = nt_stdev
        else:
            self.upper_error = nt_stdev

    def plot_error_bars(self, ax, error_size, error_bar_y, heightscale = None, N7 = False, bar_height_scale = 1.0):
        """Add error bars to corresponding subplot (top or bottom).
           Currently only used to plot standard deviation for NT
           locations informed by multiple profiles."""

        # plotProfile does not plot full bar if greater than the max colorThresh. Only plots up to max of colorThresh.
        # Need to account for this so error bars plotted at proper location in graph
        if N7:
            cutoff = 3.3
        elif self.reactprofileType == 'DMS':
            cutoff = 1
        else:
            cutoff = 3

        error_bar_y[error_bar_y > cutoff] = cutoff


        # Transform error_bar_y locations as well as error bar scaling depending on type
        # of profile that is being annotated.
        if N7:
            error_bar_y *= -1
        if heightscale is None:
            heightscale = max(4, min(10, len(self.reactprofile)/50.))
            heightscale = min(max(self.height)/4., heightscale)
            if self.reactprofileType == 'DMS':
                heightscale *= 2
                if N7:
                    heightscale *= 0.30303030303030304


        heightscale *= bar_height_scale
        error_bar_y *= heightscale
        if not N7:
            error_bar_y += self.adjust
        error_size *= heightscale

        ax.errorbar(range(1, len(error_size) + 1), error_bar_y, xerr = [0] * len(error_size), yerr = error_size, linestyle='none', linewidth=0.5, ecolor = 'black')


class ANNO():
    """ANNO object to process functional annotations"""
    
    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.df = pd.DataFrame()
    
    def get_rna_info(self, bedfile):
        self.start_coors = []
        self.end_coors = []
        self.exons = [0]
        
        with open(bedfile, 'r') as f:
            for line in f:
                spl = line.rstrip().split('\t')
                self.chr = spl[0]
                self.start_coors.append(int(spl[1]))
                self.end_coors.append(int(spl[2]))
                self.strand = spl[3]

                exon = int(spl[2]) - int(spl[1])
                self.exons.append(exon) # use this exon length to update transcript coordinate
         
        self.exons = self.exons[:-1]
        
        self.start = min(self.start_coors)
        self.end = min(self.end_coors)
        self.length = self.end - self.start + 1
        
        if 'chr' not in self.chr:
            self.chr = 'chr' + self.chr
    
    # if use "bound" flag, remove peaks outside of these bounds and clip peaks within bound
    def apply_bound(self, bound):

        lower_bound = bound[0]
        upper_bound = bound[1]
            
        temp = self.df[(self.df['tend']<lower_bound)|(upper_bound<self.df['tstart'])]
        self.df = self.df[~self.df.index.isin(temp.index)].reset_index(drop=True)
        self.df['tstart'] = self.df['tstart'].clip(lower = lower_bound)
        self.df['tend'] = self.df['tend'].clip(upper = upper_bound)
        self.df = self.df.dropna()
    
    def read_anno(self, anno, bound = None):
        
        names = ['tstart', 'tend', 'label', 'group']
        temp = pd.read_csv(anno, sep = '\t', names = names)
        temp.fillna('', inplace = True) # fill N/A values with an empty string
        self.df = pd.concat([self.df, temp])
    
    # function to find features that are mapped to the provided genomic coordinates
    def intersect(self, df, start, end, added_length):
        temp = df[df['chr']==self.chr] #check for same chromosome
        if self.strand:
            temp = temp[temp['strand']==self.strand] #check for same strand
            
        not_intersect = temp[(end<temp['start'])|(temp['end']<start)]
        intersectdf = temp[~temp.index.isin(not_intersect.index)].reset_index(drop=True)
        if not intersectdf.empty:
            intersectdf['tstart'] = intersectdf['start'] - start
            intersectdf['tstart'] = intersectdf['tstart'] + added_length

            intersectdf['tend'] = intersectdf['end'] - start
            intersectdf['tend'] = intersectdf['tend'] + added_length
        
        return intersectdf

    
    def read_annodir(self, annodir, cutoff, annolist):
        col_names = ['chr', 'start', 'end', 'info', '1000', 'strand', 'score1', 'score2']
        
        df = pd.DataFrame()
        added_length = 0
        for start, end, exon in zip(self.start_coors, self.end_coors, self.exons):
            added_length += exon
            for f in os.listdir(annodir):
                if '.bed' in f:
                    fn = '{}/{}'.format(annodir,f) # f'{annodir}/{f}'  
                    if '//' in fn:
                        fn = fn.replace('//', '/')

                    temp = pd.read_csv(fn, sep = '\t', names = col_names, 
                                       usecols = [i for i in range(8)],
                                      dtype = {'start': int, 'end': int})
                    
                    # deal with eCLIP files that does not have a proper label
                    skip = False
                    if '.' in temp['info'].unique():
                        skip = True
                        
                    if skip:
                        continue
                            
                    out = self.intersect(temp, start, end, added_length) # apply intersect function
                    if not out.empty:
                        df = pd.concat([df, out])

        if len(df) == 0:
            print('No rows intersect with {} in the reference file(s). Exiting'.format(self.name))
            exit()
            
        # apply score filter
        if cutoff:
            score1, score2 = [float(i) for i in cutoff.split(',')]
            df = df[(df['score1']>=score1)&(df['score2']>=score2)]
            
        df['tstart'] = df['tstart'].clip(lower = 0) # replace all negative values to lower bound
        df['tend'] = df['tend'].clip(upper = self.length) # set maximum value to upper bound

        if self.strand == '-':
            newstart = self.length - df['tend']
            newend = self.length - df['tstart']
            df['tstart'] = newstart
            df['tend'] = newend

        df.reset_index(drop=True, inplace = True)
        df.to_csv('{}_map.csv'.format(self.name), sep = '\t', index = False)

        count = len(df)
        print('{} features in reference files(s) mapped to provided coordinates. See output file {}_map.csv'.format(count, self.name))

        def split_info(info):
            spl = info.split('_')
            if len(spl) > 1:
                return spl[0], spl[1]
            elif len(spl) == 1:
                return spl[0], 'group'

        df['label'], df['hue'] = zip(*df['info'].apply(lambda g: split_info(g)))

        def merge_peak(group):
            rbp = group['label'].values[0]
            group['diff'] = group['tstart'] - group['tend'].shift(1)

            labels = []
            no = 0
            for idx, row in group.iterrows():
                if row['diff'] is np.nan:
                    labels.append(0)
                elif row['diff'] <= 0:
                    labels.append(no)
                else:
                    no += 1
                    labels.append(no)

            group['diff'] = labels
            return group

        df = df.sort_values(by = ['info', 'tstart'])
        df = df.groupby('info').apply(lambda g: merge_peak(g))
        df = df.groupby(['info', 'diff']).agg({'tstart': 'min', 'tend': 'max'}).reset_index()
        df = df.sort_values(by = 'info')
        df['label'], df['group'] = zip(*df['info'].apply(lambda g: split_info(g)))

        # annolist flag; only plot provided labels
        if annolist:
            df['index'] = df['label']
            df = df.set_index('index')
            df = df.loc[annolist.split(',')]

        self.df = pd.concat([self.df, df])




        # end of ANNO object

#############################################################################




def parseArgs():

    prs = argparse.ArgumentParser()
    prs.add_argument("--outputPDF",type=str, help="Name of the output PDF graphic", required = True)
    prs.add_argument("--ct", type=str, help="Base CT file to plot. By default, the first (0th) structure is plotted. Alternative structures can be selected by passing ,# after the CT file name (e.g. --ct ctfile.ct,3)")
    prs.add_argument("--fasta", type=str, help="Fasta sequence file")
    prs.add_argument("--refct", type=str, help="Reference CT to compare base ct against. By default the first (0th) structure is plotted. Alternative structures can be selected by passing ,# after the CT file name (e.g. --refct ctfile.ct,3)")
    prs.add_argument("--probability", type=str, help="Pairing probability file in dotplot format. By default, arcs are drawn for the following probability intervals: [0.03,0.1], [0.1,0.3], [0.3,0.8], [0.8,1.0]. These intervals can be modified by passing thresholds as comma-separated values. For example, --prob file.dp,0.03,0.1,0.3,0.8,1.0 specifies the default intervals. At most 4 intervals can be plotted, but fewer intervals are allowed (e.g. --prob file.dp,0.1,0.3,0.8,1.0 will plot 3 intervals).")

    prs.add_argument("--ringalpha", type=str, help="Plot alpha-scores from ringmapper correlation file. Default color thresholds are [2,5]. Can be modified by passing thresholds as comma-separated values (e.g. --corralpha corrs.txt,2,5)")
    
    prs.add_argument("--ringz", type=str, help="Plot Z-scores from ringmapper correlation file. Default color thresholds are [1,5]. Can be modified by passing thresholds as comma-separated values (e.g. --corrz corrs.txt,1,5)")
    
    prs.add_argument("--ringsig", type=str, help="Plot statistical significance from ringmapper file. Default color thresholds are [20,100]. Can be modified by passing thresholds as comma-separated values (e.g. --corrsig corrs.txt,20,500)")

    prs.add_argument("--catring",type=str, help='Plot N13, N7, and N137 rings from an input concatenated ring file.')

    prs.add_argument("--pairmap", type=str, help="Plot pairmap signals from pairmap file. By default plot principal & minor correlations. Can plot all complementary correlations by passing ,all (e.g. --pairmap pairmap.txt,all)")

    prs.add_argument("--compare_pairmap", type=str, help="Plot pairmap signals from second pairmap file. By default plot principal & minor correlations. Can plot all complementary correlations by passing ,all (e.g. --compare_pairmap pairmap.txt,all)")


    prs.add_argument("--ntshape", type=str, help="Color nucs by shape reactivty in provided shape/map file")  

    prs.add_argument("--profile", type=str, nargs = "+", help="Plot reactivity profile on top from shape/map file. If more than one profile is submitted via this flag, will plot mean normalized reactivity as well as stdev per nucleotide. ")
    prs.add_argument("--dmsprofile", type=str, nargs = "+", help='Normalize and plot DMS reactivity from profile file. If more than one profile is submitted via this flag, will plot mean normalized reactivity as well as stdev per nucleotide.')
    prs.add_argument("--N7profile",type=str, nargs = "+", help='Plots N7 reactivity from profile file. Note, you must specify a profile file with --profile or --dmsprofile in order for this command to work. If more than one profile is submitted via this flag, will plot mean normalized reactivity as well as stdev per nucleotide.')
    prs.add_argument("--bar_height", type=float, default=1.0, help="Scales N13 profile bar heights to a float in the range 0.0 < bar_height <= 1.00 so height *= bar_height. (eg if bar_height = .66 then the height will be scaled to 66 percent of the default value.) ")

    prs.add_argument("--bottom", action='store_true', help="Plot arcs on the bottom")

    prs.add_argument("--title", type=str, default='', help='Figure title')
    prs.add_argument("--showGrid", action="store_true", help="plot a grid on axis ticks")
    prs.add_argument("--bound", type=str, help="comma separated bounds of region to plot (e.g. --bound 511,796)")

    prs.add_argument("--filternc", action="store_true", help="filter out non-canonical pairs in ct")
    prs.add_argument("--filtersingle", action="store_true", help="filter out singleton pairs in ct")

    prs.add_argument("--contactfilter", type=int, help="filter rings by specified contact distance (int value)")
    
    prs.add_argument("--filternegcorrs", action="store_true", help='filter out negative correlations') 
    # annotation arguments
    prs.add_argument("--annotations", type=str, help="Path to a tab-delimited text file specifying features to plot.\
                                                     The first 2 columns contain the starting and ending coordinate. \
                                                     The third and fourth columns contain the feature label and group.")
    prs.add_argument("--name", type=str, default='query_RNA', help="RNA name. Default is query_RNA.")
    prs.add_argument("--coordinate", type=str, help='Path to a tab-delimited text file specifying RNA genomic coordinates.\
                                                Require 4 columns: chromosome start end strand. \
                                                Each line corresponds to one exon.')
    prs.add_argument("--annodir", type=str, help="Path to directory containing annotations files.")
    prs.add_argument("--cutoff", type=str, default = None, help='Cut off value for column 8 and 9. In eCLIP experiments, \
                    these columns contains p-value and enrichment scores (recommended threshold is 3 for both). Enter \
                    values as comma-separated list (e.g. 3,3).')
    prs.add_argument("--annolist", type=str, default=None, help='An ordered list of labels to plot separated by comma. \
                                                    (e.g. PUM1,PUM2 will only plot features with PUM1 and PUM2 labels.')
    prs.add_argument("--annocols", type=str, help='Define a color to plot annotation. If plotting two or more groups, \
                                                provide a list of colors separated by comma (e.g. darkblue,darkgreen)')



    args = prs.parse_args()
 

    if args.refct and not args.ct:
        exit("--refct is invalid without --ct")


    numplots = int(args.ct is not None)


    # subparse the prob argument
    args.probability_bins = None
    if args.probability:
        numplots += 1

        spl = args.probability.split(',')
        try:
            if len(spl) > 1:
                args.probability_bins = [float(x) for x in spl[1:]]
                args.probability = spl[0]
        except:
            raise TypeError('Incorrectly formatted --probability argument {}'.format(args.probability))
    

    # Make sure bar_height argument is valid
    if not 0.0 < args.bar_height <= 1.0:
        raise ValueError(" --bar_height argument must fall in range 0.0 < bar_height <= 1.0 ")
        

    def subparse3(arg, name):
        
        outarg = arg
        bins = None

        spl = arg.split(',')
        try:
            if len(spl) == 3:
                bins = [float(x) for x in spl[1:]]
                outarg = spl[0]
            elif len(spl) != 1:
                raise TypeError
        except:
            raise TypeError('Incorrectly formatted {0} argument {1}'.format(name, arg))

        return outarg, bins


    if args.ringz:
        numplots += 1
        args.ringz, args.ringz_bins = subparse3(args.ringz, '--ringz')
    
    if args.ringsig:
        numplots += 1
        args.ringsig, args.ringsig_bins = subparse3(args.ringsig, '--ringsig')

    if args.ringalpha:
        numplots += 1
        args.ringalpha, args.ringalpha_bins = subparse3(args.ringalpha, '--ringalpha')

    args.pairmap_all = False
    if args.pairmap:
        numplots += 1
        spl = args.pairmap.split(',')
        if len(spl)==1:
            pass
        elif len(spl)==2 and spl[1] == 'all':
            args.pairmap_all = True
            args.pairmap = spl[0]
        else:
            raise TypeError('Incorrectly formatted --pairmap argument {}'.format(args.pairmap))
    

    args.ringpairsuper = False
    if args.pairmap and (args.ringz or args.ringsig) and args.ct:
        args.ringpairsuper = True

    if numplots > 2 and not args.ringpairsuper:
        exit('Too many plots! Please select at maximum 2 of [--ct, --probability, --ringz, --ringsig, --pairmap --compare_pairmap]')


    if args.contactfilter and not args.ct:
        exit('Cannot perform contact filtering without --ct file')

    if args.ringpairsuper and args.contactfilter is None:
        args.contactfilter = 20

    # subparse the bounds argument
    if args.bound:
        args.bound = map(int, args.bound.split(','))
 
    
    # subparse the ct arguments
    if args.ct:
        spl = args.ct.split(',')
        if len(spl)==1:
            args.ctstructnum = 0
        else:
            args.ctstructnum = int(spl[1])
            args.ct = spl[0]

    if args.refct:
        spl = args.refct.split(',')
        if len(spl)==1:
            args.refctstructnum = 0
        else:
            args.refctstructnum = int(spl[1])
            args.refct = spl[0]


    return args





if __name__=="__main__":
 
    args = parseArgs()
    
    msg = None
    CT1=None
    
    aplot = ArcPlot(title = args.title, fasta=args.fasta)
    
    if args.annotations:
        anno = ANNO(name = args.name, length = aplot.length)
        anno.read_anno(args.annotations, args.bound)
    
    if args.coordinate and args.annodir:
        anno = ANNO(name = args.name, length = aplot.length)
        anno.get_rna_info(args.coordinate)
        anno.read_annodir(args.annodir, args.cutoff, args.annolist)
    
    # deal with bound 
    if args.bound:
        anno.apply_bound(args.bound)

    panel = 1
    if args.bottom:
        panel = -1

    if args.ct:
        
        CT1 = RNAtools.CT(args.ct, structNum=args.ctstructnum, filterNC=args.filternc, filterSingle=args.filtersingle)

        if args.refct:
            refCT = RNAtools.CT(args.refct, structNum=args.refctstructnum, filterNC=args.filternc, filterSingle=args.filtersingle)
            aplot.compareCTs( refCT, CT1, panel=panel)
        
        else:
            alpha = 0.7
            if args.ringpairsuper:
                alpha=0.2
            
            aplot.addCT( CT1, panel=panel, alpha=alpha)

        panel *= -1


    if args.probability:
        aplot.addPairProb( RNAtools.DotPlot(args.probability), panel=panel, bins=args.probability_bins)
        panel *= -1


    if args.pairmap:

        from pairmap_analysis import PairMap
        
        if args.ringpairsuper:
            aplot.addPairMap( PairMap(args.pairmap), panel=1, plotall=args.pairmap_all)

        else:
            aplot.addPairMap( PairMap(args.pairmap), panel=panel, plotall=args.pairmap_all)
            panel *= -1
        

    if args.ringz:
        aplot.addRings(args.ringz, panel=panel, metric='z', bins=args.ringz_bins,
                       contactfilter=(args.contactfilter, CT1), filterneg=args.filternegcorrs)
        panel *= -1
    
    if args.ringsig:
        aplot.addRings(args.ringsig, panel=panel, metric='sig', bins=args.ringsig_bins,
                       contactfilter=(args.contactfilter, CT1), filterneg=args.filternegcorrs)
        panel *= -1

    if args.ringalpha:
        aplot.addRings(args.ringalpha, panel=panel, metric='alpha', bins=args.ringalpha_bins,
                       contactfilter=(args.contactfilter, CT1), filterneg=args.filternegcorrs)
        panel *= -1

    if args.compare_pairmap:
        
        from pairmap_analysis import PairMap
        
        aplot.addPairMap( PairMap(args.compare_pairmap), panel=panel, plotall=args.pairmap_all)
        panel *= -1

    if args.annodir or args.annotations:
        aplot.addANNO(anno.df, length = anno.length, panel = -1, colors = args.annocols)
    if args.ntshape:
        aplot.colorSeqByMAP(args.ntshape)
    if args.profile:  
        if len(args.profile) > 1:
            avg_profile = average_profile(args.profile)
            nt_stdev = calc_stdev(args.profile)
            aplot.readProfile(avg_profile)
            aplot.add_error_bars(nt_stdev) 
        elif len(args.profile) == 1:
            aplot.readProfile(args.profile[0])
    if args.N7profile:
        if len(args.N7profile) > 1:
            avg_profile = average_profile(args.N7profile)
            nt_stdev = calc_stdev(args.N7profile)
            aplot.readN7Profile(avg_profile)
            aplot.add_error_bars(nt_stdev, lower = True) 
        elif len(args.N7profile) == 1:
            aplot.readN7Profile(args.N7profile[0])
    if args.dmsprofile:
        if len(args.dmsprofile) > 1:
            avg_profile = average_profile(args.dmsprofile)
            nt_stdev = calc_stdev(args.dmsprofile)
            aplot.readProfile(avg_profile, dms = True)
            #print("Finding std dev: ", nt_stdev)
            aplot.add_error_bars(nt_stdev) 
        elif len(args.dmsprofile) == 1:
            aplot.readProfile(args.dmsprofile[0], dms=True)

    if args.catring:
        aplot.addCatRings(args.catring, filterneg=args.filternegcorrs, metric='alpha', N7=True, N17=True, contactfilter=(args.contactfilter, CT1))
    if args.showGrid:
        aplot.grid = True

    aplot.writePlot( args.outputPDF, bounds = args.bound, msg=msg, bar_height_scale = args.bar_height) 

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
