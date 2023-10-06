import matplotlib.pyplot as plt
import numpy as np
import os
import sys

import padeopsio.budgetkey as budgetkey


class PlotIO(): 
    """
    Interface class for plotting and saving figures using a BudgetIO object. 
    """

    keylab = budgetkey.key_labels()

    x_lab = '$x$'
    y_lab = '$y$'
    z_lab = '$z$'
    fs = 12
    plt.rcParams.update({'font.size': fs})

    cbar = "coolwarm"
    dpi = 100
    # add other constants here? Note: these are mutable by the user

    def set_fontsize(fontsize): 
        """
        Sets the font size
        """
        PlotIO.fs = fontsize
        plt.rcParams.update({'font.size': fontsize})


    def plot_xy(io, budget_terms, z, xlim=None, ylim=None, 
                ax=None, xlabel=None, ylabel=None): 
        """
        Plots a slice in the xy-plane. 
        
        If `budget_terms` is a list, then plot_xy will recursively call itself for each term individually
        and produce figures of each separately. 
        
        Arguments
        ---------
        io (BudgetIO obj) : BudgetIO object linked to padeops data
        budget_terms (list or str) : terms to plot
        z (float) : z-location to pull slices from
        xlim, ylim : slice bounds, see BudgetIO.slice()
        
        Returns 
        -------
        
        """
        slices = io.slice(budget_terms, xlim=xlim, ylim=ylim, zlim=z)

        for term in budget_terms: 
            fig, ax = plt.subplots()
            im = ax.imshow(slices[term].T, extent=slices['extent'], origin='lower')

            common_cbar(fig, im, ax, label=PlotIO.keylab[term])

            ax.set_xlabel(PlotIO.x_lab)
            ax.set_ylabel(PlotIO.y_lab)

            plt.show()

    
    def xz_slice(): 
        """
        Plots a slice in the xz-plane. 
        """
        pass

    def yz_slice(): 
        """
        Plots a slice in the yz-plane. 
        """
        pass

 
# ----------- additional helper functions ------------

def common_cbar(fig, image, ax=None, location='right', label=None, height=1., width=0.02): 
    """
    Helper function to add a colorbar
    
    colorbar: https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph 
 
    parameters
    ----------
    fig : figure object to add a colorbar to
    image : AxesImage object or similar that has a color scale
    ax : (optional) Axes object to reference for size of colorbar
    location : string, # TODO
    height : height of colorbar as a fraction of the given axis's height (or selected, if no axis given). 
        Default 1
    width : width of colorbar as a fraction of the give axis's height
        Default 0.02
    """
    if ax is None:  # need to pick an axis from the figure...
        # ax = fig.axes[0]  
        # TODO: make this better... maybe choose the largest one? Maybe make a new one each time? ugh
        ax = common_axis(fig)
    
    h = ax.get_position().height  # axis height
    cax = fig.add_axes([ax.get_position().x1+h*width, 
                        ax.get_position().y0+h/2*(1-height), 
                        h*width, 
                        height*h])  # [start x, start y, width, height]
    cbar = fig.colorbar(image, cax=cax)
    
    if label is not None: 
        cbar.set_label(label)
    
    return cbar  # this probably will never be used


def common_axis(fig, xlabel=None, ylabel=None, title=None): 
    """
    Helper function to format common axes

    format common axes: 
    https://stackoverflow.com/questions/16150819/common-xlabel-ylabel-for-matplotlib-subplots
    """
    # create a new axis spanning the whole figure
    ax_new = fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    return ax_new  # return the newly created axis

