from starter2 import *

def four_on_the_floor():
    fig,ax=plt.subplots(2,2,figsize=(12,12))
    fig.subplots_adjust(wspace=0, hspace=0)
    ax[0][0].xaxis.tick_top()
    ax[0][1].xaxis.tick_top()
    ax[0][1].xaxis.set_label_position('top')
    ax[0][0].xaxis.set_label_position('top')
    ax[1][1].yaxis.tick_right()
    ax[0][1].yaxis.tick_right()
    ax[1][1].yaxis.set_label_position('right')
    ax[0][1].yaxis.set_label_position('right')

    return fig, ax


def three_way_bean(figsize=(8,8), left=None,width=None,bottom=None,height=None,histdepth=None):
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    if left is None: left = 0.15
    if width is None: width = 0.62
    if bottom is None: bottom = 0.11
    if height is None: height = 0.62
    if histdepth is None: histdepth = 0.2
    #left, width = 0.15, 0.62
    #.bottom, height = 0.11, 0.62
    bottom_h = left_h = left + width + histdepth

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    fig=plt.figure(1)
    fig.set_size_inches(figsize)

    axScatter = plt.axes(rect_scatter)
    axHistx =   plt.axes(rect_histx)
    axHisty =   plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    return fig, axScatter,axHistx, axHisty


