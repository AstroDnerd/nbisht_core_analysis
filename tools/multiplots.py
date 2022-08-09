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


def three_way_bean():
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    fig=plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx =   plt.axes(rect_histx)
    axHisty =   plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    return fig, axScatter,axHistx, axHisty


