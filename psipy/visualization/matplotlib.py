from matplotlib.projections.polar import ThetaFormatter
import numpy as np


def clear_axes_labels(ax):
    """
    Remove labels from both x and y axes.
    """
    ax.set_xlabel('')
    ax.set_ylabel('')


def set_theta_formatters(ax):
    """
    Set both x and y axes to have theta formatters (ie. degrees)
    """
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ThetaFormatter())
        axis.set_minor_formatter(ThetaFormatter())


def setup_polar_ax(ax):
    if ax is None:
        ax = plt.gca(projection='polar')
    if ax.name != 'polar':
        raise ValueError('ax must have a polar projection')
    return ax

