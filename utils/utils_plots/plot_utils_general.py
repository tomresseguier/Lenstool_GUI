import numpy as np
from matplotlib.colors import hsv_to_rgb

#import sys
#import os
#module_dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(module_dir)



def make_palette(hue_range, sat_range, alpha=1) :
    hues = np.linspace(0, 1-1/hue_range, hue_range)
    sats = np.linspace(1, 0, sat_range)
    colors = []
    for sat in sats :
        for hue in hues :
            colors.append( np.append(hsv_to_rgb([hue, sat, 1.]), alpha) )
    return colors


def make_palette_dict(alpha=1) :
    colors_keys = ['rgb1', 'rgb2', 'rgb3'] + ['cmy1', 'cmy2', 'cmy3']
    #colors_keys_circle = ['r1', 'g1', 'b1', 'r2', 'g2', 'b2', 'r3', 'g3', 'b3'] + ['c1', 'm1', 'y1', 'c2', 'm2', 'y2', 'c3', 'm3', 'y3']
    colors = {}
    for key in colors_keys :
        colors[key] = [ [0,0,0,alpha] for i in range(3) ]
    for i, values in enumerate([ [0,0.5], [0,1], [0.5,1] ]) :
        for j in range(3) :
            colors[ colors_keys[i] ][j][j] = values[1]
            colors[ colors_keys[i+3] ][j][j] = values[0]
            indexes = [0,1,2]
            indexes.remove(j)
            for index in indexes :
                colors[ colors_keys[i] ][j][index] = values[0]
                colors[ colors_keys[i+3] ][j][index] = values[1]
    return colors




