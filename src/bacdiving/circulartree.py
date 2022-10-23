#!/usr/bin/python3


import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib
from matplotlib.lines import Line2D

# The following code is copied (with small adaptations) from https://github.com/koonimaru/radialtree
# Copyright of koonimaru

colormap_list = ["nipy_spectral", "terrain", "gist_rainbow", "CMRmap", "coolwarm", "gnuplot", "gist_stern", "brg",
                 "rainbow"]


def ct_plot(Z2, fontsize=8, figsize=None, pallete="gist_rainbow", addlabels=True, show=True, sample_classes=None,
         colorlabels=None,
         colorlabels_legend=None, output_dir = "./"):

    if figsize == None and colorlabels != None:
        figsize = [10, 5]
    elif figsize == None and sample_classes != None:
        figsize = [10, 5]
    elif figsize == None:
        figsize = [5, 5]
    linewidth = 0.5
    R = 1
    width = R * 0.1
    space = R * 0.05
    if colorlabels != None:
        offset = width * len(colorlabels) / R + space * (len(colorlabels) - 1) / R + 0.05
        #print(offset)
    elif sample_classes != None:
        offset = width * len(sample_classes) / R + space * (len(sample_classes) - 1) / R + 0.05
        #print(offset)
    else:
        offset = 0
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['svg.fonttype'] = 'none'
    xmax = np.amax(Z2['icoord'])
    ymax = np.amax(Z2['dcoord'])

    ucolors = sorted(set(Z2["color_list"]))
    cmp = cm.get_cmap(pallete, len(ucolors))
    if type(cmp) == matplotlib.colors.LinearSegmentedColormap:
        cmap = cmp(np.linspace(0, 1, len(ucolors)))
    else:
        cmap = cmp.colors
    fig, ax = plt.subplots(figsize=figsize)
    i = 0
    label_coords = []
    for x, y, c in sorted(zip(Z2['icoord'], Z2['dcoord'], Z2["color_list"])):
        _color = cmap[ucolors.index(c)]
        if c == "C0":
            _color = "black"

        # transforming original x coordinates into relative circumference positions and y into radius
        # the rightmost leaf is going to [1, 0]
        r = R * (1 - np.array(y) / ymax)
        _x = np.cos(2 * np.pi * np.array(
            [x[0], x[2]]) / xmax)  # transforming original x coordinates into x circumference positions
        _xr0 = _x[0] * r[0]
        _xr1 = _x[0] * r[1]
        _xr2 = _x[1] * r[2]
        _xr3 = _x[1] * r[3]
        _y = np.sin(2 * np.pi * np.array(
            [x[0], x[2]]) / xmax)  # transforming original x coordinates into y circumference positions
        _yr0 = _y[0] * r[0]
        _yr1 = _y[0] * r[1]
        _yr2 = _y[1] * r[2]
        _yr3 = _y[1] * r[3]
        # plt.scatter([_xr0, _xr1, _xr2, _xr3],[_yr0, _yr1, _yr2,_yr3], c="b")

        # if y[0]>0 and y[3]>0:
        # _color="black"
        # plotting radial lines
        plt.plot([_xr0, _xr1], [_yr0, _yr1], c=_color, linewidth=linewidth)
        plt.plot([_xr2, _xr3], [_yr2, _yr3], c=_color, linewidth=linewidth)

        # plotting circular links between nodes
        if _yr1 > 0 and _yr2 > 0:
            link = np.sqrt(r[1] ** 2 - np.linspace(_xr1, _xr2, 100) ** 2)
            plt.plot(np.linspace(_xr1, _xr2, 100), link, c=_color, linewidth=linewidth)
        elif _yr1 < 0 and _yr2 < 0:
            link = -np.sqrt(r[1] ** 2 - np.linspace(_xr1, _xr2, 100) ** 2)

            plt.plot(np.linspace(_xr1, _xr2, 100), link, c=_color, linewidth=linewidth)
        elif _yr1 > 0 and _yr2 < 0:
            _r = r[1]
            if _xr1 < 0 or _xr2 < 0:
                _r = -_r
            link = np.sqrt(r[1] ** 2 - np.linspace(_xr1, _r, 100) ** 2)
            plt.plot(np.linspace(_xr1, _r, 100), link, c=_color, linewidth=linewidth)
            link = -np.sqrt(r[1] ** 2 - np.linspace(_r, _xr2, 100) ** 2)
            plt.plot(np.linspace(_r, _xr2, 100), link, c=_color, linewidth=linewidth)

        # Calculating the x, y coordinates and rotation angles of labels
        if y[0] == 0:
            label_coords.append([(1.05 + offset) * _xr0, (1.05 + offset) * _yr0, 360 * x[0] / xmax])
            # plt.text(1.05*_xr0, 1.05*_yr0, Z2['ivl'][i],{'va': 'center'},rotation_mode='anchor', rotation=360*x[0]/xmax)
            i += 1
        if y[3] == 0:
            label_coords.append([(1.05 + offset) * _xr3, (1.05 + offset) * _yr3, 360 * x[2] / xmax])
            # plt.text(1.05*_xr3, 1.05*_yr3, Z2['ivl'][i],{'va': 'center'},rotation_mode='anchor', rotation=360*x[2]/xmax)
            i += 1

    if addlabels == True:
        assert len(Z2['ivl']) == len(label_coords), "Internal error, label numbers " + str(
            len(Z2['ivl'])) + " and " + str(len(label_coords)) + " must be equal!"

        # Adding labels
        for (_x, _y, _rot), label in zip(label_coords, Z2['ivl']):
            plt.text(_x, _y, label, {'va': 'center'}, rotation_mode='anchor', rotation=_rot, fontsize=fontsize)

    if colorlabels != None:
        assert len(Z2['ivl']) == len(label_coords), "Internal error, label numbers " + str(
            len(Z2['ivl'])) + " and " + str(len(label_coords)) + " must be equal!"

        j = 0
        outerrad = R * 1.05 + width * len(colorlabels) + space * (len(colorlabels) - 1)
        #print(outerrad)
        intervals = []
        for i in range(len(label_coords)):
            _xl, _yl, _rotl = label_coords[i - 1]
            _x, _y, _rot = label_coords[i]
            if i == len(label_coords) - 1:
                _xr, _yr, _rotr = label_coords[0]
            else:
                _xr, _yr, _rotr = label_coords[i + 1]
            d = ((_xr - _xl) ** 2 + (_yr - _yl) ** 2) ** 0.5
            intervals.append(d)
        colorpos = intervals
        labelnames = []
        for labelname, colorlist in colorlabels.items():
            colorlist = np.array(colorlist)[Z2['leaves']]
            outerrad = outerrad - width * j - space * j
            innerrad = outerrad - width
            patches, texts = plt.pie(colorpos, colors=colorlist,
                                     radius=outerrad,
                                     counterclock=True,
                                     startangle=label_coords[0][2] * 0.5)
            circle = plt.Circle((0, 0), innerrad, fc='whitesmoke')
            plt.gca().add_patch(circle)
            labelnames.append(labelname)
            j += 1

        if colorlabels_legend != None:
            for i, labelname in enumerate(labelnames):
                #print(colorlabels_legend[labelname]["colors"])
                colorlines = []
                for c in colorlabels_legend[labelname]["colors"]:
                    colorlines.append(Line2D([0], [0], color=c, lw=4))
                leg = plt.legend(colorlines,
                                 colorlabels_legend[labelname]["labels"],
                                 bbox_to_anchor=(1.5 + 0.3 * i, 1.0),
                                 title=labelname)
                plt.gca().add_artist(leg)
    elif sample_classes != None:
        assert len(Z2['ivl']) == len(label_coords), "Internal error, label numbers " + str(
            len(Z2['ivl'])) + " and " + str(len(label_coords)) + " must be equal!"

        j = 0
        outerrad = R * 1.05 + width * len(sample_classes) + space * (len(sample_classes) - 1)
        #print(outerrad)
        intervals = []
        for i in range(len(label_coords)):
            _xl, _yl, _rotl = label_coords[i - 1]
            _x, _y, _rot = label_coords[i]
            if i == len(label_coords) - 1:
                _xr, _yr, _rotr = label_coords[0]
            else:
                _xr, _yr, _rotr = label_coords[i + 1]
            d = ((_xr - _xl) ** 2 + (_yr - _yl) ** 2) ** 0.5
            intervals.append(d)
        colorpos = intervals
        labelnames = []
        colorlabels_legend = {}
        for labelname, colorlist in sample_classes.items():
            ucolors = sorted(list(np.unique(colorlist)))
            type_num = len(ucolors)
            _cmp = cm.get_cmap(colormap_list[j], type_num)
            _colorlist = [_cmp(ucolors.index(c) / (type_num - 1)) for c in colorlist]
            _colorlist = np.array(_colorlist)[Z2['leaves']]
            outerrad = outerrad - width * j - space * j
            innerrad = outerrad - width
            patches, texts = plt.pie(colorpos, colors=_colorlist,
                                     radius=outerrad,
                                     counterclock=True,
                                     startangle=label_coords[0][2] * 0.5)
            circle = plt.Circle((0, 0), innerrad, fc='whitesmoke')
            plt.gca().add_patch(circle)
            labelnames.append(labelname)
            colorlabels_legend[labelname] = {}
            colorlabels_legend[labelname]["colors"] = _cmp(np.linspace(0, 1, type_num))
            colorlabels_legend[labelname]["labels"] = ucolors
            j += 1

        if colorlabels_legend != None:
            for i, labelname in enumerate(labelnames):
                #print(colorlabels_legend[labelname]["colors"])
                colorlines = []
                for c in colorlabels_legend[labelname]["colors"]:
                    colorlines.append(Line2D([0], [0], color=c, lw=4))
                leg = plt.legend(colorlines,
                                 colorlabels_legend[labelname]["labels"],
                                 bbox_to_anchor=(1.5 + 0.3 * i, 1.0),
                                 title=labelname)
                plt.gca().add_artist(leg)

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.spines.left.set_visible(False)
    ax.spines.bottom.set_visible(False)
    plt.xticks([])
    plt.yticks([])
    if colorlabels != None:
        maxr = R * 1.05 + width * len(colorlabels) + space * (len(colorlabels) - 1)
    elif sample_classes != None:
        maxr = R * 1.05 + width * len(sample_classes) + space * (len(sample_classes) - 1)
    else:
        maxr = R * 1.05
    plt.xlim(-maxr, maxr)
    plt.ylim(-maxr, maxr)
    if show == True:
        plt.show()
    else:
        plt.savefig(output_dir + 'overview_treeplot.pdf', bbox_inches='tight')