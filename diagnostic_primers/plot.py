#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plot.py

Code to generate plots of pdp output

(c) The James Hutton Institute 2018

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2018 The James Hutton Institute
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os

import pandas as pd
import plotly
import plotly.graph_objs as go


def markerscatter(fname, outdir):
    """Plot distance information for amplicons

    - fname        Path to the input data file
    - outdir       Path to the output directory for writing the plot

    The plot is a scatterplot of mean pairwise distance vs number of
    non-unique amplicons. This gives a pareto curve on two selection
    parameters analogous to a ROC curve in that the 'best' selection
    is expected to be towards the top left of the plot.
    """
    data = pd.read_csv(fname, sep="\t")
    markers = go.Scatter(
        x=data.nonunique, y=data.dist_mean, text=data.primer, mode="markers"
    )
    layout = go.Layout(
        xaxis=dict(
            range=(-0.1, 1.2 * max(data.unique)),
            showgrid=True,
            zeroline=True,
            showline=False,
            autotick=True,
            ticks="",
            showticklabels=True,
        ),
        yaxis=dict(
            range=(-0.01, 1.2 * max(data.dist_mean)),
            showgrid=True,
            zeroline=True,
            showline=False,
            autotick=True,
            ticks="",
            showticklabels=True,
        ),
        hovermode="closest",
    )
    data = [markers]
    outfname = os.path.join(outdir, "markerscatter.html")
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outfname, validate=False, auto_open=False)
