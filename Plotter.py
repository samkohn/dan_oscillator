#!/usr/bin/env python
#
# Useful plotting tools
#
# Note: This package requires ROOT and PyROOT to display graphs
#
# Created by: dadwyer@lbl.gov 2013/03/15
# Last Modified: dadwyer@lbl.gov 2013/03/16
# 

#from LinearAlgebra import array
from array import array as py_array
from types import ListType

try:
    __import__("ROOT")
except ImportError:
    print "ERROR: You must have ROOT and PyROOT installed to display ROOT plots"
    print "Sorry"

from ROOT import TCanvas, TGraph, TGraphErrors, TH1F, TH2F, gStyle, kBlue, TColor, TLegend
gStyle.SetOptStat(0)
gStyle.SetLineColor(kBlue)

plotterColors = [1,2,4,8,9,46]

class Palette(object):
    "Simple color palette class"
    def __init__(self, size, steps, r, g, b):
        "Constructor"
        self._size = size
        self._colors = py_array('i',[0]*size)
        palIdx = TColor.CreateGradientColorTable(len(steps),
                                                 py_array('d',steps),
                                                 py_array('d',r),
                                                 py_array('d',g),
                                                 py_array('d',b),
                                                 size)
        for i in range(size):
            self._colors[i] = palIdx + i
        return
    def enable(self):
        "Enable this palette"
        gStyle.SetPalette(self._size, self._colors)

ContourLevels = 255
Palette_DarkBody = Palette(ContourLevels, 
                           [0.0, 0.05, 0.20, 0.95, 1.0],
                           [1.0, 1.0, 1.0,  0.5, 0.0],
                           [1.0, 1.0, 0.55, 0.0, 0.0],
                           [1.0, 0.0, 0.0,  0.0, 0.0])

Palette_Oscillation = Palette(ContourLevels, 
                              [0.0, 0.05, 0.20, 0.95, 1.0],
                              reversed([1.0, 1.0, 1.0,  0.5, 0.0]),
                              reversed([1.0, 1.0, 0.55, 0.0, 0.0]),
                              reversed([1.0, 0.0, 0.0,  0.0, 0.0]))

Palette_Symmetric = Palette(ContourLevels, 
                            [0.0, 0.5, 1.0],
                            [0.0, 1.0, 1.0],
                            [0.0, 1.0, 0.0],
                            [1.0, 1.0, 0.0])

Palette_SymmetricHighContrast = Palette(ContourLevels, 
                                        [0.0, 0.48, 0.5, 0.52, 1.0],
                                        [0.0, 0.5,  1.0, 1.0,  1.0],
                                        [0.0, 0.5,  1.0, 0.5,  0.0],
                                        [1.0, 1.0,  1.0, 0.5,  0.0])


gStyle.SetNumberContours(255)
#Palette_DarkBody.enable()
Palette_Symmetric.enable()

class Plotter(object):
    "Simple class for plotting graphs and histograms"
    def __init__(self):
        "Constructor"
        self._canvases = []
        self._graphs = []
        self._histograms = []
        self._keepers = []
        self._currentCanvas = None

    def makeCanvas(self):
        "Make a new canvas"
        canvas = TCanvas()
        self._currentCanvas = canvas
        self._canvases.append(canvas)
        return canvas

    def makeGraph(self, xData, yData, xDataErr=None, yDataErr=None):
        "Make a new graph from data"
        graph = None
        if xDataErr or yDataErr:
            # Has errors
            if xDataErr: 
                xDataErr = py_array('f',xDataErr)                
            else:
                xDataErr = py_array('f',[0]*len(xData))
            if yDataErr: 
                yDataErr = py_array('f',yDataErr)                
            else:
                yDataErr = py_array('f',[0]*len(yData))
            graph = TGraphErrors(len(xData), 
                                 py_array('f',xData),
                                 py_array('f',yData),
                                 xDataErr,
                                 yDataErr)
        else:
            graph = TGraph(len(xData), 
                           py_array('f',xData),
                           py_array('f',yData))
        graph.SetLineWidth(2)
        self._graphs.append( graph )
        return graph

    def makeHistogram1D(self, xBinEdges, data):
        "Make a 1-D histogram from data and bins"
        nBinsX = len(xBinEdges)-1
        hist = TH1F("","",nBinsX, py_array('f',xBinEdges))
        for binIdx in range(1,nBinsX+1):
            hist.SetBinContent( binIdx, data[binIdx-1] )
        self._histograms.append( hist )
        return hist

    def makeHistogram2D(self, xBinEdges, yBinEdges, data):
        "Make a 2-D histogram from data and bins"
        nBinsX = len(xBinEdges)-1
        nBinsY = len(yBinEdges)-1
        hist = TH2F("","",nBinsX, py_array('f',xBinEdges),
                    nBinsY, py_array('f',yBinEdges))
        for binIdxX in range(1,nBinsX+1):
            for binIdxY in range(1,nBinsY+1):
                if type(data) is ListType:
                    hist.SetBinContent( binIdxX, binIdxY, 
                                        data[binIdxX-1][binIdxY-1] )
                else:
                    hist.SetBinContent( binIdxX, binIdxY, 
                                        data[binIdxX-1,binIdxY-1] )
        self._histograms.append( hist )
        return hist

    def keep(self, item):
        "Keep reference to item to avoid auto-deletion"
        self._keepers.append(item)

gPlotter = Plotter()

def drawGraphs(graphs, canvas=None, padIdx=None, labels=None):
    "Draw multiple 1-D graphs into canvas"
    if canvas:
        if padIdx:
            canvas.cd(padIdx)
        else:
            canvas.cd()
    curColorIdx = 0
    firstGraph = True
    for graph in graphs:
        graph.SetLineColor( plotterColors[curColorIdx] )
        if firstGraph:
            graph.Draw("APL")
            firstGraph = False
        else:
            graph.Draw("PLsame")
        curColorIdx += 1
        if curColorIdx >= len(plotterColors):
            curColorIdx = 0
    if labels:
        # Make a legend
        leg = TLegend(0.6,0.6,0.88,0.88)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        for gIdx in range(len(graphs)):
            leg.AddEntry(graphs[gIdx],labels[gIdx],"lp")
        leg.Draw()
        gPlotter.keep(leg)
    return

def draw2DGraphs(graphs, canvas=None):
    "Draw multiple 2-D graphs into canvas"
    if canvas: 
        canvas.cd()
    else:
        canvas = pltr.makeCanvas()
    canvas.Divide(2,3)
    canvasIdx = 1
    for graph in graphs:
        canvas.cd(canvasIdx)
        graph.Draw("colz")
        canvasIdx += 1
    return

