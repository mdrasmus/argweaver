"""
 file: plotting.py 
 authors: Matt Rasmussen
 date: 11/30/05

 Plotting classes and functions: R plotting, GNUPLOT wrapper, svg, heatmaps
"""

#from rasmus import util # will not allow from plotting import * inside util
from rasmus import svg


import sys, os, copy
import tempfile as temporaryfile


#=============================================================================
# Color maps


# common colors
red    = ( 1,  0,  0,  1)
orange = ( 1, .5,  0,  1)
yellow = ( 1,  1,  0,  1)
green  = ( 0,  1,  0,  1)
blue   = ( 0,  0,  1,  1)
purple = ( 1,  0,  1,  1)
black  = ( 0,  0,  0,  1)
grey   = (.5, .5, .5,  1)
white  = ( 1,  1,  1,  1)


class ColorMap:
    """ColorMap maps values on the real line to colors"""
    
    
    def __init__(self, table=[]):
        """
        'table' should be the following format:
        
        [
          [val1, color1],
          [val2, color2],
          [val3, color3],
          ...etc..
        ]
        
        Values bewteen val1 and val2 will be assigned a blend of 
        color1 and color2.  value-color pairs can be specified in any order 
        within table.
        
        """
        self.table = table
        
        self.table.sort(key=lambda x: x[0])
    
    
    def get(self, value):
        """Returns values in [0, 1]"""
    
        # determine where color falls in table
        for i in xrange(len(self.table)):
            if value <= self.table[i][0]:
                break
        if i > 0:
            i -= 1
        
        
        if value <= self.table[i][0]:
            # return lower bound color
            return self.table[i][1]
        elif value >= self.table[i+1][0]:
            # return upper bound color
            return self.table[i+1][1]
        else:
            # blend two nearest colors
            part = value - self.table[i][0]
            tot = float(self.table[i+1][0] - self.table[i][0])
            weight1 = (tot-part)/tot
            weight2 = part/tot
            
            newcolor = []
            color1 = self.table[i][1]
            color2 = self.table[i+1][1]
            for j in range(len(color1)):
                newcolor.append(weight1 * color1[j] + 
                                weight2 * color2[j])
            
            return newcolor
    
    
    def get_int(self, value):
        return [int(x*255) for x in self.get(value)]
    getInt = get_int
    

def get_webcolor(color, maxval=1):
    
    colstr = "#"
    for i in color:
        h = hex(int(i * 255.0 / maxval))[2:]
        if len(h) == 1:
            h = "0" + h
        colstr += h
    return colstr


def rainbow_color_map(data=None, low=None, high=None):
    if data != None:
        low = min(data)
        high = max(data)
    assert low != None and high != None
    
    return ColorMap([[low, blue],
                     [.5*low+.5*high, green],
                     [.25*low + .75*high, yellow],
                     [high, red]])
rainbowColorMap = rainbow_color_map


#=============================================================================
# svg plotting
    
def plothist2(x, y, ndivs1=20, ndivs2=20, width=500, height=500):
    from rasmus import util
    l, h = util.hist2(x, y, ndivs1, ndivs2)
    bwidth = util.bucket_size(x)
    bheight = util.bucket_size(y)
    
    #width *= bwidth/bheight
    
    heatmap(h, width/ndivs1, height/ndivs2)



def make_color_legend(filename, colormap, start, end, step, 
                    width=100, height=10):
    from rasmus import util
    s = svg.Svg(util.open_stream(filename, "w"))    
    s.beginSvg(width, height)
    
    xscale =  float(width) / (end + step - start)
    
    for i in util.frange(start, end + step, step):
        color = colormap.get(i)
        s.rect((i-start) * xscale, 
               0, 
               step*xscale, height, 
               color, color)
    
    s.endSvg()
makeColorLegend = make_color_legend




def heatmap(matrix, width=20, height=20, colormap=None, filename=None,
            rlabels=None, clabels=None, display=True, 
            xdir=1, ydir=1, 
            xmargin=0, ymargin=0,
            labelPadding=2,
            labelSpacing=4,
            mincutoff=None,
            maxcutoff=None,
            showVals=False,
            formatVals=str,
            valColor=black,
            clabelsAngle=270,
            clabelsPadding=None,
            rlabelsAngle=0,
            rlabelsPadding=None):

    from rasmus import util
    
    # determine filename
    if filename == None:
        filename = util.tempfile(".", "heatmap", ".svg")
        temp = True
    else:
        temp = False
    
    # determine colormap
    if colormap == None:
        colormap = rainbowColorMap(util.flatten(matrix))
    
    # determine matrix size and orientation
    nrows = len(matrix)
    ncols = len(matrix[0])
    
    if xdir == 1:
        xstart = xmargin
        ranchor = "end"
        coffset = width
    elif xdir == -1:
        xstart = xmargin + ncols * width
        ranchor = "start"
        coffset = 0
    else:
        raise Exception("xdir must be 1 or -1")
            
    if ydir == 1:
        ystart = ymargin
        roffset = height
        canchor = "start"
    elif ydir == -1:
        ystart = ymargin + nrows * width
        roffset = 0
        canchor = "end"
    else:
        raise Exception("ydir must be 1 or -1")
    
    
    # begin svg
    infile = util.open_stream(filename, "w")
    s = svg.Svg(infile)
    s.beginSvg(ncols*width + 2*xmargin, nrows*height + 2*ymargin)
    
    # draw matrix
    for i in xrange(nrows):
        for j in xrange(ncols):

            if mincutoff and matrix[i][j] < mincutoff: continue
            if maxcutoff and matrix[i][j] > maxcutoff: continue
            
            color = colormap.get(matrix[i][j])
            s.rect(xstart + xdir*j*width, 
                   ystart + ydir*i*height, 
                   xdir*width, ydir*height, color, color)
    
    # draw values
    if showVals:
        # find text size
        
        fontwidth = 7/11.0
        
        textsize = []
        for i in xrange(nrows):
            for j in xrange(ncols):

                if mincutoff and matrix[i][j] < mincutoff: continue
                if maxcutoff and matrix[i][j] > maxcutoff: continue
                
                strval = formatVals(matrix[i][j])
                if len(strval) > 0:
                    textsize.append(min(height,
                                        width/(float(len(strval)) * fontwidth)))
        textsize = min(textsize)


        yoffset = int(ydir == -1)
        for i in xrange(nrows):
            for j in xrange(ncols):

                if mincutoff and matrix[i][j] < mincutoff: continue
                if maxcutoff and matrix[i][j] > maxcutoff: continue
                
                strval = formatVals(matrix[i][j])
                s.text(strval, 
                       xstart + xdir*j*width, 
                       ystart + ydir*(i+yoffset)*height + 
                       height/2.0 + textsize/2.0, 
                       textsize,
                       fillColor=valColor)
    
    # draw labels
    if rlabels != None:
        assert len(rlabels) == nrows, \
            "number of row labels does not equal number of rows"

        if rlabelsPadding is None:
            rlabelsPadding = labelPadding
        
        for i in xrange(nrows):
            x = xstart - xdir*rlabelsPadding
            y = ystart + roffset + ydir*i*height - labelSpacing/2.
            s.text(rlabels[i], x, y, height-labelSpacing, anchor=ranchor,
                   angle=rlabelsAngle)
    
    if clabels != None:
        assert len(clabels) == ncols, \
            "number of col labels does not equal number of cols"

        if clabelsPadding is None:
            clabelsPadding = labelPadding
        
        for j in xrange(ncols):
            x = xstart + coffset + xdir*j*width - labelSpacing/2.
            y = ystart - ydir*clabelsPadding
            s.text(clabels[j], x, y, width-labelSpacing, anchor=canchor,
                   angle=clabelsAngle)
    
    # end svg
    s.endSvg()
    s.close()
    
    
    # display matrix
    if display:
        #if temp:
            os.system("display %s" % filename)
        #else:
        #    os.spawnl(os.P_NOWAIT, "display", "display", filename)
    
    # clean up temp files
    if temp:
        os.remove(filename)





       
    

