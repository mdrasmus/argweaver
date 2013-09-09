# python libs
import sys

# rasmus libs
from rasmus import util


class TextCanvas:
    """Draw ascii art on a automatically growing matrix"""
    
    def __init__(self, default=' '):
        self.mat = util.Dict(dim=2, default=default)
        self.default = default
    
    
    def set(self, x, y, char):
        self.mat[int(y)][int(x)] = char
    
    
    def line(self, x1, y1, x2, y2, char='*'):
        # swap coords if needed
        if x1 > x2:
            x1, x2 = x2, x1
        if y1 > y2:
            y1, y2 = y2, y1
        
        nsamples = int(max(x2 - x1, y2 - y1, 1))
        dx = (x2 - x1) / float(nsamples)
        dy = (y2 - y1) / float(nsamples)
        
        for i in xrange(nsamples):
            self.set(x1 + i*dx, y1 + i*dy, char)
    
    
    def text(self, x, y, text, dir="horizontal", width=10000):
        x2 = 0
        y2 = 0
    
        if dir == "horizontal":
            for i in xrange(len(text)):
                if text[i] == "\n":
                    x2 = 0
                    y2 += 1
                elif x2 < width:
                    x2 += 1
                    self.set(x+x2, y+y2, text[i])
        elif dir == "vertical":
            for i in xrange(len(text)):
                if text[i] == "\n" or x2 > width:
                    y2 = 0
                    x2 += 1
                elif x2 < width:
                    y2 += 1
                    self.set(x+x2, y+y2, text[i])
        else:
            raise Exception("unknown text direction '%s'" % dir)

    
    def display(self, out=sys.stdout):
        ykeys = util.sort(self.mat.keys())
        
        y = min(ykeys)
        for ykey in ykeys:
            while y < ykey:
                y += 1
                out.write("\n")
            
            row = self.mat[ykey]
            xkeys = util.sort(row.keys())
            x = 0
            for xkey in xkeys:
                while x < xkey:
                    x += 1
                    out.write(self.default)
                out.write(row[xkey])
                x += 1
        out.write("\n")

        
