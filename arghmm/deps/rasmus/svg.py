import os


def color2string(color):
    return "rgb(%d,%d,%d)" % (int(255 * color[0]),
                              int(255 * color[1]),
                              int(255 * color[2]))

def colorFields(strokeColor, fillColor):

    txt = ""

    if strokeColor:
        if isinstance(strokeColor, str):
            stroke_str = strokeColor
            stroke_op = None
        else:
            stroke_str = color2string(strokeColor)
            if len(strokeColor) > 3:
                stroke_op = strokeColor[3]
            else:
                stroke_op = None

        txt += "stroke='%s' " % stroke_str
        if stroke_op is not None:
            txt += "stroke-opacity='%f' " % stroke_op

    if fillColor:
        if isinstance(fillColor, str):
            fill_str = fillColor
            fill_op = None
        else:
            fill_str = color2string(fillColor)
            if len(fillColor) > 3:
                fill_op = fillColor[3]
            else:
                fill_op = None

        txt += "fill='%s' " % fill_str
        if fill_op is not None:
            txt += "fill-opacity='%f' " % fill_op    
    
    
    return txt

# common colors
#          r   g   b   a
red    = ( 1,  0,  0,  1)
orange = ( 1, .5,  0,  1)
yellow = ( 1,  1,  0,  1)
green  = ( 0,  1,  0,  1)
blue   = ( 0,  0,  1,  1)
purple = ( 1,  0,  1,  1)
black  = ( 0,  0,  0,  1)
grey   = (.5, .5, .5,  1)
white  = ( 1,  1,  1,  1)
null   = ( 0,  0,  0,  0)


class Svg:
    def __init__(self, stream):
        self.out = stream
    
    def close(self):
        self.out.close()
    
    
    def beginSvg(self, width, height):
        self.out.write(
            """<?xml version='1.0' encoding='UTF-8'?> 
            <!DOCTYPE svg PUBLIC '-//W3C//DTD SVG 1.1//EN' 
            'http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd'>\n""")
        self.out.write(
            """<svg width='%d' height='%d' 
            xmlns='http://www.w3.org/2000/svg' version='1.1'>\n""" % \
            (width, height))
            
        # default style
        self.out.write("<g>")
    
    def endSvg(self, close=True):
        self.out.write("</g></svg>")
        if close:
            self.close()

    def defLinearGradient(self, name, x1, y1, x2, y2,
                          stops=[], colors=[]):
        self.out.write('''<defs>
<linearGradient id="%s" x1="%f%%" y1="%d%%" x2="%f%%" y2="%f%%">''' %
                       (name, x1, y1, x2, y2))

        for stop, color in zip(stops, colors):
            if len(color) == 3:
                color = list(color) + [1.0]

            self.out.write('<stop offset="%f%%" style="stop-color:%s; '
                           'stop-opacity:%f"/>' % (stop,
                                                   color2string(color[:3]),
                                                   color[3]))

        self.out.write('</linearGradient></defs>')


    def defRadialGradient(self, name, cx=50, cy=50, r=50, fx=50, fy=50,
                          stops=[], colors=[]):
        self.out.write('''<defs>
<radialGradient id="%s" cx="%f%%" cy="%d%%" r="%f%%" fx="%f%%" fy="%f%%">''' %
                       (name, cx, cy, r, fx, fy))


        for stop, color in zip(stops, colors):
            if len(color) == 3:
                color = list(color) + [1.0]            
            self.out.write('<stop offset="%f%%" style="stop-color:%s; '
                           'stop-opacity:%f"/>' % (stop,
                                                   color2string(color[:3]),
                                                   color[3]))

        self.out.write('</radialGradient></defs>')
        
        
    def writeAttrOptions(self, color=None, **options):
        if color:
            if len(color) > 3:
                self.out.write("stroke-opacity='%f' stroke='%s' " % 
                (color[3], color2string(color)))
            else:
                self.out.write("stroke='%s' " % (color2string(color)))

        for key, val in options.iteritems():
            self.out.write("%s='%s' " % (key, val))
    
    
    def line(self, x1, y1, x2, y2, color=None, **options):
        self.out.write(
            """<line x1='%f' y1='%f' x2='%f' y2='%f' """ % 
            (x1, y1, x2, y2))
        
        if color != None:
            options["color"] = color
        self.writeAttrOptions(**options)
        
        self.out.write(" />\n")
            
    
    def polygon(self, verts, strokeColor=black, fillColor=black):
        self.out.write(
            "<polygon %s points='" % colorFields(strokeColor, fillColor))
        
        for i in xrange(0, len(verts), 2):    
            self.out.write("%f,%f " % (verts[i], verts[i+1]))
        self.out.write("' />\n")
    
    
    def rect(self, x, y, width, height, strokeColor=black, fillColor=black):
        self.out.write(
            """<rect x='%f' y='%f' width='%f' height='%f' %s />\n""" % \
            (x, y, width, height, colorFields(strokeColor, fillColor)))
    
    
    def circle(self, x, y, radius, strokeColor=black, fillColor=black):
        self.out.write("<circle cx='%f' cy='%f' r='%f' %s />\n" % \
            (x, y, radius, colorFields(strokeColor, fillColor)))
    
    def ellispe(self, x, y, xradius, yradius, strokeColor=black, fillColor=black):
        self.out.write("<ellipse  cx='%f' cy='%f' rx='%f' ry='%f' %s />\n" %\
            (x, y, xradius, yradius, colorFields(strokeColor, fillColor)))
    
    
    def text(self, msg, x, y, size, strokeColor=null, fillColor=black,
             anchor="start", angle=0):
        
        anglestr = "transform='translate(%f,%f) rotate(%f)'" % \
                    (x, y, angle)
        
        self.out.write(
            "<g %s><text x='0' y='0' font-size='%f' %s text-anchor='%s'>%s</text></g>\n" % \
            (anglestr, size, colorFields(strokeColor, fillColor), anchor, msg))
    
    
    def text2(self, msg, x, y, size, strokeColor=null, fillColor=black,
             anchor="start", angle=0):
        
        if angle != 0:
            anglestr = "" #transform='rotate(%f,0,0)'" % angle
        else:
            anglestr = ""
        
        
        self.out.write(
            "<text x='%f' y='%f' font-size='%f' %s text-anchor='%s' %s>%s</text>\n" % \
            (x, y, size, colorFields(strokeColor, fillColor), anchor, 
            anglestr, msg))

    
    def beginTransform(self, *options):  
        self.out.write("<g transform='")
        
        for option in options:                
            key = option[0]
            value = option[1:]
            
            if key == "scale":
                self.out.write("scale(%f, %f) " % value)
            
            elif key == "translate":
                self.out.write("translate(%f, %f) " % value)

            elif key == "rotate":
                self.out.write("rotate(%f, %f, %f) " % value)

                    
            else:
                raise Exception("unknown transform option '%s'" % key)
        
        self.out.write("' >\n")
    
    def endTransform(self):
        self.out.write("</g>\n")
    
    
    def beginStyle(self, style):
        self.out.write("<g style='%s'>\n" % style)
    
    def endStyle(self):
        self.out.write("</g>\n")
    

    def write(self, text):
        self.out.write(text)


    def comment(self, msg):
        self.out.write("\n<!-- %s -->\n\n" % msg)





def convert(filename, outfilename = None):
    if outfilename == None:
        outfilename = filename.replace(".svg", ".png")
    os.system("convert " +filename+ " " +outfilename)
    #os.system("rm " + filename)
    
    return outfilename



# testing
if __name__ == "__main__":
    svg = Svg(file("out.svg", "w"))
    
    svg.beginSvg(300, 500)
    
    svg.comment("MY COMMENT")
    
    svg.beginTransform(('scale', .5, .5))
    
    svg.line(0, 0, 100, 100, red)
    svg.rect(10, 10, 80, 100, black, (0, 1, 1, .5))
    svg.polygon([80,90, 100,100, 60,100], (0, 0, 0, 1), (0, 0, 1, .3))
    
    svg.endTransform()
    
    
    svg.beginStyle("stroke-width:3")
    svg.beginTransform(('translate', 200, 0))
    
    svg.line(0, 0, 100, 100, red)
    svg.rect(10, 10, 80, 100, black, (0, 1, 1, .5))
    svg.polygon([80,90, 100,100, 60,100], (0, 0, 0, 1), (0, 0, 1, .3))
    
    svg.endTransform()    
    svg.endStyle()
    
    svg.ellispe(150, 250, 70, 50, black, (.5, .5, .9, 1))
    svg.circle(150, 250, 50, black, red)
    svg.circle(150, 250, 30, white, blue)
    
    

    svg.beginStyle("font-family: arial")
    svg.beginTransform(('translate', 0, -200))
    svg.beginTransform(('translate', 0, 400))
    svg.text("A", 0, 0, 40, blue, (1, 0, 0, .1))
    svg.endTransform()
    
    svg.beginTransform(('translate', 0, 440))
    svg.text("C", 0, 0, 40, blue, (1, 0, 0, .1))
    svg.endTransform()
    svg.endTransform()
    svg.endStyle()
    
    svg.beginStyle("font-family: helvetica")
    svg.beginTransform(('translate', 0, -200))
    svg.beginTransform(('translate', 0, 480), ('scale', 1, 1))
    svg.text("T", 3, 0, 40, blue, (1, 0, 0, .1))
    svg.endTransform()
    
    svg.beginTransform(('translate', 0, 520), ('scale', 1, .1))
    svg.text("G", 0, 0, 40, blue, (1, 0, 0, .1))
    svg.endTransform()
    svg.endTransform()
    svg.endStyle()
    
    svg.line(35, 200, 35, 400, red)
    
    
    svg.beginStyle("font-family: courier")
    svg.text("* FIXED WIDTH    *", 100, 400, 10)
    svg.text("* IS THE DEFAULT *", 100, 410, 10)
    svg.endStyle()
    
    for i in range(0, 300, 10):
        color = (i / 300.0, 0, 1 - i/300.0, 1)
        svg.rect(i, 450, 10, 50, color, color)
    
    svg.endSvg()
    
    convert("out.svg")
