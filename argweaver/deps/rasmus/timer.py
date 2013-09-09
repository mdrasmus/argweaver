"""
    Timer class for timing nested sections of code
    
    file:   rasmus/timer.py
    author: Matt Rasmussen
    date:   2/4/2005

"""


# python libs
import os
import sys
import traceback
import time



# GLOBALS
_RASMUS_TIMER = None
_GLOBAL_NOTES = None


class Timer:
    def __init__(self, stream = sys.stderr, maxdepth=1e1000):
        self.reset()
        self.streams =  [(stream, maxdepth)]
        self.showErrors = True
        self.showWarnings = True
        self.quiets = 0
    

    def start(self, msg = ""):
        """Start a new timer"""
        
        if msg != "":
            self.indent()
            self._write("BEGIN %s:\n" % msg)
        self.msg.append(msg)
        self.flush()
        self.starts.append(time.time())
        
    
    def time(self):
        """Get the current duration of the timer"""
        
        return self.starts[-1] - time.clock()
    
    def stop(self):
        """Stop the last created timer and return duration in seconds"""
    
        duration = time.time() - self.starts.pop()
        msg = self.msg.pop()
        if msg != "":
            self.indent()
            
            if duration > 3600:
                pretty = "%.1fh" % (duration / 3600.)
            elif duration > 60:
                pretty = "%.1fm" % (duration / 60.)
            else:
                pretty = "%.3fs" % duration
                      
            if duration > .1:
                secs = "%.3fs" % duration
            else:
                secs = "%.3es" % duration
            
            self.write("END   %s: %s (%s)\n" % (msg, pretty, secs))
        self.flush()
        return duration
    
    def log(self, *text):
        """Write a message to the timer stream.  Message will be written with 
           current indentation level"""
        
        self.indent()
        for i in text:
            self._write("%s " % str(i))
        self._write("\n")
        self.flush()        
    
    def logExact(self, text):
        """Write the extact string 'text' to the timer output stream with no 
           additional indentation."""
        
        self._write(text)
        self.flush()
    
    def warn(self, text, offset=0):
        """Write a warning message to the timer output stream"""
        
        filename, lineno, func, code = traceback.extract_stack()[-2-offset]
        filename = os.path.basename(filename)
        
        if self.showWarnings:
            self.indent()
            self._write("WARNING: %s, line %d: %s\n" % (filename, lineno, text))
            self.flush()
        
    def error(self, text, offset=0):
        """Write an error message to the timer output stream"""
        
        filename, lineno, func, code = traceback.extract_stack()[-2-offset]
        filename = os.path.basename(filename)
        
        if self.showErrors:
            self.indent()
            self._write("ERROR: %s, line %d: %s\n" % (filename, lineno, text))
            self.flush()
    
    
    def indent(self):
        """Write the current indentation level to the timer output stream"""
        for i in range(self.depth()):
            self._write("  ")
    
    def reset(self):
        """Stop all timers"""
        self.msg = []
        self.starts = []
    
    def depth(self):
        """Get the current number of running timers"""
        return len(self.msg)
    
    def _write(self, text):
        """Private function for writing to output stream"""
        for stream, maxdepth in self.streams:
            if self.depth() < maxdepth and \
               self.quiets == 0:
                stream.write(text)
    
    def write(self, text):
        self._write(text)
        self.flush()
    
    def flush(self):
        for stream, maxdepth in self.streams:
            stream.flush()
    
    def addStream(self, stream, maxdepth=1e1000):
        self.streams.append((stream, maxdepth))
    
    def removeStream(self, stream):
        self.streams = filter(lambda x: x[0] != stream, self.streams)

    def suppress(self):
        """Calling this function will suppress timer output messages until 
           unsuppress() is called.  
           
           If suppress() is called multiple times,  unsuppress() must be called
           an equal number of times to resume timer  output.  This is useful for
           nesting suppress/unsuppress."""
        self.quiets += 1
    
    def unsuppress(self):
        """Calling this function will resume timer output messages that were
           disabled with suppress().
           
           If suppress() is called multiple times,  unsuppress() must be called
           an equal number of times to resume timer  output.  This is useful for
           nesting suppress/unsuppress."""
        self.quiets = max(self.quiets - 1, 0)


def globalTimer():
    global _RASMUS_TIMER
    if _RASMUS_TIMER == None:
        _RASMUS_TIMER = Timer()
    return _RASMUS_TIMER
    


def log(*text):
    return globalTimer().log(*text)

def logger(*text):
    return globalTimer().log(*text)
        
def logExact(text):
    return globalTimer().logExact(text)

def tic(msg = ""):
    return globalTimer().start(msg)

def toc():
    return globalTimer().stop()

def indent():
    return globalTimer().indent()

def warn(text, offset=0):
    return globalTimer().warn(text, offset+1)

def error(text, offset=0):
    return globalTimer().error(text, offset+1)



def note(*text):
    print >>notefile(), " ".join(text)

def noteflush():
    return notfile().flush()

def notefile(out = None):
    global _GLOBAL_NOTES

    if out == None:
        out = file("/dev/null", "w")
    if _GLOBAL_NOTES == None:
        _GLOBAL_NOTES = out
    return _GLOBAL_NOTES



################################################################################
# debugging info
#

def current_file(offset=0, abbrv=True):
    filename, lineno, func, code = traceback.extract_stack()[-2-offset]
    if abbrv:
        filename = os.path.basename(filename)
    return filename
    
def current_line(offset=0):
    filename, lineno, func, code = traceback.extract_stack()[-2-offset]
    return lineno

def current_func(offset=0):
    filename, lineno, func, code = traceback.extract_stack()[-2-offset]
    return func

def current_code(offset=0):
    filename, lineno, func, code = traceback.extract_stack()[-2-offset]
    return code



