/*
  Mashome - Genome Browser Mashup library
  Matt Rasmussen 2012

  This library allows one to visualize *very* custom genome browser
  tracks within a host genome browser web page (e.g. UCSC, Ensembl,
  generic).  Custom tracks are are built and controlled using
  javascript.  Data displayed within tracks can be fetched using AJAX.
  
*/

(function (window) {

//===================================================================
// misc functions

// import javascript script url
function importScript(url, onLoad) {
    var e = document.createElement("script");
    e.setAttribute('type', 'text/javascript');
    e.setAttribute('src', url);

    if (onLoad) {
        if (e.attachEvent)
            e.attachEvent("onload", onLoad);
        else
            e.addEventListener("load", onLoad, false);
    }
    
    document.body.appendChild(e);        
}

// import multiple javascript script urls
function importScripts(urls, onLoad, serial) {
    var nwait = urls.length;
    
    if (nwait == 0) {
        if (onLoad)
            onLoad();
        return;
    }

    if (typeof serial == "undefined")
        serial = false;

    if (serial) {
        // serial loading
        var i = 0;
        function onLoadScript() {
            i++;
            if (i < nwait)
                importScript(urls[i], onLoadScript);
            else if (onLoad)
                onLoad();
        }
        importScript(urls[i], onLoadScript);

    } else {
        // parallel loading
        function onLoadScript() {
            nwait--;
            if (nwait == 0 && onLoad)
                onLoad();
        }
    
        for (var i in urls)
            importScript(urls[i], onLoadScript);
    }
}


//===================================================================
// cookies
    
function setCookie(c_name, value, exdays)
{
    var exdate = new Date();
    exdate.setDate(exdate.getDate() + exdays);
    var c_value = escape(value) + ((exdays==null) ? "" : 
                                       "; expires=" + exdate.toUTCString());
    document.cookie = c_name + "=" + c_value;
}

function getCookie(c_name)
{
    var i, x, y, ARRcookies = document.cookie.split(";");
    for (i=0; i<ARRcookies.length; i++) {
        x = ARRcookies[i].substr(0, ARRcookies[i].indexOf("="));
        y = ARRcookies[i].substr(ARRcookies[i].indexOf("=") + 1);
        x = x.replace(/^\s+|\s+$/g, "");
            if (x == c_name)
                return unescape(y);
    }
}

function deleteCookie(c_name)
{
    setCookie(c_name, "", -1);
}


//=================================================================
// mashome

// constants
var updateInterval = 500;
var urlCookie = "mashome-url";
var cookieTime = 30; // a month

// check whether mashone already exists
if (window.mashome) {
    window.mashome.init();
    return window.mashome;
}


// setup mashome config
if (!window.mashomeConfig)
    window.mashomeConfig = {};

    
//=====================================================================
// Tracks
    
// base Track class
// options:
//   name   -- track display name
//   height -- track height
function Track(options) {

    this.mashome = null;
    this.elm =  $("<div></div>");
    this.elm.css({"overflow": "hidden",
                  "width": "100%"});
    this.width = 0;
    this.sidebarWidth = 0;
    this.mainWidth = 0;

    if (typeof options == "undefined")
        options = {};
    if (typeof options.name == "undefined")
        options.name = "user track";
    if (typeof options.height == "undefined")
        options.height = 100;
    
    this.name = options.name;
    this.height = options.height;
        
        
    // create default track layout
    this.create = function () {
        this.elm.empty();
        
        this.width = this.mashome.width;
        this.sidebarWidth = this.mashome.sidebarWidth;
        this.mainWidth = this.mashome.mainWidth;
        
        this.sidebar = $("<div></div>");
        this.sidebar.text(this.name);
        this.sidebar.css({"width": this.sidebarWidth -  1,
                          "background-color": "#fff",
                          "border-right": "solid 1px #fcc",
                          "height": this.height,
                          "float": "left"});
        this.main = $("<div></div>");
        this.main.css({"width": this.mainWidth,
                       "background-color": "#fff",
                       "height": this.height,
                       "float": "left"});
            
        this.elm.append(this.sidebar);
        this.elm.append(this.main);
    };
    
    // change the height of a track
    this.setHeight = function (height) {
        this.height = height;
        if (this.main)
            this.main.css("height", height);
        if (this.sidebar)
            this.sidebar.css("height", height);
    };

    //-------------------------------------------
    // callbacks

    // callback for when track is added to mashome
    this.onAddTrack = function(view) {};
        
    // callback for when genomic position changes
    this.onViewChange = function(view) {};  
};
   

// basic ruler track
// options:
//   name    -- track display name
//   height  -- track height
//   select  -- callback for select event
function RulerTrack(options) {
    // call super
    this.Track = Track;
    this.Track(options);
    
    // callback for when a region is clicked
    this.onSelect = options.select;
    
        
    // display the current position
    this.showPosition = function (pos, hover) {
        this.posText.text(this.mashome.positionToString(pos+1));
        if (hover)
            this.main.css("color", "#800");
        else
            this.main.css("color", "#000");
    };

    // change cursor position to position 'base'
    // base -- a chromosome position
    this.setCursor = function(base) {
        if (typeof base == "undefined")
            base = this.selectedPos;
        else
            this.selectedPos = base;
        
        var x = this.mashome.positionToPixel(base, this.view);
        this.cursor.css("left", x);
    };

    //-------------------------------------------------
    // mashome events

    // callback for when track is added to mashome
    this.onAddTrack = function (view) {
        var that = this;
        this.view = view;
        this.selectedPos = Math.round((view.end + view.start) / 2);
        var fontSize = 12

        this.main.css({"border-bottom": "1px solid #555",
                       "border-top": "1px solid #555",
                       "background-color": "#eee", 
                       "position": "relative",
                       "overflow": "hidden",
                       "font-size": fontSize + "px"});
        this.main.click(function (e) {that.onClick(e)});
        this.main.mousemove(function (e) {that.onMouseMove(e)});
        this.main.mouseout(function (e) {that.onMouseOut(e)});

        this.posText = $("<div></div>");
        this.posText.css({"text-align": "center",
                          "position": "relative",
                          "top": (this.height - 2 - fontSize) / 2});
        this.main.append(this.posText);

        this.cursor = $("<div></div>");
        this.cursor.css({"width": "0px",
                         "height": this.height - 2,
                         "border": "1px solid #f00",
                         "position": "absolute",
                         "top": 0,
                         "left": this.mainWidth / 2});
        this.main.append(this.cursor);
        this.setCursor();
    };
    
    // callback for when genomic position changes
    this.onViewChange = function (view) {
        this.view = view;
        this.showPosition(this.selectedPos);
        this.setCursor();
    };
    
    //----------------------------------------------
    // mouse events
        
    this.onClick = function (event) {
        var x = event.pageX - this.main.position().left;
        var p = this.mashome.pixelToPosition(x, this.view);
        
        this.showPosition(p);
        this.setCursor(p);
        this.onSelect(p);
    };

    this.onMouseMove = function (event) {
        var x = event.pageX - this.main.position().left;
        var p = this.mashome.pixelToPosition(x, this.view);
        this.showPosition(p, true);
    };
        
    this.onMouseOut = function (event) {
        this.showPosition(this.selectedPos);
    };
};
RulerTrack.prototype = new Track;


// A track containing a single pre-transformed canvas
// options:
//   name    -- track display name
//   height  -- track height
function CanvasTrack(options) {
    // call super
    this.Track = Track;
    this.Track(options);
        
    this.onAddTrack = function (view) {
        this.canvas = $("<canvas></canvas>");
        this.canvas.attr("width", this.mainWidth);
        this.canvas.attr("height", this.height);
        this.ctx = this.canvas.get(0).getContext("2d");

        this.main.append(this.canvas);
    };
        
    this.beginTransform = function(view) {
        var c = this.ctx;
        var scale = this.mainWidth / (view.end - view.start);

        c.clearRect(0, 0, this.mainWidth, this.height);
        c.save();
        c.scale(scale, 1);
        c.translate(-view.start, 0);
        c.lineWidth /= scale;            
    };
    
    this.endTransform = function() {
        this.ctx.restore();
    }
}
CanvasTrack.prototype = new Track;
    


//=====================================================================
// Host genome browsers

// example view object
// view = {chrom: "chr1", start: 1, end: 100};

function HostBrowser () {
    
    // initializes HostBrowser information
    this.init = function () {
        this.tracks = null;
        this.width = 400;
        this.sidebarWidth = 0;
        this.mainWidth = 400;

        return false;
    };

    // returns true if host browser is detected
    this.isDetected = function () { return false; }

    // get the current genome view
    this.getView = function () {};

    // set the genome view
    this.setView = function (view) {};

    // returns true if genome view has changed since last call to getView()
    this.hasViewChanged = function () { return false; };
};


// Interface for UCSC Browser
function UCSCBrowser () {
    this.HostBrowser = HostBrowser;
    this.HostBrowser();

    var imgTableId = "#imgTbl";
    var lastPostext = null;
    var lastView = null;
    
    // initializes HostBrowser information
    this.init = function () {
        // get UCSC variables
        try {
            this.tracks = $(imgTableId);
            this.width = this.tracks.width();
            this.sidebarWidth = window.hgTracks.insideX;
            this.mainWidth = this.width - this.sidebarWidth;
        } catch (err) {
            return false;
        }
        return true;
    };

    // returns true if host browser is detected
    this.isDetected = function () {
        return ($(imgTableId).length > 0 &&
                window.hgTracks &&
                window.genomePos);
    };
    
    // returns current genome view
    this.getView = function () {
        var postext = window.genomePos.get();
        if (postext == lastPostext)
            return lastView;
        
        // cache view
        lastPostext = postext;
        lastView = window.parsePosition(postext);

        // convert to 0-index
        lastView.start -= 1;
        return lastView;
    };

    // returns true if genome view has changed since last call to getView()
    this.hasViewChanged = function () {
        var postext = window.genomePos.get();
        return postext != lastPostext;
    };

};
UCSCBrowser.prototype = new HostBrowser;


// Interface for generic stub browser
function GenericBrowser () {
    this.HostBrowser = HostBrowser;
    this.HostBrowser();

    var browserId = "#browser";
    var sidebarId = "#sidebar";
    var viewId = "#view";
    var lastPostext = null;
    var lastView = {"chrom": "chr", "start": 0, "end": 1000};
    var browser = null;
    
    // initializes HostBrowser information
    this.init = function () {
        // get browser variables
        try {
            this.tracks = $(browserId);
            this.width = this.tracks.width();
            this.sidebarWidth = $(sidebarId).width();
            this.mainWidth = this.width - this.sidebarWidth;

            browser = window.browser;
        } catch (err) {
            return false;
        }
        return true;
    };

    // returns true if host browser is detected
    this.isDetected = function () {
        return $(browserId).length > 0 && window.browser;
    };
    
    // returns current genome view
    this.getView = function () {
        var postext = browser.viewtext;
        if (postext == lastPostext)
            return lastView;
        
        // cache view
        lastPostext = postext;
        lastView = browser.parsePosition(postext);
        return lastView;
    };

    // set the genome view
    this.setView = function (view) {
        browser.goToView(view);
    };

    // returns true if genome view has changed since last call to getView()
    this.hasViewChanged = function () {
        return browser.viewtext != lastPostext;
    };
};
GenericBrowser.prototype = new HostBrowser;


    
//=====================================================================
// mashome

// mashome toolbar
function Toolbar(urlCookie) {
    var that = this;
    
    this.urlCookie = urlCookie;
    this.onUrlChange = null;
    
    this.create = function() {
        this.elm = $(
            "<div id='mashome-toolbar'>" +
            "<form id='mashome-toolbar-form'>" +
            "tracks URL: <input id='mashome-url' type='text'></input>" +
            "<span id='mashome-close' style='float:right'><a style='color: black' href='javascript: mashome.close()'>close</a></span>" +
            "</form></div>");
        this.elm.css({"background-color": "#eee",
                      "border-top": "1px solid #ccc",
                      "padding": "1px",
                      "padding-left": "20px",
                      "text-align": "left"});
        this.elm.find("#mashome-url").css("width", 400);
        
        this.urlInput = this.elm.find("#mashome-url");
        this.elm.find("#mashome-toolbar-form").submit(function(e){
                that.onSubmit(e)});

        this.initUrl();

        this.closeElm = this.elm.find("#mashome-close");
        this.closeElm.css('padding', '5px');
        
        return this.elm;
    }

    this.initUrl = function () {
        var url = getCookie(this.urlCookie);
        if (url) {
            this.urlInput.val(url);
        }
    }
    
    this.onSubmit = function (event) {
        event.preventDefault();

        var url = this.urlInput.val();
        setCookie(this.urlCookie, url, cookieTime);

        // perform callback
        if (this.onUrlChange)
            this.onUrlChange(url);
    }

    this.getUrl = function () {
        return this.urlInput.val();
    }

    this.setUrl = function (url) {
        this.urlInput.val(url);
    }
};


// main mashome object
var mashome = {
    init: function () {

        // mashome variables
        this.tracks = [];
        this.config = window.mashomeConfig;
        this.tracksScript = "";

        // auto-detect host genome browser
        this.defaultBrowsers = [UCSCBrowser, GenericBrowser];
        if (!this.config.browsers)
            this.config.browsers = [];
        for (var i in this.defaultBrowsers)
            this.config.browsers.push(this.defaultBrowsers[i]);

        // configure genome host browser
        this.hostBrowser = this.autoDetectHostBrowser();
        if (!this.hostBrowser)
            return;
        if (!this.configHostBrowser())
            return;
            
        // setup gui
        if (!this.elm) {
            if (this.config.tracksPosition == "before") 
                this.hostTracks.before(this.create());
            else 
                // default
                this.hostTracks.after(this.create());
        }

        // check for initial url
        if (this.config.tracksUrl)
            this.toolbar.setUrl(this.config.tracksUrl);
        var url = this.toolbar.getUrl();
        if (url)
            this.loadTrackScript(url);
        
        // ensure loop is active
        this.startLoop();
    },

    // remove mashome UI from host browser
    close: function () {
        this.clearTracks();
        this.elm.remove();
        delete window.mashome;
    },

    // auto-detect which host browser is present
    autoDetectHostBrowser: function (browsers) {
        if (!browsers)
            browsers = this.config.browsers;

        for (var i in browsers) {
            var browser = new browsers[i];
            if (browser.isDetected())
                return browser;
        }
        
        alert("NONE");

        return null;
    },

    // configures mashome for host browser
    configHostBrowser: function () {
        if (!this.hostBrowser.init())
            return false;
        
        this.hostTracks = this.hostBrowser.tracks;
        this.width = this.hostBrowser.width;
        this.sidebarWidth = this.hostBrowser.sidebarWidth;
        this.mainWidth = this.hostBrowser.mainWidth;

        return true;
    },

    // create GUI
    create: function () {
        var that = this;
            
        this.elm = $("<div id='mashome'></div>");
        this.elm.css({"border": "solid 1px black",
                      "width": this.width,
                      "font-size": "12px"});

        // create track element
        this.tracksElm = $("<div id='mashome-track'></div>");
        this.tracksElm.css({"background-color": "white"});
        this.elm.append(this.tracksElm);
            
        // create toolbar
        this.toolbar = new Toolbar(urlCookie);
        this.toolbar.onUrlChange = function (url) { 
            that.loadTrackScript(url); };
        this.elm.append(this.toolbar.create());
            
        return this.elm;
    },

    //============================================================
    // main loop
        
    startLoop: function() {
        var that = this;
        if (!this.updateId)
            this.updateId = setInterval(function() {that.update();}, 
                                        updateInterval);
    },

    stopLoop: function() {
        clearInterval(this.updateId);
        this.updateId = null;
    },

    update: function() {
        // check for view change
        if (this.hasViewChanged())
            this.onViewChange(this.getView());
    },

    //===========================================================
    // tracks

    addTrack: function (track) {
        // init track variables
        track.mashome = this;
        
        // add track
        this.tracks.push(track);
        
        // add element
        track.create();
        this.tracksElm.append(track.elm);
        
        // notify track
        var view = this.getView();
        track.onAddTrack(view);
        track.onViewChange(view);
    },

    removeTrack: function (track) {
        for (var i in this.tracks) {
            if (this.tracks[i] == track) {
                this.tracks.splice(i, 1);
                return;
            }
        }
    },

    clearTracks: function () {
        this.tracks = [];
        this.tracksElm.empty();
    },

    loadTrackScript: function (url, onLoad) {
        this.clearTracks();
        this.tracksScript = url;
        importScript(url, onLoad);
    },

    reloadTrackScript: function (onLoad) {
        this.loadTrackScript(this.tracksScript, onLoad);
    },
    

    // builtin track classes
    Track: Track,
    RulerTrack: RulerTrack,
    CanvasTrack: CanvasTrack,

    //===========================================================
    // callbacks
    
    onViewChange: function (view) {
        // notify all tracks of view change
        for (var i in this.tracks)
            this.tracks[i].onViewChange(view);
    },
    

    //============================================================
    // view and position functions
    
    getView: function () {
        return this.hostBrowser.getView();
    },

    hasViewChanged: function () {
        return this.hostBrowser.hasViewChanged();
    },
        
    pixelToPosition: function (x, view) {
        if (!view)
            view = this.getView();
        var basePerPixel = (view.end - view.start) / this.mainWidth;
        return view.start + Math.round(basePerPixel * x);
    },

    positionToPixel: function (pos, view) {
        if (!view)
            view = this.getView();
        var pixelPerBase = this.mainWidth / (view.end - view.start);
        return Math.round(pixelPerBase * (pos - view.start));
    },

    pad: function(str, char, len) {
        while (str.length < len)
            str = char + str;
        return str;
    },

    positionToString: function(pos) {
        var s = "";
        while (pos > 1000) {
            s = "," + this.pad("" + (pos % 1000), "0", 3) + s;
            pos = Math.floor(pos / 1000);
        }
        s = pos + s;
        return s;
    },

    viewToString: function(view) {
        return view.chrom + ":" + this.positionToString(view.start+1) + "-" + 
        this.positionToString(view.end);
    },
    
    
    //==========================================================
    // misc functions
    
    importScript: importScript,
    importScripts: importScripts,
    getCookie: getCookie,
    setCookie: setCookie,

};

// load dependencies and start mashome module
var scripts = [];
if (!window.jQuery)
    scripts.push('//ajax.googleapis.com/ajax/libs/jquery/1.7.2/jquery.min.js');

importScripts(scripts, function () {
    if (!window.mashomeConfig.noautoload)
        $(document).ready(function () {mashome.init() });
});

window.mashome = mashome;
return mashome;
})(window);

