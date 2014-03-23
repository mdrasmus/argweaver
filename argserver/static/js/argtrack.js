(function () {


//==================================================================
// logging

// log a json object to the error console
function logJSON(json)
{
    console.log(JSON.stringify(json));
}

//=============================================================================
// Misc functions

function min(a, b) {
    return (a < b) ? a : b;
}
function max(a, b) {
    return (a > b) ? a : b;
}

function formatView(view, oneIndex) {
    var start = view.start;
    if (oneIndex)
        start += 1
    return view.chrom+':'+start+'-'+view.end;
}

//=================================================================
// sites functions


function getHighFreq(col) {
    var counts = {};

    for (var i=0; i<col.length; i++) {
        if (!(col[i] in counts))
            counts[col[i]] = 0;
        counts[col[i]]++;
    }

    var maxcount = 0;
    var maxchar;
    for (var c in counts) {
        if (counts[c] > maxcount) {
            maxcount = counts[c];
            maxchar = c;
        }
    }

    return maxchar;
}


function isParsimonySite(col, tree)
{
    var dna2int = {"A": 0, "C": 1, "G": 2, "T": 3};
    var maxcost = 1000000;

    function walk(name) {
        var node = tree.nodes[name];
        var costs = new Array(4);

        if (node.children.length == 0) {
            var i = parseInt(name);
            for (var a=0; a<4; a++)
                costs[a] = maxcost;
            costs[dna2int[col[i]]] = 0;
        } else {
            var left_costs = walk(node.children[0]);
            var right_costs = walk(node.children[1]);

            for (var a=0; a<4; a++) {
                var left_min = maxcost;
                var right_min = maxcost;
                for (var b=0; b<4; b++) {
                    left_min = min(left_min, Number(a != b) + left_costs[b]);
                    right_min = min(right_min, Number(a != b) + right_costs[b]);
                }
                costs[a] = left_min + right_min;
            }
        }
        return costs;
    }
    var costs = walk(tree.root);

    var root_min = maxcost;
    for (var a=0; a<4; a++)
        root_min = min(root_min, costs[a]);

    return root_min <= 1;
}

//=============================================================================
// tree functions

function getTreeAges(tree) {
    var ages = {};

    function walk(name) {
        var node = tree.nodes[name];

        if (node.children.length == 0)
            ages[name] = 0.0;
        else {
            for (var i in node.children)
                walk(node.children[i]);

            var child = tree.nodes[node.children[0]];
            ages[name] = ages[node.children[0]] + child.dist;
        }
    }
    walk(tree.root);

    return ages;
}


function getTreeLeaves(tree, name) {
    // get leaf order
    var order = [];
    function walk(name) {
        var node = tree.nodes[name];
        if (node.children.length == 0)
            order.push(name);
        for (var i in node.children)
            walk(node.children[i]);
    }
    if (typeof name == "undefined")
        name = tree.root;
    walk(name);
    return order;
}


function layoutTree(tree, options)
{
    // setup options
    if (typeof options == "undefined")
        options = {};
    if (typeof options.x == "undefined")
        options.x = 0;
    if (typeof options.y == "undefined")
        options.y = 0;
    if (typeof options.yscale == "undefined")
            options.yscale = 1.0;
    var minlen = 0;
    var maxlen = 1e12;

    var sizes = {};
    var nodept = {}; // distance between node x and left bracket x
    var layout = {};


    function walk(name) {
        var node = tree.nodes[name];

        // compute node sizes
        sizes[name] = 0;
        for (var i in node.children) {
            var cname = node.children[i];
            sizes[name] += walk(cname);
        }

        if (node.children == 0) {
            sizes[name] = 1;
            nodept[name] = 1;
        } else {
            var n = node.children.length;
            var top = nodept[node.children[0]];
            var bot = (sizes[name] - sizes[node.children[n-1]]) +
                nodept[node.children[n-1]];
            nodept[name] = (top + bot) / 2.0;
        }

        return sizes[name];
    }
    walk(tree.root);

    // determine x, y coordinates
    function walk2(name, x, y) {
        var node = tree.nodes[name];
        var ychildren = y + min(max(node.dist, minlen), maxlen) *
            options.yscale;
        layout[name] = [x + nodept[name], ychildren];

        if (!node.children.length == 0) {
            var xchild = x;
            for (var i in node.children) {
                var cname = node.children[i];
                walk2(cname, xchild, ychildren);
                xchild += sizes[cname];
            }
        }
    }
    walk2(tree.root, options.x, options.y);

    return layout;
}


//=============================================================================
// tree and spr drawing

function drawTree(tree, layout, labels)
{
    var ages = getTreeAges(tree);
    if (typeof layout == "undefined")
        layout = layoutTree(tree, {x:0, y:ages[tree.root], yscale: -1.0});

    var g = group();

    function walk(name) {
        var size = 1;
        var node = tree.nodes[name];

        var x = layout[name][0];
        var y = layout[name][1];

        if (node.parent) {
            var px = layout[node.parent][0];
            var py = layout[node.parent][1];
            g.push(lines(x, y, x, py));
        }

        if (node.children.length > 0) {
            var n = node.children.length;
            var x1 = layout[node.children[0]][0];
            var x2 = layout[node.children[n-1]][0];

            g.push(lines(x1, y, x2, y));

            for (var i in node.children)
                walk(node.children[i]);
        }
    }
    walk(tree.root);

    // draw labels
    var w = 10;
    var h = 12;
    for (var name in tree.nodes) {
        if (tree.nodes[name].children.length == 0) {
            var x = layout[name][0];
            var y = layout[name][1];
            var label = (labels && name in labels ? labels[name] : ""+name);
            g.push(translate(x, y, zoomClamp({},
                    textLabel(label, -w/2, -h, w/2, 0))));
        }
    }

    return g;
}


function drawSpr(spr, layout)
{
    var rx = layout[spr.recomb_node][0];
    var ry = spr.recomb_time;
    var cx = layout[spr.coal_node][0];
    var cy = spr.coal_time;
    var w = 5;

    return group(style({fill: "#f00"}),
                 translate(rx, ry, zoomClamp({},
                     quads(-w, -w, -w, w,
                           w, w, w, -w))),
                 style({stroke: "#c00"}),
                 lines(rx, ry, cx, cy));
}


//=====================================================================
// create tracks

function setup() {

    // setup argtrack url
    var argtrackurl = mashome.getCookie("argtrack-url") ||
        "http://localhost:8080";

    // create ruler track
    var ruler = new mashome.RulerTrack({name: "ruler",
                                        height: 20,
                                        select: function (pos) {
            treetrack.showTree(ruler.view.chrom, pos);
            sitetrack.showSites(ruler.view.chrom, pos); }});
    mashome.addTrack(ruler);


    // plot recombinations
    var recombs = new mashome.CanvasTrack({name: "recombs", height: 20});
    recombs.onViewChange = function (view) {
        var that = this;
        this.view = view;

        // limit range
        if (this.view.end - this.view.start > 4e6) {
            this.ctx.clearRect(0, 0, this.mainWidth, this.height);
            return;
        }

        $.ajax({
            dataType: 'jsonp',
            url: argtrackurl + '/sprs/' + formatView(view)
        }).done(function (result) {
            that.plot(JSON.parse(result));
        });
    };
    recombs.plot = function (sprs) {
        var c = this.ctx;
        this.beginTransform(this.view);

        c.strokeStyle = "#f00";

        c.beginPath();
        for (var i in sprs) {
            var spr = sprs[i];
            var x = spr.pos;
            c.moveTo(x, 0);
            c.lineTo(x, this.height);
        }
        c.stroke();
        c.closePath();

        this.endTransform();
    };
    mashome.addTrack(recombs);


    // plot mutations
    var muts = new mashome.CanvasTrack({name: "sites", height: 50});
    muts.onViewChange = function (view) {
        var that = this;
        this.view = view;
        this.fullSites = true;

        // limit range
        if (this.view.end - this.view.start > 4e6) {
            this.ctx.clearRect(0, 0, this.mainWidth, this.height);
            return;
        }

        $.ajax({
            dataType: 'jsonp',
            url: argtrackurl + '/sites/' + formatView(view)
        }).done(function (result) {
            that.plot(JSON.parse(result));
        });
    };
    muts.plot = function (sites) {
        var c = this.ctx;
        this.beginTransform(this.view);

        if (sites.length == 0)
            return;
        var nseqs = sites[0].col.length;
        var scale = this.height / nseqs;

        c.strokeStyle = "#00f";

        var hscale = this.mainWidth / (this.view.end - this.view.start);
        if (hscale > 1.0)
            c.lineWidth *= hscale;

        c.beginPath();
        if (this.fullSites) {
            // draw full site pattern
            for (var i in sites) {
                var site = sites[i];
                var highFreq = getHighFreq(site.col);
                var x = site.pos;

                for (var j=0; j<site.col.length; j++) {
                    if (site.col[j] != highFreq) {
                        c.moveTo(x, j * scale);
                        c.lineTo(x, (j+1) * scale);
                    }
                }
            }
        } else {
            // just indicate polymorphic site
            for (var i in sites) {
                var x = sites[i].pos;
                c.moveTo(x, 0);
                c.lineTo(x, this.height);
            }
        }
        c.stroke();
        c.closePath();

        this.endTransform();
    };
    mashome.addTrack(muts);


    // plot reordered mutations
    var muts2 = new mashome.CanvasTrack({name: "reordered sites", height: 50});
    muts2.onViewChange = function (view) {
        var that = this;
        this.view = view;
        this.main.css("border-top", "1px solid #ccc");
        this.main.css("border-bottom", "1px solid #ccc");

        // limit range
        if (this.view.end - this.view.start > 1e6) {
            this.ctx.clearRect(0, 0, this.mainWidth, this.height);
            return;
        }

        var region = formatView(view);
        var sites;
        $.ajax({
            dataType: 'jsonp',
            url: argtrackurl + '/sites/' + region
        }).then(function (result) {
            sites = JSON.parse(result);
            return $.ajax({
                dataType: 'jsonp',
                url: argtrackurl + '/trees/'+ region
            });
        }).then(function (result) {
            var trees = JSON.parse(result);
            that.plot(trees, sites);
        });
    };
    muts2.plot = function (trees, sites) {
        var c = this.ctx;
        this.beginTransform(this.view);

        // filter trees
        var trees2 = [];
        for (var i in trees)
            if (trees[i].tag == "TREE")
                trees2.push(trees[i]);
        trees = trees2;


        if (sites.length == 0 || trees.length == 0)
            return;
        var nseqs = sites[0].col.length;
        var scale = this.height / nseqs;
        var treei = 0;

        var tree = parseNewick(trees[treei].tree);
        var order = getTreeLeaves(tree);
        for (var i in order)
            order[i] = parseInt(order[i]);

        c.strokeStyle = "#00f";

        c.beginPath();
        // draw full site pattern
        for (var i in sites) {
            var site = sites[i];
            var highFreq = getHighFreq(site.col);
            var x = site.pos;
            while (treei < trees.length && site.pos > trees[treei].end) {
                treei++;
                tree = parseNewick(trees[treei].tree);
                order = getTreeLeaves(tree);
                for (var j in order)
                    order[j] = parseInt(order[j]);
            }

            if (isParsimonySite(site.col, tree)) {
                if (c.strokeStyle != "#00f") {
                    c.stroke();
                    c.closePath();
                    c.beginPath();
                }
                c.strokeStyle = "#00f";
            } else {
                if (c.strokeStyle != "#f00") {
                    c.stroke();
                    c.closePath();
                    c.beginPath();
                }
                c.strokeStyle = "#f00";
            }

            for (var j=0; j<site.col.length; j++) {
                var k = order[j];
                if (site.col[k] != highFreq) {
                    c.moveTo(x, j * scale);
                    c.lineTo(x, (j+1) * scale);
                }
            }
        }
        c.stroke();
        c.closePath();

        this.endTransform();
    };
    //mashome.addTrack(muts2);


    // create arg track
    var argtrack = new mashome.Track({name: "ARG", height: 200});
    argtrack.onAddTrack = function (view) {
        var that = this;
        this.main.css({"border-top": "1px solid #ccc",
                       "border-bottom": "1px solid #ccc"});

        this.text = $("<div></div>");
        this.main.append(this.text);

        this.canvas = $('<canvas></canvas>');
        this.canvas.attr("width", this.mainWidth-2);
        this.canvas.attr("height", this.height-2);
        this.main.append(this.canvas);

        this.view = null;
        this.scanvas = new Summon.Canvas(this.canvas.get(0));
        this.firstDisplay = true;
        this.labels = null;

        //var updateInterval = 500;
        //this.updateId = null;
        //if (!this.updateId)
        //    this.updateId = setInterval(function() {that.update();}, 
        //                                updateInterval);
    };
    /*argtrack.onRemoveTrack = function (view) {
        clearInterval(this.updateId);
        this.updateId = null;
    }
    argtrack.update = function() {
        if (this.view) {
            var visible = this.scanvas.getVisible();
            this.view.start = visible[0];
            this.view.end = visible[2];
        }
    };*/
    argtrack.onViewChange = function (view) {
        // don't redraw if view hasn't changed
        this.view = view;
        this.show(view);
    };
    argtrack.show = function(view) {
        var that = this;
        var layout;
        var region = formatView(view);
        $.ajax({
            dataType: 'jsonp',
            url: argtrackurl + '/arg-layout/' + region
        }).then(function (resultLayout) {
            layout = JSON.parse(resultLayout);
            that.draw(layout);
            return $.ajax({
                dataType: 'jsonp',
                url: argtrackurl + '/sites/'+ region
            });
        }).then(function (resultSites) {
            var sites = JSON.parse(resultSites);
            that.drawSites(layout, sites);
            that.scanvas.draw();
        });
    }
    argtrack.draw = function(layout) {
        this.scanvas.clear();
        this.scanvas.add(this.drawARG(layout));
        this.scanvas.draw();

        if (this.firstDisplay) {
            this.scanvas.home('exact');
        } else {
            var vis = this.scanvas.getVisible();
            this.scanvas.setVisible(this.view.start, vis[1], 
                                    this.view.end, vis[3], 'exact');
        }
        this.firstDisplay = false;
    }
    argtrack.drawARG = function(layout) {
        var gap = 5;
        var g = group(style({stroke: "rgba(0, 0, 0, .5)"}));
        
        var lastLeafLayout = null;
        for (var i=0; i<layout.length; i++) {
            var block = layout[i].block;
            var leafLayout = layout[i].leafLayout;
            var x1 = block[1] - 1 + gap;
            var x2 = block[2];
            var l = lines();

            for (var name in leafLayout) {
                if (leafLayout.hasOwnProperty(name)) {
                    var y = leafLayout[name];
                    l.extend([x1, y, x2, y]);
                    if (lastLeafLayout) {
                        var y0 = lastLeafLayout[name];
                        l.extend([x1 - gap, y0, x1, y]);
                    }
                }
            }
            g.push(l);
            lastLeafLayout = leafLayout;
        }
        
        return g;
    }
    argtrack.drawSites = function(layout, sites) {
        for (var i in sites) {
            var site = sites[i];
            var highFreq = getHighFreq(site.col);

            for (var j in site.col) {
                var a = site.col[j];
            }
        }
    }
    mashome.addTrack(argtrack);


    // create arg track
    var treetrack = new mashome.Track({name: "local tree", height: 300});
    treetrack.onAddTrack = function (view) {
        this.main.css({"border-top": "1px solid #ccc",
                       "border-bottom": "1px solid #ccc"});

        this.text = $("<div></div>");
        this.main.append(this.text);

        this.canvas = $('<canvas></canvas>');
        this.canvas.attr("width", this.mainWidth-2);
        this.canvas.attr("height", this.height-2);
        this.main.append(this.canvas);

        this.scanvas = new Summon.Canvas(this.canvas.get(0));
        this.firstDisplay = true;
        this.labels = null;

        // show initial tree
        var mid = Math.round((view.end + view.start) / 2);
        this.showTree(view.chrom, mid);
    };
    treetrack.onViewChange = function (view) {
    };
    treetrack.showTree = function(chrom, pos) {
        var that = this;
        $.ajax({
            dataType: 'jsonp',
            url: argtrackurl + '/tree-spr/'+chrom+':'+pos
        }).done(function (result) {
            that.displayTreeSpr(JSON.parse(result))
        });
    }
    treetrack.displayTree = function(treeInfo) {
        var tree = parseNewick(treeInfo.tree);

        this.scanvas.clear();
        this.scanvas.add(drawTree(tree));
        this.scanvas.draw();

        if (this.firstDisplay)
            this.scanvas.home("exact");
        this.firstDisplay = false;
    }
    treetrack.displayTreeSpr = function(treeSpr) {
        var tree = parseNewick(treeSpr[0]["tree"]);
        var spr = treeSpr[1];

        // layout tree
        var ages = getTreeAges(tree);
        var layout = layoutTree(tree, {x:0, y:ages[tree.root], yscale: -1.0});

        // draw tree and spr
        this.scanvas.clear();
        this.scanvas.add(drawSpr(spr, layout));
        this.scanvas.add(drawTree(tree, layout, this.labels));

        // center view
        if (this.firstDisplay) {
            this.scanvas.home("exact");
            this.scanvas.focusCenter();
            this.scanvas.zoom(.9, .9);
        }
        this.firstDisplay = false;

        this.scanvas.draw();
    }
    treetrack.setLabels = function(labels) {
        this.labels = labels;
    }
    mashome.addTrack(treetrack);


    // create site viewer
    var sitetrack = new mashome.Track({name: "local sites", height: 200});
    sitetrack.onViewChange = function (view) {
        this.view = view;
        this.main.css("overflow", "auto");
    };
    sitetrack.showSites = function(chrom, pos) {
        var that = this;
        var tree;

        // fetch local tree.
        $.ajax({
            dataType: 'jsonp',
            url: argtrackurl + '/tree/' + chrom + ':' + pos
        }).then(function(result) {
            var item = JSON.parse(result);
            tree = parseNewick(item.tree);
            var region = chrom + ':' + item.start + '-' + item.end;

            // fetch sites in region of local tree.
            return $.ajax({
                dataType: 'jsonp',
                url: argtrackurl + '/sites/' + region
            });
        }).then(function (result) {
            var sites = JSON.parse(result);
            console.log(tree);
            console.log(sites);

            // draw tree and sites.
            that.plot(tree, sites);
        });
    }
    sitetrack.plot = function (tree, sites) {
        var text = "<pre>";

        var order = getTreeLeaves(tree);
        for (var i in order)
            order[i] = parseInt(order[i]);

        for (var i in sites) {
            var site = sites[i];
            var highFreq = getHighFreq(site.col);

            text += mashome.pad(""+ site.pos, " ", 9) + " ";
            for (var j=0; j<site.col.length; j++) {
                var a = site.col[order[j]];

                if (a == highFreq)
                    text += a;
                else
                    text += "<span style='color: red'>" + a + "</span>";
            }
            text += "\n";
        }

        text += "</pre>";

        this.main.html(text);

    };
    mashome.addTrack(sitetrack);


    // tree track config
    var treeConfigTrack = new mashome.Track({name: "ARG URL",
                                      height: 20});
    treeConfigTrack.onAddTrack = function (view) {
        var that = this;
        this.labelsCookie = "argtrack-tree-labels";
        this.elm.html(
            "<b>tree labels:</b>" +
            "<form id='argtrack-tree-config-form'>" +
            "<textarea id='argtrack-tree-labels'></textarea>" +
            "<input type='submit' value='set labels'>" +
            "</form>");
        this.elm.css({"width": this.width - 21, // offset padding
                      "background-color": "#eee",
                      "border-top": "1px solid #ccc",
                      "padding": "1px",
                      "padding-left": "20px",
                      "text-align": "left"});
        this.elm.find("#argtrack-tree-labels").css({width: 400, height: 20});

        this.labelsInput = this.elm.find("#argtrack-tree-labels");
        this.elm.find("#argtrack-tree-config-form").submit(function(e){
                that.onSubmit(e)});

        this.initLabels();
    };
    treeConfigTrack.onSubmit = function (event) {
        event.preventDefault();
        var text = this.labelsInput.val();
        var cookieTime = 30;
        mashome.setCookie(this.labelsCookie, text, cookieTime);
        this.onLabelsChange(text);
    };
    treeConfigTrack.initLabels = function () {
        var text = mashome.getCookie(this.labelsCookie);
        if (text) {
            this.labelsInput.val(text);
            treetrack.setLabels(text.split("\n"));
        }
    };
    treeConfigTrack.onLabelsChange = function(text){
        treetrack.setLabels(text.split("\n"));
    };
    mashome.addTrack(treeConfigTrack);


    // create url track
    var urlTrack = new mashome.Track({name: "ARG URL",
                                      height: 20});
    urlTrack.onAddTrack = function (view) {
        var that = this;
        this.urlCookie = "argtrack-url";

        this.elm.html(
            "<form id='argtrack-toolbar-form'>" +
            "ARG URL: <input id='argtrack-url' type='text'></input>" +
            "</form>");
        this.elm.css({"width": this.width - 21, // offset padding
                      "background-color": "#eee",
                      "border-top": "1px solid #ccc",
                      "padding": "1px",
                      "padding-left": "20px",
                      "text-align": "left"});
        this.elm.find("#argtrack-url").css("width", 400);

        this.urlInput = this.elm.find("#argtrack-url");
        this.elm.find("#argtrack-toolbar-form").submit(function(e){
                that.onSubmit(e)});

        this.initUrl();
    };
    urlTrack.onSubmit = function (event) {
        event.preventDefault();
        var url = this.urlInput.val();
        var cookieTime = 30;
        mashome.setCookie(this.urlCookie, url, cookieTime);
        this.onUrlChange(url);
    };
    urlTrack.initUrl = function () {
        var url = mashome.getCookie(this.urlCookie);
        if (url) {
            this.urlInput.val(url);
            argtrackurl = url;
        } else {
            var cookieTime = 30;
            this.urlInput.val(argtrackurl);
            mashome.setCookie(this.urlCookie, argtrackurl, cookieTime);
        }
    };
    urlTrack.onUrlChange = function (url) {
        argtrackurl = url;
        mashome.reloadTrackScript();
    };
    mashome.addTrack(urlTrack);


}

// import dependencies
mashome.importScripts(
    [window.arghost + "/static/js/summon.js",
     window.arghost + "/static/js/newick.js"],
    setup);
})();
