/**
* Newick format parser in JavaScript.
*
*
*
* Part of this code is derived from Jason Davies's newick.js which has the 
* following copyright notice.
*
* Copyright (c) Jason Davies 2010.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*/


// Returns a tree in a flat JSON format:
//
// tree = {nodes: {A: {name: "A",
//                     dist: .2,
//                     children: ["B", "C"],
//                     parent: null},
//                 B: {name: "B",
//                     dist: .4,
//                     children: [],
//                     parent: "A"},
//                 C: {name: "C",
//                     dist: .5,
//                     children: [],
//                     parent: "A"}},
//          root: "A"};
//
function parseNewick(text) {
    var ancestors = [];
    var node = {parent: null,
                children: [],
                dist: 0.0};
    var tree = {nodes: {},
                root: node};
    var nodes = [node];
    var tokens = text.split(/\s*(;|\(|\)|,|:|\[|\])\s*/);
    var nextName = 1;

    for (var i=0; i<tokens.length; i++) {
        var token = tokens[i];
        
        switch (token) {
        case '(': // new branchset
            var child = {parent: node, 
                         children: [],
                         dist: 0.0};
            nodes.push(child);
            node.children.push(child);
            ancestors.push(node);
            node = child;
            break;
        case ',': // another branch
            var parent = ancestors[ancestors.length-1];
            var child = {parent: parent,
                         children: [],
                         dist: 0.0};
            nodes.push(child);
            parent.children.push(child);
            node = child;
            break;
        case ')': // optional name next
            node = ancestors.pop();
            break;
        case ':': // optional length next
            break;
        case '[': // optional comment next
            // advance to closing ]
            var j = i+1;
            while (tokens[i] != ']') i++;
            // TODO: parse comment
            //node.comment = "".join(tokens[j:i-1])
            i++;
            break;
        default:
            var x = tokens[i-1];
            if (x == ')' || x == '(' || x == ',') {
                node.name = token;
            } else if (x == ':') {
                node.dist = parseFloat(token);
            }
        }
    }


    // setup node names
    for (var i in nodes) {
        var node = nodes[i];

        // ensure name is set
        if (typeof node.name == "undefined" || node.name === null ||
            node.name == "")
            node.name = nextName++;

        tree.nodes[node.name] = node;
    }
    
    // rewite node names
    for (var name in tree.nodes) {
        var node = tree.nodes[name];
        if (node.parent)
            node.parent = node.parent.name;
        for (var i in node.children) {
            node.children[i] = node.children[i].name;
        }
    }
    tree.root = tree.root.name;
    
    return tree;
}



/*
* Example tree (from http://en.wikipedia.org/wiki/Newick_format):
*
* Newick format:
* (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;
*
* Converted to JSON:
* {
* name: "F",
* branchset: [
* {name: "A", length: 0.1},
* {name: "B", length: 0.2},
* {
* name: "E",
* length: 0.5,
* branchset: [
* {name: "C", length: 0.3},
* {name: "D", length: 0.4}
* ]
* }
* ]
* }
*
* Converted to JSON, but with no names or lengths:
* {
* branchset: [
* {}, {}, {
* branchset: [{}, {}]
* }
* ]
* }
*/
function parseNewickNested(text) {
    var ancestors = [];
    var tree = {};
    var tokens = text.split(/\s*(;|\(|\)|,|:)\s*/);
    
    for (var i=0; i<tokens.length; i++) {
        var token = tokens[i];
        
        switch (token) {
        case '(': // new branchset
            var subtree = {};
            tree.branchset = [subtree];
            ancestors.push(tree);
            tree = subtree;
            break;
        case ',': // another branch
            var subtree = {};
            ancestors[ancestors.length-1].branchset.push(subtree);
            tree = subtree;
            break;
        case ')': // optional name next
            tree = ancestors.pop();
            break;
        case ':': // optional length next
            break;
        default:
            var x = tokens[i-1];
            if (x == ')' || x == '(' || x == ',') {
                tree.name = token;
            } else if (x == ':') {
                tree.length = parseFloat(token);
            }
        }
    }
    return tree;
}