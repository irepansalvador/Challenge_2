'use strict';


var infile = process.argv[2];
var out = infile.substring(0, infile.length - 3); // "12345.0"
out = out + ".json";
// Function to read a txt file (the newick tree)
var XMLHttpRequest = require("xmlhttprequest").XMLHttpRequest;
var xhr = new XMLHttpRequest();
function readTextFile(file)
    {
    var rawFile = new XMLHttpRequest();
    rawFile.open("GET", file, false);
    rawFile.onreadystatechange = function ()
        {
        if(rawFile.readyState === 4)
            {
            if(rawFile.status === 200 || rawFile.status == 0)
                {newick = rawFile.responseText;
                console.log(newick);}
        }}
    rawFile.send(null);
    }
/*--- Function to parse a newick file and converto to json --*/
const myModule = require('./newick_mod.js');
// Functions required to export the json file
const fs = require('fs');
//var newick = "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;";
var newick; // create an empty variable to store the txt file
readTextFile("file:///home/irepan/Desktop/Github/Challenge_2/training_set/" + infile)

var json = myModule.parse(newick);
console.log(json);

let data = JSON.stringify(json);
fs.writeFileSync(out, data);
