// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, true, false, true, true, false, false, false, true, false, false, false, false, false, true, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, true, true, false, false, false, false, false, true, true, true, false, false, false, false, false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false, false ];
var arrayMetadata    = [ [ "1", "GSM472023", "1" ], [ "2", "GSM472024", "2" ], [ "3", "GSM472025", "3" ], [ "4", "GSM472026", "4" ], [ "5", "GSM472027", "5" ], [ "6", "GSM472028", "6" ], [ "7", "GSM472029", "7" ], [ "8", "GSM472030", "8" ], [ "9", "GSM472031", "9" ], [ "10", "GSM472032", "10" ], [ "11", "GSM472033", "11" ], [ "12", "GSM472034", "12" ], [ "13", "GSM472035", "13" ], [ "14", "GSM472036", "14" ], [ "15", "GSM472037", "15" ], [ "16", "GSM472038", "16" ], [ "17", "GSM472039", "17" ], [ "18", "GSM472040", "18" ], [ "19", "GSM472041", "19" ], [ "20", "GSM472042", "20" ], [ "21", "GSM472043", "21" ], [ "22", "GSM472044", "22" ], [ "23", "GSM472045", "23" ], [ "24", "GSM472046", "24" ], [ "25", "GSM472047", "25" ], [ "26", "GSM472048", "26" ], [ "27", "GSM472049", "27" ], [ "28", "GSM472050", "28" ], [ "29", "GSM472051", "29" ], [ "30", "GSM472052", "30" ], [ "31", "GSM472053", "31" ], [ "32", "GSM472054", "32" ], [ "33", "GSM472055", "33" ], [ "34", "GSM472056", "34" ], [ "35", "GSM472057", "35" ], [ "36", "GSM472058", "36" ], [ "37", "GSM472059", "37" ], [ "38", "GSM472060", "38" ], [ "39", "GSM472061", "39" ], [ "40", "GSM472062", "40" ], [ "41", "GSM472063", "41" ], [ "42", "GSM472064", "42" ], [ "43", "GSM472065", "43" ], [ "44", "GSM472066", "44" ], [ "45", "GSM472067", "45" ], [ "46", "GSM472068", "46" ], [ "47", "GSM472069", "47" ], [ "48", "GSM472070", "48" ], [ "49", "GSM472071", "49" ], [ "50", "GSM472072", "50" ], [ "51", "GSM472073", "51" ], [ "52", "GSM472074", "52" ], [ "53", "GSM472075", "53" ], [ "54", "GSM472076", "54" ], [ "55", "GSM472077", "55" ], [ "56", "GSM472078", "56" ], [ "57", "GSM472079", "57" ], [ "58", "GSM472080", "58" ], [ "59", "GSM472081", "59" ], [ "60", "GSM472082", "60" ], [ "61", "GSM472083", "61" ], [ "62", "GSM472084", "62" ], [ "63", "GSM472085", "63" ], [ "64", "GSM472086", "64" ], [ "65", "GSM472087", "65" ], [ "66", "GSM472088", "66" ], [ "67", "GSM472089", "67" ], [ "68", "GSM472090", "68" ], [ "69", "GSM472091", "69" ], [ "70", "GSM472092", "70" ], [ "71", "GSM472093", "71" ], [ "72", "GSM472094", "72" ], [ "73", "GSM472095", "73" ], [ "74", "GSM472096", "74" ], [ "75", "GSM472097", "75" ], [ "76", "GSM472098", "76" ], [ "77", "GSM472099", "77" ], [ "78", "GSM472100", "78" ], [ "79", "GSM472101", "79" ], [ "80", "GSM472102", "80" ], [ "81", "GSM472103", "81" ], [ "82", "GSM472104", "82" ], [ "83", "GSM472105", "83" ], [ "84", "GSM472106", "84" ], [ "85", "GSM472107", "85" ], [ "86", "GSM472108", "86" ], [ "87", "GSM472109", "87" ], [ "88", "GSM472110", "88" ], [ "89", "GSM472111", "89" ], [ "90", "GSM472112", "90" ], [ "91", "GSM472113", "91" ], [ "92", "GSM472114", "92" ], [ "93", "GSM472115", "93" ], [ "94", "GSM472116", "94" ], [ "95", "GSM472117", "95" ], [ "96", "GSM472118", "96" ], [ "97", "GSM472119", "97" ], [ "98", "GSM472120", "98" ], [ "99", "GSM472121", "99" ], [ "100", "GSM472122", "100" ], [ "101", "GSM472123", "101" ], [ "102", "GSM472124", "102" ], [ "103", "GSM472125", "103" ], [ "104", "GSM472126", "104" ], [ "105", "GSM472127", "105" ], [ "106", "GSM472128", "106" ], [ "107", "GSM472129", "107" ], [ "108", "GSM472130", "108" ], [ "109", "GSM472131", "109" ], [ "110", "GSM472132", "110" ], [ "111", "GSM472133", "111" ], [ "112", "GSM472134", "112" ], [ "113", "GSM472135", "113" ], [ "114", "GSM472136", "114" ], [ "115", "GSM472137", "115" ], [ "116", "GSM472138", "116" ], [ "117", "GSM472139", "117" ], [ "118", "GSM472140", "118" ], [ "119", "GSM472141", "119" ], [ "120", "GSM472142", "120" ], [ "121", "GSM472143", "121" ], [ "122", "GSM472144", "122" ], [ "123", "GSM472145", "123" ], [ "124", "GSM472146", "124" ], [ "125", "GSM472147", "125" ], [ "126", "GSM472148", "126" ], [ "127", "GSM472149", "127" ], [ "128", "GSM472150", "128" ], [ "129", "GSM472151", "129" ], [ "130", "GSM472152", "130" ], [ "131", "GSM472153", "131" ], [ "132", "GSM472154", "132" ], [ "133", "GSM472155", "133" ], [ "134", "GSM472156", "134" ], [ "135", "GSM472157", "135" ], [ "136", "GSM472158", "136" ], [ "137", "GSM472159", "137" ], [ "138", "GSM472160", "138" ], [ "139", "GSM472161", "139" ], [ "140", "GSM472162", "140" ], [ "141", "GSM472163", "141" ], [ "142", "GSM472164", "142" ], [ "143", "GSM472165", "143" ], [ "144", "GSM472166", "144" ], [ "145", "GSM472167", "145" ], [ "146", "GSM472168", "146" ], [ "147", "GSM472169", "147" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
