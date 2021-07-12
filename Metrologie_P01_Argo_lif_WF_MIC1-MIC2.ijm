/**
 * 20210702_Metrologie_MIC1-MIC2.ijm
 * UTILISABLE POUR MIC1 et MIC2
 * Cette macro ImageJ met en oeuvre le protocole de metrologie P01 pour le microscope
 * a champ large MIC2 de la plateforme de microscopie photonique de l'IGBMC.
 * Images necessaires :
 * - Lame Argolight Argo-M : Mosaique de 4 x 4 carres 
 * - Image de correction de l'illumination 
 * acquises avec LASX, Objectif 40x ou 20x, cube GFP ou L5, Camera Hamamatsu Orca-Flash 4.0 LT
 * Noms des images dans le fichier LIF: 
 * 		¤ Images du motif de 4x4 carres : 
 * 			- le nom doit contenir P01_40x ou P01_20x
 * 			- s'il y a plusieurs images du motif par grossissement, seule la premiere est 
 * 			  traitee.
 * 		¤ Image de correction de l'illumination (flatField) : 
 * 			- le nom doit contenir P02_40x ou P02_20x
 * 			- s'il y a plusieurs images flatField par grossissement, seule la premiere est 
 * 			  utilisee.
 * Par defaut, Le carre le plus intense est suppose être en haut a gauche 
 * de l'image. S'il est a droite, utiliser scanDirection = "RightToLeft"
 * 
 * Marcel Boeglin June-July 2021
 * boeglin@igbmc.fr
 */

/**
 * 20210702_Metrologie_MIC1-MIC2.ijm
 * WORKS WITH MIC1 & MIC2
 * This ImageJ macro implements the Light Microscopy Facility IGBMC Metrology PO1 protocol 
 * for wide-field microscope MIC2.
 * It uses Argolight Argo-M 4 x 4 squares Mosaique and flatField images taken with the LASX
 * software, objective 40x or 20x, fluorescence cube GFP or L5,
 * Camera Hamamatsu Orca-Flash 4.0 LT
 * Image names in the LIF file: 
 * 		¤ 4x4 squares motif : 
 * 			- the name should contain P01_40x or P01_20x
 * 			- if more than one image per magnification, only first one is processed
 * 		¤ flatField correction Image : 
 * 			- the name should contain P02_40x ou P02_20x
 * 			- if more than one image per magnification, only first one is used
 * By default, the brightest square is supposed to be the top-left one;
 * If it's  top-right, use scanDirection = "RightToLeft"
 * 
 * Marcel Boeglin June-July 2021
 * boeglin@igbmc.fr
 */


var inputPath;
var inputDir;
var inputFileName;
var extension;
var argoLightID;
var backgroundID;
var unit, pixelWidth, pixelHeight;
var flatFieldID;
var flatFieldMean = 1;
var nSeries;
var seriesIndex;//index of image in Bio-Format importer
var seriesName;
var outpuDir;
var outputPath;
var outputName;

var microscope = "MIC1";
var objective = "20x";
var argolightImageFilter;
var darkImageFilter;
var flatFieldImageFilter;

var width, height;
var unit;
var pixelWidth, pixelHeight;

var dotsWidth = 50; //dimension nominale des carres en microns;
var dotsWidth = 51.5;//dimension mesuree MIC2
var dotsWidth = 53;//dimension mesuree MIC2 apres process et threshold 

var dotsArrayWdth = 178.65;//MIC2
var dotsArrayWdth = 180.5;//MIC1

/* arrayParameter : mesure comme  montre dans
 * 2021_04_05_Métrologie_MIC2.lif - P01_40x_gfp_900ms_100_Definition_PasDuReseau.png */
var arrayParameter = dotsArrayWdth / 3; // 59.55 microns
var elementarySquareArea = 2652;//micron square, calculated
//var elementarySquareArea = 2700;//micron square (mesure sur image au 40x)
var elementarySquareArea; //moyenne aire des 16 carres segmentes par Triangle

//computation of Rois-grid
var nDots = 16;
var dotCentersX = newArray(nDots);

var dotsGap = arrayParameter - dotsWidth;//microns
var dotsGapPixels;// = dotsGap/pixelWidth;//pixels

var gridCols = 4;
var gridRows = 4;
var gridCenterX, gridCenterY;
var gridWidth;
var gridHeight;
var roisWidth = dotsWidth * 0.85;//physical units, for grid computation

var sampleTilt;

var xStep = dotsArrayWdth / gridCols;
var yStep = dotsArrayWdth / gridRows;
var TopLeft = newArray(2);
var TopRight = newArray(2);
var BottomRight = newArray(2);
var BottomLeft = newArray(2);

var brightestSquarePosition = "TopRight";//"TopLeft" pour MIC2
var scanDirection;


//scanDirection = "RightToLeft" si le carre le plus intense est en haut a droite
//Autre cas (carre le plus intense en bas) non prevus

run("Close All");
print("\\Clear");
run("Bio-Formats Macro Extensions");

getInputPath();
print("inputPath = "+inputPath);
print("\ninputDir :");
print(inputDir);
print("inputFileName = "+inputFileName);
print("extension = "+extension);
inputFile = inputFileName + extension;
print("inputFile = "+inputFile);

outputDir = inputDir + inputFile;
outputDir += "_Mesures";
outputDir += File.separator;
print("outputDir : "+outputDir);
if (!File.exists(outputDir)) {
	print("\nCreating output directory\n ");
	File.makeDirectory(outputDir);
}
else print("\nOutput directory already exists\n ");

selectMicroscope();
getParams();

roisWidth = dotsWidth * 0.85;//physical units, for grid computation

dotsGap = arrayParameter - dotsWidth;//microns

printParams();

//IMAGE D'OBSCURITE
darkImageFilter = "P01D_"+objective;
darkImageName = openSeries(inputDir, inputFile, darkImageFilter);
darkSignal = 0;
if (darkImageName != "") {
	getStatistics(area, darkSignal, min, max, std, histogram);
	close();
}

//IMAGE ARGOLIGHT
argolightImageFilter = "P01_"+objective;
flatFieldImageFilter = "P02_"+objective;
argolightImageName = openSeries(inputDir, inputFile, argolightImageFilter);
if (argolightImageName=="") exit("Couldn't open Argolight Image");
argoLightID = getImageID();
//Ssoustraction fond constant
print("\ndarkSignal = "+darkSignal+"\n ");
if (darkSignal != 0)
	run("Subtract...", "value="+darkSignal);
width = getWidth();
height = getHeight();
getPixelSize(unit, pixelWidth, pixelHeight);
print("inputFileName = "+inputFileName);
print("\nargolightImageName = "+argolightImageName+"\n ");
resetMinAndMax();
outputName = inputFile+" - "+argolightImageName;
print("\noutputName = "+outputName+"n ");
rename("Argolight");
getPixelSize(unit, pixelWidth, pixelHeight);

//IMAGE DE CORRECTION DE L'ILLUMINATION
print("flatFieldImageFilter = "+flatFieldImageFilter);
flatFieldName = openSeries(inputDir, inputFile, flatFieldImageFilter);
//exit;
flatfieldoutputName = inputFile+" - "+flatFieldName;
rename("Flatfield");
run("Gaussian Blur...", "sigma=2");
getStatistics(area, flatFieldMean, min, flatFieldMax, std, histogram);
print("");
print("flatFieldMean = "+flatFieldMean);
flatFieldID = getImageID();
print("");

//Create corrected ArgoLight image for intensity measurements
run("Calculator Plus",
	"i1=Argolight i2=Flatfield operation=[Divide: i2 = (i1/i2) x k1 + k2] k1="+
	flatFieldMean+" k2=0 create");
setVoxelSize(pixelWidth, pixelHeight, 0, unit);
rename("CorrectedArgoLight");
correctedArgoLightID = getImageID();

//Create processed Argolight image for squares detection
run("Duplicate...", " ");
rename("ProcessedArgoLight");
processedArgoLightID = getImageID();

outliersRadiusPixels = 1.29 / pixelWidth;
if (outliersRadiusPixels<4) outliersRadiusPixels = 4;
print("outliersRadiusPixels = "+outliersRadiusPixels);
run("Remove Outliers...", "radius="+outliersRadiusPixels+" threshold=0 which=Bright");
run("Remove Outliers...", "radius="+outliersRadiusPixels+" threshold=0 which=Dark");

dotsGapPixels = dotsGap/pixelWidth;
print("dotsGap = "+dotsGap);
print("pixelWidth = "+pixelWidth);
print("dotsGapPixels = "+dotsGapPixels);
rollingPixels = dotsWidth/pixelWidth/2;
print("rollingPixels = "+rollingPixels);
run("Subtract Background...", "rolling="+(rollingPixels/2)+" sliding disable");

run("Gamma...", "value=0.50");

/*
outliersRadiusPixels = 1.29 / pixelWidth / 2;
print("outliersRadiusPixels = "+outliersRadiusPixels);
run("Remove Outliers...", "radius="+outliersRadiusPixels+" threshold=0 which=Bright");
run("Remove Outliers...", "radius="+outliersRadiusPixels+" threshold=0 which=Dark");
*/
//exit;
getDotsAsRois(processedArgoLightID);

//Create Rois grid and adjust it's position on processedArgoLight image
getDotsArrayCornersFromRoiManager();
sampleTilt = computeSampleTilt(processedArgoLightID);
getGridCenterFromDotsArrayCorners();
createRoisGrid(processedArgoLightID);
aroundImageCenter = true;
rotateRois(sampleTilt, aroundImageCenter);
translation = computeXYShift(processedArgoLightID);
tx = translation[0];
ty = translation[1];
translateRois(tx, ty);//Translate Roi-grid to fit dot positions
suffix = "_measuredImage";
getDimensions(w, h, channels, slices, frames);
print("channels = "+channels);
run("Set Measurements...", "area mean median redirect=None decimal=3");
run("Clear Results");
selectImage(correctedArgoLightID);
for (c=1; c<=channels; c++) {
	if (channels>1) Stack.setChannel(c);
	print("measuring channel"+c);
	Roi.remove;
	roiManager("deselect");
	for (i=0; i<roiManager("count"); i++) {
		updateResults();
		roiManager("select", i);
		roiname = Roi.getName;
		run("Measure");
	}
}
for (i=0; i<nResults; i++) {
	roiManager("select", i);
	setResult("Roi", i, i+1);
	updateResults();
}
saveAs("Results", outputDir+outputName+"_Corrected_Results.txt");

//Add measurement squares to corrected Argolight image
roiManager("Show None");
nRois = roiManager("count");
for (i=0; i<nRois; i++) {
	roiManager("select", i);
	Overlay.addSelection;
	Roi.getBounds(x, y, width, height);
	size = dotsWidth / (pixelWidth*2);
	options = "bold scale";
	Overlay.setLabelFontSize(size, options);
	Overlay.drawLabels(true);
}
Roi.remove;
//run("Fire");
//save corrected Argolight image
resetMinAndMax();
saveAs("Tiff", outputDir+outputName+"_Corrected.tif");


selectImage(flatFieldID);
for (i=0; i<nRois; i++) {
	roiManager("select", i);
	Overlay.addSelection;
	Roi.getBounds(x, y, width, height);
	size = dotsWidth / (pixelWidth*2);
	options = "bold scale";
	Overlay.setLabelFontSize(size, options);
	Overlay.drawLabels(true);
}
	Roi.remove;
	run("Fire");
	saveAs("Tiff", outputDir+flatfieldoutputName+".tif");

if (isOpen(processedArgoLightID)) {
	selectImage(processedArgoLightID);
	//close();
}
selectWindow("Results");
Plot.create("Plot of Results", "Roi", "Median");
Plot.add("Circle",
	Table.getColumn("Roi", "Results"), Table.getColumn("Median", "Results"));
Plot.setStyle(0, "blue,#a0a0ff,4.0,Circle");
Plot.addLegend(outputName+"_Corrected", "Bottom-Left");
Plot.show();
saveAs("PNG", outputDir+outputName+"_Corrected_Results_Median.png");

//Measure uncorrected Argolight image
selectImage(argoLightID);
run("Set Measurements...", "area mean median redirect=None decimal=3");
run("Clear Results");
for (c=1; c<=channels; c++) {
	if (channels>1) Stack.setChannel(c);
	print("measuring channel"+c);
	Roi.remove;
	roiManager("deselect");
	for (i=0; i<roiManager("count"); i++) {
		updateResults();
		roiManager("select", i);
		roiname = Roi.getName;
		run("Measure");
	}
}
for (i=0; i<nResults; i++) {
	roiManager("select", i);
	setResult("Roi", i, i+1);
	updateResults();
}
saveAs("Results", outputDir+outputName+"_Uncorrected_Results.txt");
//Add measurement squares to uncorrected Argolight image
roiManager("Show None");
nRois = roiManager("count");
for (i=0; i<nRois; i++) {
	roiManager("select", i);
	Overlay.addSelection;
	Roi.getBounds(x, y, width, height);
	size = 60;
	size = dotsWidth / (pixelWidth*2);
	options = "bold scale";
	Overlay.setLabelFontSize(size, options);
	Overlay.drawLabels(true);
}
Roi.remove;
//save uncorrected Argolight image
resetMinAndMax();
saveAs("Tiff", outputDir+outputName+"_Uncorrected.tif");
selectWindow("Results");
Plot.create("Plot of Results", "Roi", "Median");
Plot.add("Circle",
	Table.getColumn("Roi", "Results"), Table.getColumn("Median", "Results"));
Plot.setStyle(0, "blue,#a0a0ff,4.0,Circle");
Plot.addLegend(outputName+"_Uncorrected", "Bottom-Left");
Plot.show();
saveAs("PNG", outputDir+outputName+"_Uncorrected_Results_Median.png");
//run("Close All");

print("");
printParams();
//END MACRO

function selectMicroscope() {
	Dialog.create("Metrologie_MIC1 - MIC2");
	microscopes = newArray("MIC1", "MIC2");
	Dialog.addChoice("Microscope", microscopes, microscope);
	Dialog.show();
	microscope = Dialog.getChoice();
	if (microscope=="MIC1") {
		dotsWidth = 53.5;
		dotsArrayWdth = 180.5;
		gridWidth = 180.5;
		gridHeight = 180.5;
		scanDirection = "RightToLeft";
		objective = "20x";
		brightestSquarePosition = "TopRight";
	}
	else if (microscope=="MIC2") {
		dotsWidth = 53;
		dotsArrayWdth = 178.65;
		gridWidth = 178.65;
		gridHeight = 178.65;
		scanDirection = "LeftToRight";
		objective = "40x";
		brightestSquarePosition = "TopLeft";
	}
}

function printParams() {
	print("microscope : "+microscope);
	print("brightestSquarePosition : "+brightestSquarePosition);
	print("scanDirection : "+scanDirection);
	print("objective : "+objective);

	print("arrayParameter = "+arrayParameter);
	print("elementarySquareArea = "+elementarySquareArea);
	print("dotsWidth = "+dotsWidth);
	print("dotsGap = "+dotsGap);
	print("roisWidth = "+roisWidth);


	print("dotsWidth : "+dotsWidth);
	print("gridWidth : "+gridWidth);
	print("gridHeight : "+gridHeight);

	print("inputDir : "+inputDir);
	print("inputFileName : "+inputFileName);
	print("extension : "+extension);
	print("argolightImageFilter : "+argolightImageFilter);
	print("darkImageFilter : "+darkImageFilter);
	print("flatFieldImageFilter : "+flatFieldImageFilter);
	print("pixelWidth : "+pixelWidth);
	print("pixelHeight : "+pixelHeight);
	print("unit : "+unit);
	print("flatFieldMean : "+flatFieldMean);
	print("outputName : "+outputName);
}

function getParams() {
	Dialog.create(microscope);
	brightestSquarePositions = newArray("TopLeft", "TopRight");
	Dialog.addChoice("Brightest Square Position", brightestSquarePositions, brightestSquarePosition);
	objectives = newArray("10x", "20x", "40x");
	Dialog.addChoice("Objective", objectives, objective);
	Dialog.show();
	brightestSquarePosition = Dialog.getChoice();
	if (brightestSquarePosition=="TopLeft") scanDirection = "LeftToRight";
	else scanDirection = "RightToLeft";
	objective = Dialog.getChoice();
}

function getInputPath()  {
	inputPath = File.openDialog("Select input images file");
	inputDir = File.getDirectory(inputPath);
	inputFileName = File.getNameWithoutExtension(inputPath);
	extension = substring(inputPath, lastIndexOf(inputPath, "."));
}

function openSeries(dir, file, seriesFilter) {
	/* Adapted from macro "extract_series_from_leica_lif.ijm"
	* July 2019
	* Erwan Grandgirard  & Bertand Vernay
	* grandgie@igbmc.fr & vernayb@igbmc.fr
	*/
	//print("openSeries(dir, file, seriesFilter),\nfile = "+file+"\nseriesFilter = "+seriesFilter);
	msg = "openSeries(dir, file, seriesFilter) failed";
	path=dir+file;
	Ext.setId(path);
	Ext.getCurrentFile(file);
	Ext.getSeriesCount(seriesCount);//gets the number of series in input file
	//print("Processing the file = " + file);
// See:
//http://imagej.1557.x6.nabble.com/multiple-series-with-bioformats-importer-td5003491.html
	//while next size is a fraction of current (1/3 for CZI, 1/2 for BDP2 HDF5
	//the images belong to the same pyramidal series with different resolutions
	Ext.getCurrentFile(currentFile);
	print("file = "+currentFile);
	seriesname = "";
	for (j=0; j<seriesCount; j++) {
		//print("Extracting Series "+j+1+" / "+seriesCount);
		Ext.setSeries(j);
		Ext.getSeriesName(seriesName);
		//print("series "+j+" :  seriesName = "+seriesName);
		Ext.getUsedFileCount(count);
//		print("usedFileCount = "+count);
		Ext.getSizeX(sizeX);
        Ext.getSizeY(sizeY);
//		print("sizeX = "+sizeX+"    sizeY = "+sizeY);
		//print("raw seriesName = "+seriesName);
		str = replace(seriesName, "\"", "");
		STR = toUpperCase(str);
		SERIESNAME = toUpperCase(seriesName);
		SERIESFILTER = toUpperCase(seriesFilter);
		if (indexOf(SERIESNAME, SERIESFILTER) >= 0) {
			run("Bio-Formats Importer",
				"open=&path color_mode=Default view=Hyperstack stack_order=XYCZT series_"+j+1);
			seriesname = seriesName;
			break;
		}
	}
	print("Opened seriesName : "+seriesName);
	//print("openSeries(dir, file, seriesFilter) END");
	return seriesname;
}

/** Processes, segments and analyzes active channel of image 'id'
	and stores the dots as Rois in RoiManager.
	This is the critical part of the procedure. May fail if dots are irregular
	or have unequal intensities. */
function getDotsAsRois(id) {
	print("getDotsAsRois(id):");
	selectImage(id);
	Roi.remove;
	run("Remove Overlay");
	setOption("BlackBackground", true);
	resetMinAndMax;
	setAutoThreshold("Triangle dark");
	roiManager("reset");
	run("Set Measurements...", "centroid display redirect=None decimal=3");
	tolerance = 20;// %
	tolerance = 25;// %
	//tolerance = 15;// < 15 => problemes de detection

	if (isOpen("ROI Manager")) roiManager("reset");
	minArea = elementarySquareArea*(100-tolerance)/100;
	maxArea = elementarySquareArea*(100+tolerance)/100;
//	minArea = elementarySquareArea;
//	minArea -= minArea*tolerance/100;
//	maxArea = elementarySquareArea;
//	maxArea += maxArea*tolerance/100;
	print("minArea = "+minArea);
	print("maxArea = "+maxArea);
	print("objective = "+objective);

	run("Threshold...");
	//setAutoThreshold("Huang dark");// MIC1 et MIC2
	//setAutoThreshold("Li dark");// MIC1 et MIC2
	//setAutoThreshold("Yen dark");// MIC1 et MIC2
	//setAutoThreshold("Mean dark");// MIC1
	setAutoThreshold("RenyiEntropy dark");// MIC1 et MIC2
	//setAutoThreshold("Triangle dark");// MIC2
	setTool("rectangle");
	//waitForUser("Draw a rectangle around the array of fluorescent squares");
	waitForUser("You can draw a Roi, \nchange the threshold or\nedit image with pencil to separate or close squares");

	//exit;
//	wait(1000);
	print("minDotArea = "+minArea+"   maxDotArea = "+maxArea);
	run("Analyze Particles...",
		"size="+minArea+"-"+maxArea+" exclude display clear include add");
	detectedDots = roiManager("count");
	print("Detected "+detectedDots+" dots");
	resetThreshold;
	n=gridRows*gridCols;
	if (detectedDots != gridRows*gridCols) {
		print("An error occured in dots detection: should find "+n+"dots");
		close();
		print("End getDotsAsRois(id):");
		return 0;
	}
	print("End getDotsAsRois(id):");
	return getImageID();
}

function getGridCenterFromDotsArrayCorners() {
	print("getGridCenterFromDotsArrayCorners()");
	if (roiManager("count")!=nDots) {
		print("dotsDetectionFailed");
		print("End getGridCenterFromDotsArrayCorners()");
		return;
	}
	//corners coordinates are in physical units
	gridCenterX = (TopLeft[0]+TopRight[0]+BottomRight[0]+BottomLeft[0])/4;
	gridCenterY = (TopLeft[1]+TopRight[1]+BottomRight[1]+BottomLeft[1])/4;
	print("gridCenterX = "+gridCenterX+" "+unit);
	print("gridCenterY= "+gridCenterY+" "+unit);
	print("End getGridCenterFromDotsArrayCorners()");
}

/* Creates a grid of Rois */
function createRoisGrid(imageID) {
//function createCircularRoisGridForResolution3(imageID) {
	print("createCircularRoisGridForResolution3(imageID)");
	print("gridCols = "+gridCols);
	print("gridRows = "+gridRows);
	xStep = gridWidth/(gridCols-1);
	yStep = gridHeight/(gridRows-1);
	print("xStep = "+xStep);
	print("yStep = "+yStep);
	print("scanDirection = "+scanDirection);
	print("roisWidth = "+roisWidth+" "+unit);
	getGridCenterFromDotsArrayCorners();
	print("gridWidth = "+gridWidth);
	print("gridHeight = "+gridHeight);
	physicalUnits = true;

	//anciennement createRoiGridForResolution3
	createRoiGrid(gridCenterX, gridCenterY,
		gridCols, gridRows, xStep, yStep,
		scanDirection, roisWidth, physicalUnits);
	print("End createRoiGrid(imageID)");
}

/*	Creates an grid of rois and stores them in the roiManager;
	centerX, centerY: coordinates of grid-center
	xStep, yStep: x and y periods
	scanSirection: "RightToLeft" or "LeftToRight"
	roisWidth: size  of the rois
	If !physicalUnits, lengths and positions must be passed in pixels */
function createRoiGrid(centerX, centerY, cols, rows,
		xStep, yStep, scanDirection, roisWidth, physicalUnits) {
	roiManager("reset");
	//cx, cy, xs, ys: center coordinates and steps in pixels
	cx = centerX; cy = centerY;
	xs = xStep; ys = yStep;
	roisWidthPixels = roisWidth / pixelWidth;

	//convert lengths and positions to pixels
	cx /= pixelWidth; cy /= pixelHeight;
	xs /= pixelWidth; ys /= pixelHeight;
	print("roisWidthPixels = "+roisWidthPixels);
	
	str = ""+cols*rows;
	digits = str.length;
	i=0;
	//roiCenterX, roiCenterY in pixels
	for (r=0; r<rows; r++) {
		roiCenterY = cy - ys*(rows-1)/2 + r*ys;
		if (scanDirection=="RightToLeft") {
			for (c=cols-1; c>=0; c--) {
				roiCenterX = cx - xs*(cols-1)/2 + c*xs;
				makeRectangle(roiCenterX-roisWidthPixels,
					roiCenterY-roisWidthPixels,
					1*roisWidthPixels, 1*roisWidthPixels);
				setSelectionName(String.pad(++i, digits));
				roiManager("add");
			}
		}
		else if (scanDirection=="LeftToRight") {
			for (c=0; c<cols; c++) {
				roiCenterX = cx - xs*(cols-1)/2 + c*xs;
				makeRectangle(roiCenterX-roisWidthPixels,
					roiCenterY-roisWidthPixels,
					1*roisWidthPixels, 1*roisWidthPixels);
				setSelectionName(String.pad(++i, digits));
				roiManager("add");
			}
		}
	}
	roiManager("deselect");
	Roi.remove;
	print("End createRoiGrid()");
}

/* Computes corners coordinates in physical units (microns) */
function getDotsArrayCornersFromRoiManager() {
	print("getDotsArrayCornersFromRoiManager()");
	xmin = width*pixelWidth; ymin = height*pixelHeight;//TopLeft dot
	xmax = 0; ymax = 0;//BottomRight dot
	run("From ROI Manager");
	nRois = getValue("results.count");
	print("nRois = "+nRois);
	print("nDots = "+nDots);
	if (nRois!=nDots) {
		exit(""+nRois+" Rois detectees au lieu de 16");
	}
	for (i=0; i<nDots; i++) {
		x = getResult("X", i);
		y = getResult("Y", i);
		if (x+y < xmin+ymin) {
			xmin=x; ymin=y;
		}
		if (x+y > xmax+ymax) {
			xmax=x; ymax=y;
		}
	}
	TopLeft[0]=xmin; TopLeft[1]=ymin;
	BottomRight[0]=xmax; BottomRight[1]=ymax;
	print("TopLeft[0]="+xmin+"  TopLeft[1]="+ymin+
		"  BottomRight[0]="+xmax+"  BottomRight[1]="+ymax);
	xmax=0; ymin=height*pixelHeight;//TopRight dot
	xmin=width*pixelWidth; ymax=0;//BottomLeft dot
	nDots = getValue("results.count");
	for (i=0; i<nDots; i++) {
		x = getResult("X", i);
		y = getResult("Y", i);
		if (x-y < xmin-ymax) {
			xmin=x; ymax=y;
		}
		if (x-y > xmax-ymin) {
			xmax=x; ymin=y;
		}
	}
	TopRight[0]=xmax; TopRight[1]=ymin;
	BottomLeft[0]=xmin; BottomLeft[1]=ymax;
	print("BottomLeft[0]="+xmin+"  TopRight[1]="+ymin+
			"  TopRight[0]="+xmax+"  BottomLeft[1]="+ymax);
	print("End getDotsArrayCornersFromRoiManager()");
}

/** Computes the angle by which rotate the Roi-grid to fit the dots array
	in image to be analyzed 
	For that, detects centers of TopLeft, BottomRight, TopRight and BottomLeft dots*/
function computeSampleTilt(imgID) {
	print("computeSampleTilt(imgID)");
	diagonalAngle = 45;
	if (gridRows!=gridCols)
		diagonalAngle = Math.atan(gridRows/gridCols);//to be verified
	print("BEFORE getDotsArrayCornersFromRoiManager()");
	getDotsArrayCornersFromRoiManager();
	print("After getDotsArrayCornersFromRoiManager()");
	fitSquare = true;
	run("Remove Overlay");
	print("TopLeftX="+TopLeft[0]+"  TopLeftY="+TopLeft[1]);
	print("BottomRightX="+BottomRight[0]+"  BottomRightY="+BottomRight[1]);
	//makeSelection("angle",newArray(x1,x2,x3),newArray(y1,y2,y3));
	makeSelection("angle",newArray(BottomRight[0]/pixelWidth,TopLeft[0]/pixelWidth,width),
		newArray(BottomRight[1]/pixelWidth,TopLeft[1]/pixelWidth,TopLeft[1]/pixelWidth));
	run("Add Selection...");

	run("Clear Results");
	run("Measure");
	tilt1 = getResult("Angle", 0) - diagonalAngle;

	print("tilt1="+tilt1);
	//2nd measurement of tilt, to be averaged with 1st:
	print("TopRightX="+TopRight[0]+"  TopRightY="+TopRight[1]);
	print("BottomLeftX="+BottomLeft[0]+" BottomLeftY="+BottomLeft[1]);
	x1 = width*pixelWidth; x2 = BottomLeft[0]; x3 = TopRight[0];
	makeSelection("angle",newArray(width,BottomLeft[0]/pixelWidth,TopRight[0]/pixelWidth),
		newArray(BottomLeft[1]/pixelWidth,BottomLeft[1]/pixelWidth,TopRight[1]/pixelWidth));
	run("Add Selection...");
	run("Remove Overlay");
	run("Clear Results");
	run("Measure");
	tilt2 = diagonalAngle - getResult("Angle", 0);
	print("tilt2="+tilt2);
	print("End computeSampleTilt(imgID)");
	tilt = (tilt1+tilt2)/2;
	print("tilt = "+tilt);
	print("End computeSampleTilt(imgID)");
	return tilt;
}

/**
 * Rotates all rois in roiManager by 'angle' degrees:
 * around image center if aroundImageCenter is true;
 * aroud roi center otherwise.
	//test-code:
		aroundImageCenter = true;
		angle = 12; //degrees
		rotateRois(12, aroundImageCenter);
 */
function rotateRois(angle, aroundImageCenter) {
	nrois = roiManager("count");
	param = "";
	if (aroundImageCenter) param = "rotate ";
	for (i=0; i<nrois; i++) {
		roiManager("select", i);
		run("Rotate...", param+" angle="+angle);
		roiManager("update");
	}
	roiManager("deselect");
	Roi.remove;
}

/** Translates all rois in roiManager by 'tx', 'ty'
	//test-code:
		tx=10; ty=20;
		translateRois(tx, ty);
*/
function translateRois(tx, ty) {
	nrois = roiManager("count");
	for (i=0; i<nrois; i++) {
		roiManager("select", i);
		getSelectionBounds(x, y, w, h);
		Roi.move(x+tx, y+ty);
		roiManager("update");
	}
	roiManager("deselect");
	Roi.remove;
}

/** Returns the translation to be applied to the Roi-grid after it has been
	rotated around the image center to fit the dots in brightfield image.
	Uses the four corners of the grid for better precision
	Components of translation are expressed in pixels */
function computeXYShift(imageID) {
	print("computeXYShift()");
	selectImage(imageID);
	getPixelSize(unit, pixWidth, pixHeight);
	if (nDots!=gridRows*gridCols) {
		print("An error occured in dots detection");
		print("tx="+0+"  ty="+0);
		print("End computeXYShift()");
		return newArray(0,0);
	}
	roiManager("select", 0);//TopLeft
	getSelectionBounds(x0, y0, w0, h0);
	Roi.remove;
	//Corner coordinates are returned in physical units
	//while getSelectionBounds returns values in pixels
	tx0 = TopLeft[0]/pixWidth - (x0+w0/2);
	ty0 = TopLeft[1]/pixHeight - (y0+h0/2);
	print("tx0="+tx0+"\nty0="+ty0);

	rows = sqrt(nDots);
	roiManager("select", rows-1);//TopRight
	getSelectionBounds(x1, y1, w1, h1);
	Roi.remove;
	tx1 = TopRight[0]/pixWidth - (x1+w1/2);
	ty1 = TopRight[1]/pixHeight - (y1+h1/2);
	print("tx1="+tx1+"\nty1="+ty1);

	roiManager("select", nDots-1);//BottomRight
	getSelectionBounds(x2, y2, w2, h2);
	Roi.remove;
	tx2 = BottomRight[0]/pixWidth - (x2+w2/2);
	ty2 = BottomRight[1]/pixHeight - (y2+h2/2);
	print("tx2="+tx2+"\nty2="+ty2);

	roiManager("select", nDots-rows);//BottomLeft
	getSelectionBounds(x3, y3, w3, h3);
	Roi.remove;
	tx3 = BottomLeft[0]/pixWidth - (x3+w3/2);
	ty3 = BottomLeft[1]/pixHeight - (y3+h3/2);
	print("tx3="+tx3+"t\ny3="+ty3);

	//average to increase precision
	tx = (tx0+tx1+tx2+tx3)/4;
	ty = (ty0+ty1+ty2+ty3)/4;
	print("tx="+tx+"\nty="+ty);

	print("End computeXYShift()");
	return newArray(tx,ty));
}
