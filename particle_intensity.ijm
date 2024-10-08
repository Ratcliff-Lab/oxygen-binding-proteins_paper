// This macro processes all the images in a folder and any subfolders.

  extension = ".tif"; ///.tif
  dir1 = getDirectory("Choose Source Directory ");
  dir2 = getDirectory("Choose Destination Directory ");
  setBatchMode(true); // batchmode, images displayed
  n = 0;
  processFolder(dir1);

  function processFolder(dir1) { 
     list = getFileList(dir1); // list of files in directory
     for (i=0; i<list.length; i++) { // go through list of files in directory by increments of 1, list.length
          if (endsWith(list[i], "/")) // 
              processFolder(dir1+list[i]); // if index is a folder, go through function again to access within it
          else if (endsWith(list[i], extension)) 
             processImage(dir1, list[i]); // if index is a .tif file, analyze file with function below
      }
  }
  
  function processImage(dir1, name) {
     open(dir1+name); //open file
     print(n++, name); //print a table of index and file name
     /// reorderding slices in stack
     	// open image as a stack in default mode, image>stacks>tools>stack sorter>reverse, then save 
     ////run("Stack to Images");
     // BF image
     fname = name + " - C=0"; //
     newname = replace(dir2+name, ".tif", "");
     selectWindow(fname); // - C=0
     run("Enhance Contrast...", "saturated=0.5");
     run("Despeckle");
	 setOption("ScaleConversions", true); //
     run("8-bit");
     setAutoThreshold("MaxEntropy dark");
     setThreshold(0 , 110);
     setOption("BlackBackground", true);
     run("Threshold..."); 
     run("Convert to Mask"); 
     run("Erode");
     run("Open");
     run("Erode");
     //run("Erode");
     run("Fill Holes");
     //run("Open");
     run("Dilate");
     run("Analyze Particles...", "size=20-3500 circularity=0.18-1 show=Overlay exclude include add");
     roiManager("Save", newname + ".zip");
     
     // to expand ROIs https://forum.image.sc/t/extract-each-roi-dilate-each-roi-and-measure-area-individually/26841/2
    n = roiManager("count");
	for (i = 0; i < n; i++) {
    roiManager("select", i);
    run("Enlarge...", "enlarge=2");
    roiManager("update")//Replaces the selected ROI on the list with the current selection.
    roiManager("select", i);
    //roiManager("Measure");
	}
	//run("Select All");
	roiManager("Show All with labels");
	roiManager("Save", newname + "_dilate.zip");
	roiManager("Delete");
	//run("Close");
    
    selectWindow(name + " - C=1");
    run("Enhance Contrast...", "saturated=0.75");
    setOption("ScaleConversions", true); //
    run("8-bit");
    roiManager("Delete");
    roiManager("Open", newname + "_dilate.zip");
    //roiManager("Show All");
    roiManager("Measure");
    selectWindow("Results"); 
    saveAs("Results", newname + ".csv");
    selectWindow("Results"); //
    run("Close"); //
    run("Clear Results");
    roiManager("Delete");
    //run("Close");
	close(); 
    run("Collect Garbage");
  }
  
  
  
  // comm -13 <(sort -u file_1) <(sort -u file_2) {run this on command line to check if output files are the same}
  