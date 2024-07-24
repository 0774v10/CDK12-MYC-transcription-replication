//------this macro measures the intensity on TxRed channel in the nuclear area-------
// This script was written by Adrian Andronache, PhD

datadir = getDirectory("Choose the data directory");
filelist = getFileList(datadir);
savedir = getDirectory("Choose the save directory");
print(savedir);

run("Bio-Formats Macro Extensions");

run("Set Measurements...", "area mean integrated display redirect=None decimal=0");

for (i = 0; i < filelist.length; i++) {

	path = datadir + filelist[i];

	if (endsWith(filelist[i], "lif")) {
	
		print(path);
		Ext.setId(path);
		Ext.getCurrentFile(filename);
		Ext.getSeriesCount(seriesCount);
		
		for (j = 1; j <= seriesCount; j++) {
		//for (j = 1; j <= 49; j++) {
	
		dapi = "C2";
		ch_an = "C1";
		roiManager("reset");
		run("Bio-Formats Importer", "open=&path color_mode=Colorized rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+j);
		filename = File.nameWithoutExtension;
		print(filename);
		title = getTitle();
	
		getDimensions(width, height, channels, slices, frames);
	
		if (channels == 2) {
			run("Duplicate...", "duplicate");
			rename(title + "_dup");
			run("Split Channels");
			selectWindow(dapi + "-" + title + "_dup");
			run("Median...", "radius=5 stack");
			setAutoThreshold("Otsu dark");
			setOption("BlackBackground", true);
			run("Convert to Mask", "method=Otsu background=Dark calculate black");
			run("Fill Holes");
			run("Watershed");
			//run("Close-");
			run("Analyze Particles...", "size=150-600 circularity=0.60-1.00 exclude add stack");
		
			selectWindow(ch_an + "-" + title + "_dup");
			run("Subtract Background...", "rolling=100");
			run("Gaussian Blur...", "sigma=1");
			roiManager("deselect");
			roiManager("measure");
			//waitForUser("");
			// Do the processing here by adding your own code.
			// Leave the print statements until things work, then remove them.
			//print("Processing: " + input + File.separator + file);
			//print("Saving to: " + output);
			close("*");
			
			}
	
		
		}

	
		savepath = savedir + filename;
		print(savepath);
		saveAs("Results", savepath + ".txt");
		selectWindow("Results");
		run("Close");
	}

}