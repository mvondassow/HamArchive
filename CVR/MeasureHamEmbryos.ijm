/* I wrote this macro to measure the size of Haminoea vesciula embryos.
It identifies dark blobs (embryos), separately from the surrounding dark ring (the 
envelope around the embryo). The images this was designed for were all grayscale
at the same magnification, with H. v. embryos, which have an envelope lifted well off
the embryo. For other embryos or magnifications, the parameters  would need to be
adjusted.

The basic steps are: 
1) identify centers of large homogenous dark blobs (embryo centers) by thresholding 
after a strong Gaussian blur in a duplicate image.
2) enhance separation of the interior embryo (of interest) from the exterior (envelope) 
to be ignored by using unsharp mask to enhance edges and valleys in brightness,
and averaging that sharpened image with the very blurred image to weight points
nearer to the center of the embryo higher than those further away.
3) then threshold and find above-threshold pixels that in regions that connect to the
embryo center.
4) finally, modify that selection by expanding/contracting the selection to remove close
some of the inlets and penninsulas.
5) measure. 
*/
macro "measure embryos 5"
{	run("Set Scale...", "distance=0");
	run("Set Measurements...", "stack display redirect=None decimal=3");
	// Duplicate the image stack and blur to make a stack for finding the center of the large
	// blobs.
	im1 = getImageID();
	run("Duplicate...", "title=A duplicate range=1-" + nSlices());
	im2 = getImageID()
	run("Gaussian Blur...", "sigma=50 stack");
	// Go back to the original image stack, use Unsharp Mask to enhance edges and valleys.
	// Then use the blurred image to weight pixels near the blob (embryo) center more highly
	// than pixels further away (e.g. the envelope).
	selectImage(im1);
	run("Unsharp Mask...", "radius=25 mask=0.50 stack");
	imageCalculator("Average stack", im1,im2);
	selectImage(im2);
	// Threshold and analyze particles to find centers.
	setAutoThreshold("Default stack");
	setOption("BlackBackground", false);
	// Modified for cleavers embryo 12 and on to calculate threshold for each frame.
	run("Convert to Mask", "method=Default background=Light calculate");
	roiManager("reset");
	// Parameters (size and circularity) chosen to match the particular images used 
	// while designing this.
	run("Analyze Particles...", "size=3000-200000 circularity=0.5-1.00 display exclude clear include add in_situ stack");
	run("Clear Results");
	// Go back to image 1 and threshold it with a threshold it.
	selectImage(im1);
	setAutoThreshold("Default");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Default background=Light calculate");
	run("Options...", "iterations=1 count=1 pad do=Nothing");
	// Use "close" to try to eliminate some missing pieces at the edge of the blob where
	// there are below threshold spots on the perimeter leading to clear spots (e.g. nuclei) 
	// near the surface (e.g. in cleavage stage embryos).
	run("Close-", "stack");
	// feret's diameter and min feret's diameter seem to be less sensitive to the error due
	// to missing bits at the edge of the embryo, so feret seems better for calculating volumes. 
	run("Set Measurements...", "area fit feret's stack display redirect=None decimal=3");
	// The following loop finds the area of the desired blob (i.e. an embryo). It goes to the
	// centers of blobs found using "analyze particles" on the blurred image (based on the 
	// ROIs), and finds the region of above threshold pixels connected to those centers.
	maxk = roiManager("count");
	for(k=0;k<maxk;k++) {
			roiManager("select", k);
			Roi.getBounds(x,y,width,height);
			doWand(x+width/2, y+height/2, 0.0, "4-connected");
			// the enlarging and shrinking the ROI corrects for many missed
			// inerior areas (inlets and narrow bays) when there are clear spots near
			// the surface of the blob (e.g. in cleavage stage embryos).
			run("Enlarge...", "enlarge=10");
			run("Enlarge...", "enlarge=-10");
			//
			roiManager("add");
		};
	// Delete ROIs from the initial step of identifying blob (embryo) centers. 
	// ROIs have to removed in reverse order to avoid changing the indices in
	// problematic ways (if delete ROI 0, 
	for(k=maxk-1;k>=0;k--) { 
			roiManager("select", k);
			roiManager("delete");
		};
	// Measure all ROIs.
	roiManager("measure");
}

// The following takes a selection and fills it with the mean color in the image
// for the whole stack.
macro "add background mask to stack"
{
for(k = 1; k<=nSlices(); k++)
	{
		setSlice(k);
		run("Make Inverse");
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		setColor(mean);
		run("Make Inverse");
		fill();		
	}
}

/* FIRST TRIES:
The following macros are commented out. 
They were initial attempts that failed to accurately identify embryos in some image sequences.
*/
/*
macro "measure embryos 1"
{	maxk = getNumber("Number of times to erode then dilate?", 3);
	run("Set Scale...", "distance=0");
	run("Set Measurements...", "area fit stack display redirect=None decimal=3");
	setAutoThreshold("Moments");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Moments background=Light calculate");
	for(k=1;k<=maxk; k++) 
		{ 	
			run("Erode", "stack");
		};
	for(k=1;k<=maxk; k++) 
		{ 
			run("Dilate", "stack");
		};
	run("Analyze Particles...", "size=10000-100000 circularity=0.10-1.00 display exclude clear include add in_situ stack");
}
*/
/*
macro "measure embryos 2"
{	run("Set Scale...", "distance=0");
	run("Set Measurements...", "area fit stack display redirect=None decimal=3");
	setAutoThreshold("Moments");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Moments background=Light calculate");
	im1 = getImageID()
	run("Duplicate...", "title=A duplicate range=1-" + nSlices());
	im2 = getImageID()
	run("Minimum...", "radius=30 stack");
	run("Maximum...", "radius=60 stack");
	imageCalculator("AND stack", im1, im2);
	run("Analyze Particles...", "size=10000-100000 circularity=0.10-1.00 display exclude clear include add in_situ stack");
}
*/
/*
macro "measure embryos 3"
{	run("Set Scale...", "distance=0");
	run("Set Measurements...", "area fit stack display redirect=None decimal=3");
	im1 = getImageID();
	run("Duplicate...", "title=A duplicate range=1-" + nSlices());
	im2 = getImageID()
	setAutoThreshold("Minimum");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Minimum background=Light calculate");
	run("Maximum...", "radius=10 stack");
	selectImage(im1);
	setAutoThreshold("MaxEntropy");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=MaxEntropy background=Light calculate");
	imageCalculator("AND stack", im1, im2);
	run("Options...", "iterations=1 count=1 pad do=Nothing");
	run("Close-", "stack");
	run("Analyze Particles...", "size=10000-100000 circularity=0.10-1.00 display exclude clear include add in_situ stack");
}
*/
/*
macro "measure embryos 4"
{	run("Set Scale...", "distance=0");
	run("Set Measurements...", "area fit stack display redirect=None decimal=3");
	im1 = getImageID();
	run("Duplicate...", "title=A duplicate range=1-" + nSlices());
	im2 = getImageID()
	run("Gaussian Blur...", "sigma=50 stack");
	setAutoThreshold("Shanbhag stack");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Shanbhag background=Light");
	roiManager("reset");
	run("Analyze Particles...", "size=5000-100000 circularity=0.00-1.00 display exclude clear include add in_situ stack");
	run("Clear Results");
	selectImage(im1);
	setAutoThreshold("MaxEntropy");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=MaxEntropy background=Light calculate");
	run("Options...", "iterations=1 count=1 pad do=Nothing");
	run("Close-", "stack");
	maxk = roiManager("count");
	for(k=0;k<maxk;k++) {
			roiManager("select", k);
			Roi.getBounds(x,y,width,height);
			doWand(x+width/2, y+height/2, 0.0, "4-connected");
			roiManager("add");
		};
	for(k=maxk-1;k>=0;k--) { 
			roiManager("select", k);
			roiManager("delete");
		};
	roiManager("measure");
}
*/
