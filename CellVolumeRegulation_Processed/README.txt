This folder contains comma-separated-value files of data from the cell volume regulation experiment.
The data were first generated using the ImageJ macro 'MeasureHamEmbryos.ijm' (saved as tab separated values, but ImageJ defaults to the .xls extension); then they were processed using the Python script 'HamSequence.py' (which depends on a 'TrackPoints.py' to group blobs from different
images together based on XY center-center distance). In addition to linking blobs (=ROIs=embryos) in different images, it groups data and removes unanalyzable images/blobs based on the manual classifications in 'CellVolumeRegulation.txt'

Columns include the following:
'IJind' : int, ImageJ measurement indices
'Label': str, ImageJ label (includes {containing folder}:{ROI name}:{image name} (image name includes [0-9]{2}_[0-9]{10} where the first part is the image name from QCam, and the last part is
the time (in seconds) the image was taken. Names were generated when the images were flattened using '48bitRGBto16bitGray.py'
'Area': float, Area of fit ellipse in pixels
'X' : float, X-center, pixels (USED FOR LINKING BLOBS BETWEEN IMAGES)
'Y' : float, y-center, pixels (USED FOR LINKING BLOBS BETWEEN IMAGES)
'Major' : float, fit ellipse major axis, pixels
'Minor' : float, fit ellipse minor axis, pixels
'Angle' : float, fit ellipse angle, degrees
'Feret' : float, maximum Feret's diameter (maximum caliper), pixels (USED FOR VOLUME CALCULATION)
'Slice' : int, slice in ImageJ stack
'FeretX' : float, Beginning of Feret's diameter, pixels (IGNORE: Position varies depending on how ImageJ determines angle)
'FeretY' : float, Beginning of Feret's diameter, pixels (IGNORE)
'FeretAngle' : float, Angle of Feret's max. diameter, degrees, (IGNORE)
'MinFeret' : float, Minimum Feret's diameter (minimum caliper), pixels (USED FOR VOLUME CALCULATION)
'Image' : int, image name (number)
'Time' : int, image acquisition time, seconds
'Volume' : float, volume estimate (4*Pi/3)*(scale*(Feret+MinFeret)/2)^3, µm^3
'Media' : int, which media (FSW or high salt) embryos were in for each image
'Moving' : boolean, whether embryos were manually classified as moving
'ImGroup' : int, groups of images to treat separately in tracking algorithm (TrackPoints.py)
'blobID' : int, keeps track of which measurements belong to the same embryo (e.g. which blobs are linked as the same blob).