# AperturePhotometry

![Screenshot of 3c273](https://github.com/dbmcclain/AperturePhotometry/assets/3160577/d3f79667-e961-4f46-87ce-e4d5935bffa9)

The screenshot shows us running an image of the area around 3C273. The cursor in the left panel image is pointing at 3C273. It has been appointed a magnitude of 12.9. All other stars are relative to this scale.

Automated Aperture Photometry engine. In conjunction with LispPlotter (available in a sister repo), it opens a FITS file, extracts the Green (or only) channel, and performs automated star detection and measurement, showing the results on screen. Move the mouse to stars in either pane and see the measured magnitude next to your cursor. 

Load a FITS file with `(PHOTOM)`. See a listing of found stars with `(REPORT-STARS img)`.

**PHOTOM** _&optional filename channel => img_

Reads in the G channel, by default, displays the image on screen, finds all the stars with SNR > 5σ, measures them. The found stars are shown in a sister panel with green outlines around each found star. Move the mouse to a star, in either panel, and see the measured magnitude. If you don't provide a filename, the system will bring up a file selection dialog for you to choose your file. If the file is from a CFA, it will be demoisaic'd from the Bayer matrix, to pull out one of the green channels, or whichever channel you specified.

**REPORT-STARS** _img &key sort =>_

Print a report of the found stars for image _img_. By default the report is shown in magnitude order, but you can specify :X or :Y ordering.

There is a facility for planting fake stars, then reaping them together with real stars. Afterward the list of known fake stars is checked againt the list of harvested stars to see how well the engine performed. Do this for a series of known magnitudes to get an estimate for the quality of measurements being performed, and its repeatability. Also allows for estimating the limiting magnitude in the image, or the probability of detecting a faint star at some magnitude.

This is very much an ongoing work in progress. Feel free to jump in there and try your own ideas!
