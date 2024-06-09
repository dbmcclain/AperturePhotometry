# AperturePhotometry

![Screenshot of 3c273](https://github.com/dbmcclain/AperturePhotometry/assets/3160577/d3f79667-e961-4f46-87ce-e4d5935bffa9)

Automated Aperture Photometry engine. In conjunction with LispPlotter (available in a sister repo), it opens a FITS file, extracts the Green (or only) channel, and performs autometed star finding and measuring, showing the results on screen. Move the mouse to stars in either pane and see the measured magnitude. 

Load a FITS file with `(PHOTOM)`. See a listing of found stars with `(REPORT-STARS img)`.

**PHOTOM** _&optional filename => img_
Reads in the G channel, displays the image on screen, finds all the stars with SNR > 5σ, measures them. The found stars are shown in a sister panel with green outlines around each found star. Move the mouse to a star, in either panel, and see the measured magnitude. If you don't provide a filename, the system will bring up a file selection dialog for you to choose your file.

**REPORT-STARS** _img &key sort =>_
Print a report of the found stars for image _img_. By default the report is shown in magnitude order, but you can specify :X or :Y ordering.
