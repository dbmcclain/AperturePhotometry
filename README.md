# AperturePhotometry

_**Automated Aperture Photometry Engine**_ 


![Screenshot of 3c273](https://github.com/dbmcclain/AperturePhotometry/assets/3160577/d3f79667-e961-4f46-87ce-e4d5935bffa9)

The screenshot shows us running 540x540 pixel subframe image of the area around 3C273. The image came from a Seestar 50S telescope which made 10x10s exposures. Those FITS frames were pulled into PixInsight for coalignment, then drizzle integrated 1:1, with Windsorized sigma clipping to remove hot/cold pixels, to produce the stacked image shown here. The cursor in the left panel view is pointing at 3C273. _(The cursor was a red crosshair, but the screen capture momentarily replaced it with the standard pointer arrow.)_ By use of a suitable magnitude offset on the engineering scale, it has been appointed a magnitude of 12.9. All other stars are relative to this scale.

In conjunction with LispPlotter (available in a sister repo), it opens a FITS file, extracts the Green (or only) channel, and performs automated star detection and measurement, showing the results on screen. Move the mouse to stars in either pane and see the measured magnitude next to your cursor. 

Load and analyze a FITS file with `(PHOTOM)`. See a listing of found stars with `(REPORT-STARS img)`.

---
**PHOTOM** _&optional filename channel => img_

Reads in the G channel, by default, displays the image on screen, finds all the stars with SNR > 5σ, measures them. The found stars are shown in a sister panel with green outlines around each found star. Move the mouse to a star, in either panel, and see the measured magnitude next to the mouse cursor. If you don't provide a filename, the system will bring up a file selection dialog for you to choose your file. If the file is from a CFA, it will be demoisaic'd from the Bayer matrix, to pull out one of the green channels, or whichever channel you specified. Channels can be one of 0-4, :R, :G, or :B.

---
**REPORT-STARS** _img &key sort =>_

Print a report of the found stars for image _img_. By default the report is shown in magnitude order, but you can specify :X or :Y ordering.

```
COM.RAL.PHOTOMETRY 7 > (report-stars *sub*)

Count  Star Pos       Mag     SNR  CoreSum  RingSD
         X    Y
--------------------------------------------------
  1     35  145      8.76   601.9  3_62736    22.2
  2    408  217      9.59   410.8  1_69051    16.3
  3     29  146      9.71   387.6  1_51182    31.1
  4     35  151      9.82   369.3  1_37276    29.7
  5    393  708      9.94   349.3  1_22304    16.3
  6     41  147      9.98   342.3  1_18164    31.1
  7    407  223     10.47   273.9    75423    20.8
  8    390  714     10.73   242.4    59271    22.2
  9    711   14     10.81   234.2    55075    14.8
 10    440  238     11.47   172.4    29945    14.8
 11    918  507     11.51   168.9    28887    19.3
 12    627  788     11.53   167.2    28208    16.3
...
```
---

There is a facility for planting fake stars, then reaping them together with real stars. Afterward the list of known fake stars is checked againt the list of harvested stars to see how well the engine performed. Do this for a series of known magnitudes to get an estimate for the quality of measurements being performed, and its repeatability. Also allows for estimating the limiting magnitude in the image, or the probability of detecting a faint star at some magnitude.

This is very much an ongoing work in progress. Feel free to jump in there and try your own ideas!
