# AperturePhotometry

_**Automated Aperture Photometry Engine**_ 


![Screenshot of 3c273](https://github.com/dbmcclain/AperturePhotometry/assets/3160577/f7bb2c68-9346-4a18-8934-4a6ca91ef306)

The screenshot shows us running a 540x540 pixel subframe image of the area around 3C273. The image came from a Seestar 50S telescope which made 10x10s exposures. Those FITS frames were pulled into PixInsight for coalignment, then drizzle integrated 1:1, with Windsorized sigma clipping to remove hot/cold pixels, to produce the stacked image shown here. The cursor in the left panel view is pointing at 3C273. _(The cursor was a red crosshair, but the screen capture momentarily replaced it with the standard pointer arrow.)_ By use of a suitable magnitude offset on the engineering scale, it has been appointed a magnitude of 12.9. All other stars measure relative to this scale.

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

**SHOW-IMG** _pane img &key binarize =>_

Show an image, _img_, in any pane of your choosing - just give it a name. BINARIZE will produce a hard contrast image for all stars above the 5σ level. By default the scaling is linear from image Median to (Median + 15*MAD).

`(show-img 'img *saved-img* :binarize t)`

---
An _**img**_ is a data structure that contains the image array shown on screen, the list of stars detected and measured, some info about the image, like its overall Median and MAD, as well as the σ-limit used during detection (by default 5σ). The SNR of a star measurement is based on the ratio of its measured summed Flux in a core region covering the star, and the RMS sum of background noise levels from nearby image MAD measurement, and the Poisson self-noise of the starlight flux. Bright stars are self-limiting, with their SNR growing as Sqrt(Flux). Faint stars are limited by the sky background, with SNR growing in direct proportion to their measured Flux.

---
There is a facility for planting fake stars, then reaping them together with real stars. Afterward the list of known fake stars is checked against the list of harvested stars to see how well the engine performed. Do this for a series of known magnitudes to get an estimate for the quality of measurements being performed, and its repeatability. Also allows for estimating the limiting magnitude in the image, or the probability of detecting a faint star at some magnitude.

---
So what is going on here? I thought that DAOPHOT was the canonical standard for Aperture Photometry?

Yes, it probably is. And I should probably have a close look at its source code. Its author, Peter Stetson, is highly regarded, and DAOPHOT has been around for at least 20 years. But it was written in Fortran. Perhaps now modernized with C, or C++, or perhaps even has a Java or Python interface. That's fine if you like those things for yourself.

But I now have the time to discover things for myself, and I thrive best in a highly interactive, extensible, programming environment like that offered by Common Lisp. I like to explore ideas right at the keyboard and see immediate results, or not. Lisp lets me do all of that. I love having Lisp Macrology at my fingertips, to bend the core language to my DSL needs. And I love having to think really hard about what the measurement process actually means. You don't get to do any of that if you just use DAOPHOT.

This is very much an ongoing work in progress. Feel free to jump in there and try your own ideas!
