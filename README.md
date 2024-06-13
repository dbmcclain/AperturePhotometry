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

**AUTO-CAL** _img => results_

Takes an _img_ and runs an automated series of fake star planting and harvesting, for magnitudes 9 to 18 in steps of 0.5 mag. It collects the results of harvested fakes and presents them as a list. The list contains, for each magnitude level, how many of the 100 planted fakes at that magnitude were recovered, what their mean magnitude measurement reported, what the standard deviation of those measurements was, and what their mean SNR was. _results_ is a list of lists with this information.

---
So what is going on here? I thought that _**DAOPHOT**_ was the canonical standard for Aperture Photometry?

Yes, it probably is. And I should probably have a close look at its source code. Its author, Peter Stetson, is highly regarded, and _**DAOPHOT**_ has been around for at least 20 years. But it was written in Fortran. Perhaps now modernized with C, or C++, or perhaps even has a Java or Python interface. That's fine if you like those things for yourself.

But I now have the time to discover things for myself, and I thrive best in a highly interactive, extensible, programming environment like that offered by Common Lisp. I like to explore ideas right at the keyboard and see immediate results, or not. Lisp lets me do all of that. I love having Lisp Macrology at my fingertips, to bend the core language to my DSL needs. And I love having to think really hard about what the measurement process actually means. You don't get to do any of that if you just use _**DAOPHOT**_.

This is very much an ongoing work in progress. Feel free to jump in there and try your own ideas!

Updates 24-06-11
---
I found that the original algorithm for star finding would cause detection clusters around the skirts of very bright stars. So I then implemented a better peak finder that walks its way to the peak of the star image. But then I was missing the very brightest stars entirely !?

Turns out I was detecting faint spillover from the brightest stars, and this would erase those peripheral areas for future searches to skip over. But then successive scans at the next y index would find another, and stamp out another region ahead of ourselves. By the time we reach the true connected peak of the bright star, it had been eroded by prior skirt finds.

So now we proceed in threshold layers, from bright to faint. The basic threshold is 5σ (or whatever you specify). But we scan the image in stages, 100 basic thresh, 50 basic thresh, 25, 12, 6, and finally 1. That way we won't get confused by bright star spillover - we detect the core of the brightest stars first, then erode them from their true peak pixel outward. This works well.

I also read Peter Stetson's original paper on DAOPHOT where he describes using matched filtering against a presumed Gaussian star profile. That's fine if you also want to discard non-stellar objects. But you are still left with the problem - how do you best find peaks? In his case, the peaks of the convolved image.

But here we don't need to get so fancy. I want to measure everything above threshold that I can see. I don't really care if it is a star or a galaxy. And I also try to solve the problem of finding the peaks, which was glossed over in Peter's paper.

But now... I relent. I see the virtue in using a Gaussian core + BG model. Least-squares solution for star amplitude is good. Does not work well on blown-out bright stars, but nothing else could do any better either. And now we are no longer "Aperture Photometry", but rather a fitted Photometry. 

I decided to switch over because it was apparent that the brightest stars were leaving too much uncounted. It seemed that the aperture core size needed to adapt to how bright the star is. Using the Gaussian fitting is a kind of naturally adaptive solution. So I expanded the Gaussian to a 15x15 overall size, with σ adapted from estimated core sizes of stars. Noise contribution in the wings is not bad, since it is weighted by very low amplitude Gaussian wings at larger distances from the core. Conversely, had I done this expansion in Aperture Photometry, the noise would be killing us at larger distance from the core.

Update 24-06-13
---
We now have nice cursor readback in place. Just move the cursor near a star and it reports the detected star magnitude. Put the cusor anywhere and click the left button to see a complete readback of detection and measurement statistics right there at the cursor position. 

We are implementing the DAOPHOT method for finding stars, using convolution by a depressed Gaussian star-core model. This is a kind of 2D high-pass filter, implemented with 2D FFT's. Then we scan the convolved image, which does a much better job of separating close stars than you can do in the base visual image, looking for image peaks. I use a bit-map mask to block out already detected stars. Being a high-pass filter, it blocks out low spatial frequency variations, to first order, from image gradients and nebular clouds.

What I'm seeing for results is nothing short of astonishing! If you can see it, then it probably got detected and measured. Even many you can't readily see. But, in pixel peeping, some faint smudges simply don't reach detection threshold levels, as shown by our mouse-click statistics feedback. 

I took a stack on M5 (the globular cluster in Hercules). Even in the most crowded regions of the cluster, the detection algorithm works wonders. We are able to separate and measure very close pairs of stars - too close to even see them clearly separated. Image stacks, of 5 mins duration, coming off a little 2-inch telescope, are showing statistically significant detections of stars down to magnutude 16.2. The camera/telescope combination has a best possible limit of 16.7 mag, if you could drive the noise to zero. Single 10s exposure images show typical limits around 14.5 mag. 1 minute stacks show 15.0 mag. The dynamic range of the system shows that stars brighter than about 7.5 mag will saturate and blow the cores of the stars. So we have an effective dynamic range of about 7.5 mag in a 1 min stack, or about 1000:1 for the 16-bit camera. (The camera sensor is most likely really only 12-14 bits.)

However, what we have implemented matches just the star detection front end of DAOPHOT. The rest of DAOPHOT performs a much more rigorous teasing apart of crowded star groups and does multi-star simultaneous fitting against semi-empirical model PSF's. It is indeed the Cadillac system for study of crowded star fields. For the moment I have no intention of proceeding any further. This was intended to be a quick and dirty system for scanning, mostly amateur, images.

Checking against AAVSO star charts shows that we are already getting better accuracy than the charts report. Our measured magnitudes are ±0.07 mag of AAVSO numbers. (And as a result, we should be using 13.2 mag for 3C273).

![M27](https://github.com/dbmcclain/AperturePhotometry/assets/3160577/4f252d35-dd09-414f-be67-b08161f66902)

![M5](https://github.com/dbmcclain/AperturePhotometry/assets/3160577/c2099e70-ba96-49b9-9a5d-570c5addad70)

