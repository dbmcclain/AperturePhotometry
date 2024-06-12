;; photom.lisp
;;
;; DM/RAL  2024/06/02 11:51:17 UTC
;; ----------------------------------

(in-package #:com.ral.photometry)

#|
(setf *print-length* 20)
(plt:set-cmap "B-W LINEAR")
|#
;; ----------------------------------

(defvar *fits-segment-length* 2880)

(defvar *fits-hdr*)
(defvar *mag-offset*   #.(+ 21.18 (- 12.9 12.0) (- 12.9 14.96)))

(defvar *fake-star*  nil)
(defvar *fake-mag*   )
(defvar *fake-radius*)
(defvar *fake-stars* )
(defvar *fake-img*   )

(defvar *ring-radius*   11)
(defvar *moat-radius*    7)
(defvar *core-radius*    2)  ;; intentional undersize to prevent overrun with noise
(defvar *erasure-radius* 5)

;; --------------------------------------------

(defvar *fake-star-75*   nil)
(defvar *fake-star-130*  nil)

(defun do-with-seestar (fn)
  (let* ((*core-radius*   2)
         (*mag-offset*    #.(+ 21.18 (- 12.9 12.0) (- 12.9 14.96)))
         (*fake-star*     (or *fake-star-130*
                              (let ((*core-radius* 7))
                                (setf *fake-star-130* (make-gaussian-fake-star :sigma 1.3))))))
    (funcall fn)))

(defun do-with-seestar-channel (fn)
  (let* ((*core-radius*   2)
         (*mag-offset*    #.(+ 9.4f0 (- 12.9f0 0.8617316f0) -0.1f0)) ;; 9.4f0)
         (*fake-star*     (or *fake-star-75*
                              (let ((*core-radius* 3))
                                (setf *fake-star-75* (make-gaussian-fake-star :sigma 0.75))))))
    (funcall fn)))

(defmacro with-seestar (&body body)
  `(do-with-seestar (lambda () ,@body)))

(defmacro with-seestar-channel (&body body)
  `(do-with-seestar-channel (lambda () ,@body)))

;; --------------------------------------------

(defstruct img
  arr med mad hdr
  thr stars
  (mag-off   *mag-offset*)
  (core      *core-radius*)
  (fake-star *fake-star*)
  s0sq)

(defstruct star
  x y mag snr core bg sd)

(defstruct fake
  krnl Δ npix ksum k2sum radius sigma box)


(defun do-with-img (img fn)
  (let* ((*core-radius*  (img-core img))
         (*mag-offset*   (img-mag-off img))
         (*fake-star*    (img-fake-star img)))
    (funcall fn)))

(defmacro with-img (img &body body)
  `(do-with-img ,img (lambda () ,@body)))

;; ---------------------------------

(defun photom (&optional fname (chan :G))
  (um:with-remembered-filename (path "Select FITS File"
                                      :photom fname
                                      :filter "*.fit;*.fts")
    (terpri)
    (format t "~%Channel ~A:" chan)
    (let+ ((img (extract-image path chan))
           (s0sq (s0sq img)))
      (setf (img-s0sq img) s0sq)
      (measure-stars img)
      (show-img 'img img)
      (report-stars img)
      img
      )))

#|
(photom)
 |#
;; ------------------------------------------------------------------------------
;; Subslices of images...
#|
(defun slice (img box)
  (let* ((ximg           (copy-img img))
         (src            (img-arr img))
         (dst            (setf (img-arr ximg)
                               (make-image-array (box-height box) (box-width box)))))
    (loop for srcy from (box-top box) below (box-bottom box)
          for dsty from 0
          do
            (loop for srcx from (box-left box) below (box-right box)
                  for dstx from 0
                  do
                    (setf (aref dst dsty dstx) (aref src srcy srcx))
                  ))
    (let* ((med  (vm:median dst))
           (mad  (vm:mad dst med)))
      (setf (img-med ximg) med
            (img-mad ximg) mad)
      ximg
      )))
|#
#|
(defvar *slice-img* )

(with-seestar
 (let* ((box  (make-box-centered 536 688 540 540))
        (simg (slice *saved-img* box)))
   (show-img 'slice simg)
   (measure-stars simg)
   (report-stars simg :sort :x)
   (with-img simg
     (auto-cal *fake-star*))
   (setf *slice-img* simg)
   (with-img simg
     (show-sub-dets 437 500))
   ))
 |#
;; -------------------------------------------------------------------

(defun show-found-stars (img fake-star-positions)
  (let* ((found-stars (sort  ;; list of found stars sorted in y index
                             (find-stars img 5)
                             #'<
                             :key 'star-y)))
    (plt:with-delayed-update ('found-sky)
      (show-img 'found-sky img)
      (hilight-stars 'found-sky found-stars  :green)
      (hilight-stars 'found-sky fake-star-positions :magenta))
    ))
   
(defun show-fake-sky+hilighting (img fake-star-positions mag)
  (let* ((fake-img (plant-fake-stars img fake-star-positions (img-fake-star img) mag nil)))
    (show-found-stars fake-img fake-star-positions)))
         
(defun fake-sky (img fake-star-positions mag viz)
  (let* ((fake-img (plant-fake-stars img fake-star-positions (img-fake-star img) mag viz))
         (ans      (recover-fake-stars fake-img fake-star-positions 5)))
    (when (consp ans)
      (destructuring-bind (_ _ nbr _ mn _ sd _ snr) ans
        (declare (ignore _))
        (when (and ;; (plusp sd)
                   ;; (> nbr 10)
                   )
          `((,mag ,nbr ,mn ,sd ,snr))
          )))
    ))
#|
(with-seestar-channel
  (let* ((img (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235426.fit"))
         (fake-star-positions (generate-fake-star-positions img t)))
    (show-img 'img img)
    (fake-sky img fake-star-positions 12 t)
    (show-fake-sky+hilighting img fake-star-positions 12)))

(with-seestar-channel
  (let* ((noise-floor 7.5))
    (labels ((fn (mag use-self-noise-p)
               (let* ((flux  (expt 10d0 (* -0.4 (- mag *mag-offset*))))
                      (noise (if use-self-noise-p
                                 (abs (complex noise-floor (sqrt flux)))
                               noise-floor))
                      (snr   (/ flux noise)))
                 snr)))
      (plt:fplot 'plt '(6 18)
                 (um:rcurry #'fn nil)
                 :clear t
                 :title "Poisson Noise vs Background Noise"
                 :xtitle "Star Magnitude [mag]"
                 :ytitle "SNR [σ]"
                 :yrange '(1 10e7)
                 :ylog t
                 :legend (format nil "No SelfNoise, ~Adu BG" noise-floor))
      (plt:fplot 'plt '(6 18)
                 (um:rcurry #'fn t)
                 :thick 2
                 :color :red
                 :legend (format nil "SelfNoise + ~Adu BG" noise-floor))
      (plt:plot 'plt '(0 20) `(5 5)
                :color :gray50
                :thick 3
                :legend "5σ Minimum")
      )))
|#

(defun auto-cal (img)
  (let* ((fake-star-positions  (generate-fake-star-positions img nil))
         (ans  (loop for mag from 9.0 to 18.0 by 0.5 nconc
                       (fake-sky img fake-star-positions mag nil))
               ))
    (when ans
      (plt:plot 'cal-mag (mapcar 'first ans)
                (mapcar 'fourth ans)
                :clear t
                :ylog  t
                :symbol :circle
                :plot-joined t
                :yrange '(0.001 1)
                :title "Accuracy vs Magnitude"
                :xtitle "Magnitude"
                :ytitle "Uncertainty [mag]")
      (plt:plot 'cal-snr (mapcar 'fifth ans)
                (mapcar 'fourth ans)
                :clear t
                :ylog  t
                :xlog  t
                :symbol :circle
                :plot-joined t
                :title "Accuracy vs SNR"
                :xtitle "SNR = Sum(Flux-Med)/σ"
                :ytitle "Uncertainty [mag]")
      (plt:plot 'cal-snrmag
                (mapcar 'first ans)
                (mapcar 'fifth ans)
                :clear t
                :ylog  t
                :symbol :circle
                :plot-joined t
                :title "SNR vs Magnitude"
                :xtitle "Magnitude"
                :ytitle "SNR = Sum(Flux-Med)/σ"))
    ans))

;; -------------------------------------------------------------------
#|
(photom)
(measure-stars *saved-img* :thresh 20)
(report-stars *saved-img*  :sort :mag)
(report-stars *saved-img*  :sort :x)
(show-img 'img *saved-img*)
(show-img 'fake-sky *fake-img*)

;; based on 25-75 percentile in measured SNRs
(gen-fake-star (lambda (snr) (< 380 snr 1680)))
(gen-fake-star (lambda (snr) (< 50 snr 287)))
(make-gaussian-fake-star 0.75)
(make-gaussian-fake-star 1.3)

(let+ ((*core-radius*  3)
       (sigma          0.75)
       ( (img radius mag) (make-gaussian-fake-star :sigma sigma)) ;; core matches 1.3, tails match 1.5
       (dists  (vops:voffset (- radius) (vm:framp (1+ (* 2 radius)))))
       (pk     (aref img radius radius)))
  (plt:tvscl  'stack img
              :clear t
              :zlog  t
              :magn  16)
  (plt:spline 'fake-slice dists (select-row img radius)
              :clear t
              :title "Central Slice of Fake Star"
              :xtitle "Dist from Center [pix]"
              :ytitle "Amplitude"
              :symbol :circle
              :legend "Horizontal")
  (plt:spline 'fake-slice dists (select-col img radius)
              :color :red
              :symbol :circle
              :legend "Vertical")
  (plt:fplot 'fake-slice `(,(- radius) ,radius)
             (lambda (x)
               (* pk (exp (* -0.5 (sqr (/ x sigma))))))
             :color :gray50
             :thick 2
             :legend (format nil "Gaussian ~D" sigma))
  )
    

(generate-fake-star-positions *saved-img* t)
(plant-fake-stars 14.0)
(show-img)
(show-img 'fake-sky *fake-img*)
(recover-fake-stars  5)

;; 100 fake stars planted randomly
;; Isolated by at least 23 pixels.
;;
;; Reported SNR = Sum(Flux-Med)/MAD
;;                                           DeBayered 1080x1920            Stacked G Chan 540x960
;;                Single 10s frame            Stack of 10 frames            Stack of 10 frames
;; Test       Recoverd      #Recov  SNR    Recov Recovered      SNR
;;  mag       mag   sd        %               %    mag  sd 
;; ----    ----------------------------     ---  ------------   -----
;;   9        9.00 ±0.0055   100   2518     100   9.00 ±0.012   12308      100  9.00 0.0011 10819
;;   9.5      9.50  0.0088   100   1589     100   9.50  0.013    7759      100  9.50 0.0019  6834
;;  10       10.00  0.014    100   1004     100  10.00  0.015    4902      100 10.00 0.0029  4311
;;  10.5     10.50  0.022    100    634     100  10.50  0.020    3090      100 10.50 0.0046  2720
;;  11       11.00  0.034    100    400     100  11.00  0.030    1950      100 11.00 0.0073  1716
;;  11.5     11.50  0.055    100    253     100  11.50  0.045    1229      100 11.50 0.012   1082
;;  12       12.00  0.087    100    160     100  12.01  0.070     773      100 12.00 0.018    682
;;  12.5     12.50  0.14     100    102     100  12.51  0.11      486      100 12.50 0.029    430
;;  13       13.00  0.22     100     65     100  13.03  0.18      305      100 13.01 0.046    271
;;  13.5     13.53  0.37      99     41     100  13.55  0.29      191      100 13.51 0.074    170
;;  14       14.04  0.72      56     29     100  14.15  0.56      118      100 14.02 0.12     107
;;  14.5     14.14  0.60       5     25      55  14.67  0.80       79      100 14.54 0.20      67
;;  15                                        9  14.44  0.59       79      100 15.09 0.36      42
;;  15.5                                      1  14.02  ----       70       55 15.56 0.57      28
;;  16                                                                       8 16.17 0.73      15 
;; --------------------------------------------------------------------------------------

;; ---------------------------------------

(defvar *ac-single-chan*
  '((9.0 100 8.998151F0 0.008582744F0 1038.2202F0)
    (9.5 100 9.497104F0 0.013594616F0 655.7461F0)
    (10.0 100 9.99549F0 0.021527859F0 414.4213F0)
    (10.5 100 10.493052F0 0.03408263F0 262.15564F0)
    (11.0 100 10.989496F0 0.053959728F0 166.08245F0)
    (11.5 100 11.484619F0 0.08551538F0 105.46438F0)
    (12.0 100 11.978812F0 0.13609277F0 67.216996F0)
    (12.5 100 12.474662F0 0.2196937F0 43.084507F0)
    (13.0 100 12.984F0 0.3743197F0 27.83852F0)
    (13.5 38 13.378751F0 0.4854413F0 20.546284F0)))
(defvar *ac-sc-stack10*
  '((9.0 100 9.0002575F0 0.002845857F0 3313.8338F0)
    (9.5 100 9.50041F0 0.0045124847F0 2090.6658F0)
    (10.0 100 10.000666F0 0.007156929F0 1318.8988F0)
    (10.5 100 10.501072F0 0.011356414F0 831.94654F0)
    (11.0 100 11.001758F0 0.01803457F0 524.7006F0)
    (11.5 100 11.502925F0 0.02868164F0 330.8415F0)
    (12.0 100 12.004995F0 0.04574219F0 208.52475F0)
    (12.5 100 12.508841F0 0.0733797F0 131.34803F0)
    (13.0 100 13.016477F0 0.119386055F0 82.652794F0)
    (13.5 100 13.533394F0 0.20318781F0 51.9282F0)
    (14.0 99 14.048965F0 0.2811226F0 32.883404F0)
    (14.5 78 14.567439F0 0.45969647F0 21.658888F0)))
(defvar *ac-sees-stack*
  '((9.0 100 9.001291F0 0.0061301254F0 5151.5186F0)
    (9.5 100 9.502063F0 0.009718817F0 3248.1638F0)
    (10.0 100 10.003311F0 0.015412031F0 2047.2265F0)
    (10.5 100 10.505354F0 0.024451752F0 1289.4882F0)
    (11.0 100 11.008755F0 0.038830463F0 811.3867F0)
    (11.5 100 11.514553F0 0.06179052F0 509.72528F0)
    (12.0 100 12.02482F0 0.09878918F0 319.38965F0)
    (12.5 100 12.5440035F0 0.15981102F0 199.29611F0)
    (13.0 100 13.083263F0 0.2677654F0 123.39049F0)
    (13.5 98 13.645895F0 0.44533956F0 76.80098F0)
    (14.0 39 14.082402F0 0.5774779F0 54.1143F0)))
(defvar *ac-pi-stack*
  '((9.0 100 9.000418F0 0.005224646F0 5780.8877F0)
    (9.5 100 9.500671F0 0.008283234F0 3646.5723F0)
    (10.0 100 10.001093F0 0.013135403F0 2299.912F0)
    (10.5 100 10.501808F0 0.020839483F0 1450.2262F0)
    (11.0 100 11.003052F0 0.033092358F0 914.11097F0)
    (11.5 100 11.505306F0 0.052651566F0 575.8449F0)
    (12.0 100 12.009622F0 0.08414309F0 362.41367F0)
    (12.5 100 12.518446F0 0.13595255F0 227.74755F0)
    (13.0 100 13.038283F0 0.22680733F0 142.77895F0)
    (13.5 100 13.594147F0 0.46312556F0 89.3583F0)
    (14.0 74 14.008674F0 0.4752034F0 61.32873F0)))

;; ------------------------------------------------------
;; Half-size Images without deBayering Interpolation

;; Single frame analysis
(with-seestar-channel
  (let ((img 
        ;; (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235118.fit")
        ;; (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235129.fit")
        ;; (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235303.fit")
        ;; (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235415.fit")
        (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235426.fit")
        ;; (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235437.fit")
        ;; (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235528.fit")
        ;; (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235551.fit")
        ;; (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235613.fit")
        ;; (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235630.fit")
        ))
    (setf *ac-single-chan*
          (auto-cal img))))

;; ---------------------------------
;; PixInsight Demosaic Chan 0 10-Stack analysis

(with-seestar-channel
  (let ((img (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Photom/Master/CH 0 Integration.fit")))
    (setf *ac-sc-stack10*
          (auto-cal img))))


;; ------------------------------------------------------
;; Full-size Images with deBayering Interpolation

;; Seestar 10-Stack analysis
(with-seestar
  (let ((img (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar Ss50 Stacked/Stacked_10_Unknown_10.0s_IRCUT_20240508-235640.fit")))
    (setf *ac-sees-stack*
          (auto-cal img))))

;; --------------------------------------------------------
(defvar *saved-img*)

;; PixInsight Drizzle deBayered 10-Stack analysis
(with-seestar
  (let ((img (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Master/drizzle_integration.fit")))
    (setf *saved-img* img)
    (setf *ac-pi-stack*
          (auto-cal img))))
;; --------------------------------------------------------

;; limiting mag based on stdev noise floor and 5σ min signal
F > 5*Sqrt(F + NF^2)
F^2 > 25*(F + NF^2)
F^2 - 25*F - 25*NF^2 > 0

F_min = 12.5 ± Sqrt(156.25 + 25*NF^2)

(defun min-flux (nfsd)
  (+ 12.5
     (sqrt (+ 156.25 (* 25 nfsd nfsd)))
     ))
(min-flux 16.3)
(magn *saved-img* 110)

(plt:histogram 'histo (mapcar #'star-sd (img-stars *saved-img*))
               :clear t)

;; ----------------------------------------------
;; show a histogram of found stars by magnitude
(let* ((arr   (make-array 24
                          :initial-element 0))
       (stars (img-stars *saved-img*)))
  (dolist (star stars)
    (let ((ix (round (+ (star-mag star) 0.5))))
      (incf (aref arr ix))))
  (plt:plot 'plt arr
            :clear t
            :xrange '(7 24)
            :thick 2
            :line-type :histo))
(progn
  (plt:fplot 'plt '(0 100) #'identity :clear t)
  (plt:fplot 'plt '(0 100) (lambda (x) (* 5 (abs (complex (sqrt x) 7.5)))) :color :red))

(with-seestar
  (let* ((mad      5)
         (nf-sigma (* mad +mad/sd+ (1+ (* 2 *core-radius*)))))
    (labels ((inv-mag (x)
               (expt 10.0 (* -0.4 (- x *mag-offset*))))
             (nf (x)
               (* 5 (abs (complex (sqrt (inv-mag x)) nf-sigma)))))
      (plt:fplot 'plt2 '(6 20) (lambda (x) (/ (inv-mag x) (nf x) 0.2))
                 :clear t
                 :thick 2
                 :ylog  t
                 :title "SNR vs Magnitude"
                 :xtitle "Star Magnitude [mag]"
                 :ytitle "SNR [σ]"
                 :legend (format nil "Noise Floor 1σ = ~4,1Fdu" nf-sigma))
      (plt:plot 'plt2 '(6 20) '(5 5)
                :color :gray50
                :thick 3
                :legend "5σ Limit")
                 
      (plt:fplot 'plt '(6 20) #'inv-mag
                 :ylog t
                 :clear t)
      (plt:fplot 'plt '(0 100) #'nf :color :red)
      )))

;; ----------------------------------------------------
;; Measure the variation in magnitude measurements for 3c273
(setf *core-radius* 3)
(make-gaussian-fake-star 0.75)
(show-img)
(report-stars *saved-img* :sort :x)
(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235118.fit")
(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235129.fit")
(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235303.fit")
(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235415.fit")
(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235426.fit")

(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235437.fit")
(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235528.fit")
(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235551.fit")
(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235613.fit")
(photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235630.fit")
(vm:stdev '(12.98 12.52 12.60 12.61 12.82
          12.80  13.07 12.79 13.27 12.70     ))
(vm:mean  '(70  96 100 92 82
                78 58 71 51 81  ))
     
;; ----------------------------------------------------
;; Report (graphically) the results from all the AUTO-CAL scans.
;; Accuracy vs Magnitude
;;
(let* ((single-chan *ac-single-chan*)
       (sc-stack10  *ac-sc-stack10*)
       ;; images with deBayering interpolation
       (sees-stack  *ac-sees-stack*)
       (pi-stack    *ac-pi-stack*))
  (plt:plot 'cal (mapcar 'first sc-stack10)
            (mapcar 'fourth sc-stack10)
            :clear t
            :ylog  t
            :symbol :circle
            :plot-joined t
            :yrange '(0.001 1)
            :title "Accuracy vs Magnitude"
            :xtitle "Magnitude"
            :ytitle "Uncertainty [mag]"
            :legend "PI SC Stack 10")
  (plt:plot 'cal (mapcar 'first single-chan)
            (mapcar 'fourth single-chan)
            :symbol :square
            :color  :red
            :plot-joined t
            :legend "Single Frame")
  (plt:plot 'calb (mapcar 'first pi-stack)
            (mapcar 'fourth pi-stack)
            :clear t
            :ylog  t
            :yrange '(0.001 1)
            :title "Accuracy vs Magnitude"
            :xtitle "Magnitude"
            :ytitle "Uncertainty [mag]"
            :symbol :triangle
            :color  :blue
            :plot-joined t
            :legend "PI DeBayered Stack 10")
  (plt:plot 'calb (mapcar 'first sees-stack)
            (mapcar 'fourth sees-stack)
            :symbol :square
            :color  :gray50
            :plot-joined t
            :legend "Seestar Stack 10") )

;; --------------------------------------------------------------------
;; Do it again for Accuracy vs SNR

(let* ((single-chan *ac-single-chan*)
       (sc-stack10  *ac-sc-stack10*)
       ;; images with deBayering interpolation
       (sees-stack  *ac-sees-stack*)
       (pi-stack    *ac-pi-stack*))
  (plt:plot 'cal (mapcar 'fifth sc-stack10)
            (mapcar 'fourth sc-stack10)
            :clear t
            :ylog  t
            :symbol :circle
            :plot-joined t
            :title "Accuracy vs SNR"
            :xtitle "SNR = Sum(Flux-Med)/σ"
            :ytitle "Uncertainty [mag]"
            :legend "PI SC Stack 10")
  (plt:plot 'cal (mapcar 'fifth pi-stack)
            (mapcar 'fourth pi-stack)
            :symbol :triangle
            :color  :blue
            :plot-joined t
            :legend "PI DeBayered Stack 10")
  (plt:plot 'calb (mapcar 'fifth single-chan)
            (mapcar 'fourth single-chan)
            :symbol :square
            :color  :red
            :plot-joined t
            :legend "Single Frame")
  (plt:plot 'calb (mapcar 'fifth sees-stack)
            (mapcar 'fourth sees-stack)
            :symbol :square
            :color  :gray50
            :plot-joined t
            :legend "Seestar Stack 10") )

;; ------------------------------------------------------------------------------
;; And again..., SNR vs Magnitude
(let* ((single-chan *ac-single-chan*)
       (sc-stack10  *ac-sc-stack10*)
       ;; images with deBayering interpolation
       (sees-stack  *ac-sees-stack*)
       (pi-stack    *ac-pi-stack*))
  (plt:plot 'cal (mapcar 'first sc-stack10)
            (mapcar 'fifth sc-stack10)
            :clear t
            :ylog  t
            :symbol :circle
            :plot-joined t
            :title "SNR vs Magnitude"
            :xtitle "Magnitude"
            :ytitle "SNR = Sum(Flux-Med)/σ"
            :legend "PI SC Stack 10")
  (plt:plot 'cal (mapcar 'first pi-stack)
            (mapcar 'fifth pi-stack)
            :symbol :triangle
            :color  :blue
            :plot-joined t
            :legend "PI DeBayered Stack 10")
  (plt:plot 'cal (mapcar 'first single-chan)
            (mapcar 'fifth single-chan)
            :symbol :square
            :color  :red
            :plot-joined t
            :legend "Single Frame")
  (plt:plot 'cal (mapcar 'first sees-stack)
            (mapcar 'fifth sees-stack)
            :symbol :square
            :color  :gray50
            :plot-joined t
            :legend "Seestar Stack 10") )
;; ----------------------------------------------------------------------------------------
;; Test repeatability - same image, 5 runs with different random fake star placements

(defvar *sav-coll-1*)
(defvar *sav-coll-stack*)

(with-seestar-channel
  (let ((coll-1
         (let ((img (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235426.fit")))
           (loop repeat 16 collect
                   (with-img img
                     (auto-cal img)))))
        (coll-stack
         (let ((img (photom "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Photom/Master/CH 0 Integration.fit")))
           (loop repeat 16 collect
                   (with-img img
                     (auto-cal img))))))

    (setf *sav-coll-1*     coll-1
          *sav-coll-stack* coll-stack)))

(labels ((show (data &rest args)
           (apply 'plt:plot 'cala
                  (mapcar 'first data)
                  (mapcar 'fourth data)
                  ;; :symbol :circle
                  ;; :plot-joined t
                  args)))
  (show (car *sav-coll-stack*)
        :clear t
        :ylog  t
        :yrange '(0.001 1)
        :title "Accuracy vs Magnitude"
        :xtitle "Magnitude"
        :ytitle "Uncertainty [mag]"
        :legend "10x10s Stack"
        )
  (dolist (data (cdr *sav-coll-stack*))
    (show data))
  (show (car *sav-coll-1*)
        :color :red
        :legend "Single 10s Frame")
  (dolist (data (cdr *sav-coll-1*))
    (show data :color :red)))

(labels ((show (data xoff &rest args)
           (let* ((mags  (mapcar 'first data))
                  (xs    (mapcar (um:rcurry #'+ xoff) mags))
                  (mmags (mapcar 'third data)))
             (apply 'plt:plot 'calb
                    xs
                    (mapcar #'- mmags mags)
                    ;; :plot-joined t
                    args))))
  (show (car *sav-coll-stack*) 0.05
        :clear t
        :yrange '(-0.2 0.2)
        :symbol :circle
        :title "Accuracy vs Magnitude"
        :xtitle "Planted Magnitude"
        :ytitle "(Meas - Plant) [mag]"
        :legend "10x10s Stack"
        )
  (dolist (data (cdr *sav-coll-stack*))
    (show data 0.05 :symbol :circle))
  (show (car *sav-coll-1*) -0.05
        :symbol :triangle
        :color :red
        :legend "Single 10s Frame")
  (dolist (data (cdr *sav-coll-1*))
    (show data -0.05 :color :red :symbol :triangle)))

;; ----------------------------------------------------------------------------------------

(let* ((xs   (vm:framp 15))
       (ys   (copy-seq xs)))
  (multiple-value-bind (xc yc slope sigma)
      (linfit:regression xs ys 1)
    (plt:plot 'tst xs ys
              :symbol :circle
              :clear t)
    (plt:fplot 'tst '(0 14)
               (lambda (x)
                 (+ yc (* slope (- x xc))))
               :thick 2)))

(let* ((mags    #(9.04 9.54 10.04 10.54 11.04 11.54 12.22 12.72 13.21 13.71 14.22 14.75 15.36 15.72)) ;;  15.79))
       (mags    #(9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5))
       (sigmas  #(0.0011 0.0019 0.0029 0.0046 0.0073 0.012 0.018 0.029 0.046 0.074 0.12 0.20 0.36 0.57 0.73))
       (wts     #(8700 5491 3466 2188 1382 875 555 352 224 143 92 60 36 28)) ;;  20))
       (mags1   #(9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 )) ;; 14.5))
       (sigmas1 #(0.0055 0.0088 0.014 0.022 0.034 0.055 0.087 0.14 0.22 0.37 0.72 )) ;;  0.60))
       (mags1a  #(9 10 11 12 13))
       (sigmas1a #(0.062 0.070 0.11 0.21 0.46))
       (lsigmas (map 'vector #'log10 sigmas)))
  (multiple-value-bind (xc yc slope σ)
      (linfit:regression mags lsigmas 1)
    (plt:plot 'uncerts
              mags sigmas
              :clear t
              :symbol :circle
              :ylog t
              :thick 2
              :title "Magnitude Accuracy vs Magnitude"
              :xtitle "Magnitude"
              :ytitle "Uncertainty [mag]"
              :legend "10 Image Stack")
    (plt:fplot 'uncerts '(9 16)
               (lambda (x)
                 (expt 10 (+ yc (* slope (- x xc)))))
               :thick 2)
    #|
    (plt:plot 'uncerts
              mags1a sigmas1a
              :alpha 0.3
              :thick 2
              :color :red)
    |#
    (multiple-value-bind (xc yc slope sigma)
        (linfit:regression mags1 (map 'vector #'log10 sigmas1) 1)
      (plt:plot 'uncerts
              mags1 sigmas1
              :symbol :triangle
              :color :red
              :legend "Single Frame")
      (plt:fplot 'uncerts '(9 16)
                 (lambda (x)
                   (expt 10 (+ yc (* slope (- x xc)))))
                 :color :red
                 :thick 2)
      )))

(progn
  (plt:spline 'uncerts
              '(8700 5491 3466 2188 1382 875 555 352 224 143 92 60 36 28 20)
              '(0.0061 0.0094 0.014 0.022 0.034 0.050 0.074 0.11 0.15 0.21 0.30 0.44 0.71 0.93 0.21)
              :clear t
              :symbol :circle
              :xlog t
              :ylog t
              :thick 2
              :title "Magnitude Accuracy vs SNR"
              :xtitle "SNR"
              :ytitle "Uncertainty [mag]"
              :legend "Trusted")
  (plt:plot 'uncerts
              '(20)
              '(0.21 )
              :symbol :circle
              :color  :red
              :legend "Untrustworthy"))

          
(progn
  (plt:spline 'uncerts
              '(22251 8856 3523 1400 556 219 138 95 52 8)
              '(0.0022 0.0055 0.014 0.035 0.087 0.23 0.40 0.55 0.69 0.33)
              :clear t
              :symbol :circle
              :xlog t
              :ylog t
              :thick 2
              :title "Magnitude Accuracy vs SNR"
              :xtitle "SNR"
              :ytitle "Uncertainty [mag]"
              :legend "Trusted")
  (plt:plot 'uncerts
              '(52 8)
              '(0.69 0.33)
              :symbol :circle
              :color  :red
              :legend "Untrustworthy"))

          
|#

(defun report-stars (img &key (sort :mag))
  (let ((sort-key (case sort
                    (:x #'star-x)
                    (:y #'star-y)
                    (t  #'star-mag))
                  )
        (stars   (img-stars img)))
    
    (plt:histogram 'ring-sd (mapcar #'star-sd stars)
                   :clear t
                   :title "Ring Stdev"
                   :xtitle "Ring Stdev [du]"
                   :ytitle "Density")
    (plt:plot 'ring-sd-vs (mapcar #'star-mag stars)
              (mapcar #'star-sd stars)
              :clear t
              :symbol :cross)

    (format t "~%Count  Star Pos       Mag   SNR    CoreSum     BG    SD")
    (format t "~%         X    Y")
    (format t "~%---------------------------------------------------------")
    (loop for star in (sort (img-stars img) #'< :key sort-key)
          for ct from 1
          do
            (with-accessors ((x   star-x)
                             (y   star-y)
                             (mag star-mag)
                             (snr star-snr)
                             (core star-core)
                             (bg   star-bg)
                             (sd   star-sd)) star
              (format t "~%~3D   ~4D ~4D    ~6,2F  ~6,1F  ~7D  ~5D  ~6,1F"
                      ;; crude mag adj based on 3c273
                      ct x y mag (float snr 1.0)
                      (round core) (round bg) sd)
              ))))

#|
(report-stars *saved-img* :sort :mag)
(report-stars *saved-img* :sort :x)
(report-stars *slice-img* :sort :y)
(with-img *slice-img*
  (show-sub-dets 435 500))

|#
;; ---------------------------------------------------------------
;; Slices - sub-images of images

#|
(defun img-slice (img xc yc radius &key binarize)
  (let+ ((new-img   (copy-img img))
         (arr       (img-arr img))
         (med       (img-med img))
         (mad       (img-mad img))
         (img-box   (make-box-of-img img))
         (len       (1+ (* 2 radius)))
         (sub-arr   (make-image-array len len
                                      :initial-element med))
         (src-box   (box-intersection img-box (make-box-of-radius xc yc radius)))
         (sub-box   (make-box-of-array sub-arr))
         (dst-box   (box-intersection sub-box
                                      (move-box src-box (- radius xc) (- radius yc))))
         (wd        (box-width src-box))
         (src-lf    (box-left src-box))
         (dst-lf    (box-left dst-box))
         (hi        (+ med (* 15 mad)))
         (lim       (+ med (* (sd-to-mad (img-thr img)) mad))))
    (loop for src-y from (box-top src-box) below (box-bottom src-box)
          for dst-y from (box-top dst-box)
          do
            (let* ((src-start (array-row-major-index arr src-y src-lf))
                   (src-end   (+ src-start wd))
                   (src-row   (vm:make-overlay-vector arr
                                                      :start src-start
                                                      :end   src-end))
                   (dst-start (array-row-major-index sub-arr dst-y dst-lf))
                   (dst-end   (+ dst-start wd))
                   (dst-row   (vm:make-overlay-vector sub-arr
                                                      :start dst-start
                                                      :end   dst-end)))
              (if binarize
                  (map-into dst-row (lambda (x)
                                      (if (>= x lim)
                                          hi
                                        med))
                            src-row)
                ;; else
                (replace dst-row src-row))))
    (setf (img-arr new-img) sub-arr)
    #|
    (when (equalp sub-box dst-box)
      (let* ((med  (vm:median sub-arr))
             (mad  (vm:mad sub-arr med)))
        (setf (img-med new-img) med
              (img-mad new-img) mad)))
    |#
    new-img
    ))
|#

(defun img-slice (img xc yc radius &key binarize)
  (let ((new-img   (copy-img img))
        (new-arr   (extract-subarray (img-arr img)
                                     (make-box-of-radius xc yc radius))))
    (when binarize
      (let* ((med (img-med new-img))
             (mad (img-mad new-img))
             (sd  (* +mad/sd+ mad))
             (lim (+ med (* sd (img-thr new-img))))
             (vec (vm:make-overlay-vector new-arr)))
        (map-into vec (lambda (x)
                        (if (>= x lim) 1f0 0f0))
                  vec)
        ))
    (setf (img-arr new-img)   new-arr
          (img-stars new-img) nil)
    new-img))
        
(defun show-sub-dets (img xc yc)
  (let* ((sub-img (img-slice img xc yc (img-core img)))
         (arr     (img-arr sub-img))
         (box     (make-box-of-array arr)))
    (loop for row from (box-top box) below (box-bottom box) do
            (format t "~%~{~5D~^ ~}"
                    (map 'list #'round (array-row arr row))))
    ))

(defun show-sub-img (img xc yc radius &key binarize)
  (let* ((sub-img  (img-slice img xc yc radius :binarize binarize))
         (med      (img-med img))
         (hi       (+ med (* 15 (img-mad img))))
         (size     (* 4 (1+ (* 2 radius)))))
    (plt:window 'subimg :xsize size :ysize size)
    #|
    (plt:plot-image 'subimg `(,lf ,rt) `(,bt ,tp) subarr
                    :clear t
                    :zrange '(400 1000))
    |#
    #||#
    (plt:tvscl 'subimg (img-arr sub-img)
               :clear t
               :magn  4
               :flipv t
               :zrange `(,med ,hi))
    #||#
    ))
            
#|
(loop for iy 0 below 4 do
        (loop for ix from 0 below 4 do
              (setf (aref (img-arr *saved-img*) iy ix) 1f6)))

(show-sub-dets 366 292)
(sub-image 366 292)

(sub-image 297 417)

(show-sub-dets *sub* 458 220)
x 459.  y 219.


(show-img 'test (img-slice *saved-img* 100 100 200))

(with-seestar-channel
  (setf *saved-img* (photom)))

(with-seestar
  (setf *saved-img* (photom)))
(defvar *sub*)
(setf *sub* (img-slice *saved-img* 520.  945. 400))
(measure-stars *sub* :thresh 5)
(show-img 'img *saved-img* :binarize t)
(show-img 'img *saved-img* :binarize nil)
(show-img 'sub *sub* :binarize t)
(show-img 'sub *sub* :binarize nil)
(report-stars *sub*)
(measure-stars *saved-img*)
(report-stars *saved-img*)

 |#
#|
 ;; maybe use 2D-FFT's for fast bulk filtering of some sort?
(defun tst ()
  (let+ ((arr (img-arr *sub*))
         (med (img-med *sub*))
         (mad (img-mad *sub*))
         ( (ht wd) (array-dimensions arr))
         (wdx (um:ceiling-pwr2 wd))
         (htx (um:ceiling-pwr2 ht))
         (wrk-arr (make-image-array htx wdx :initial-element 0f0)))
    (implant-subarray wrk-arr arr 
                      (truncate (- htx ht) 2)
                      (truncate (- wdx wd) 2))
    #|
    (let ((vec (vm:make-overlay-vector wrk-arr)))
      (map-into vec (um:rcurry #'- med) vec))
    |#
    (plt:window 'wrk :xsize wdx :ysize htx)
    (plt:tvscl 'wrk (vm:shifth wrk-arr)
               :clear t
               :zrange `(,med ,(+ med (* 15 mad))))

    ;; two ways to compute the same thing, except the first is the
    ;; square of the second

    ;; the first way
    (let+ ((farr  (fft2d:fwd wrk-arr))
           (fcorr (copy-array farr #'conjugate)))
      (map-array-into fcorr #'* fcorr farr)
      (let+ (( (nrows ncols) (array-dimensions fcorr))
             (img  (make-image-array nrows ncols)))
        (map-array-into img #'realpart fcorr)
        (plt:window 'fft :xsize wdx :ysize htx)
        (plt:tvscl 'fft (vm:shifth img)
                   :clear t
                   :zlog t)

        ;; let's see the inverse transform of the autocorrelation spectrum
        (let+ ((xcimg (fft2d:inv fcorr))
               (img   (make-image-array nrows ncols)))
          (map-array-into img #'realpart xcimg)
          (let+ ((med (vm:median img))
                 (mad (vm:mad img med)))
            (plt:window 'ifft :xsize wdx :ysize htx)
            (plt:tvscl 'ifft img
                       :clear t
                       :zrange `(,med ,(+ med (* 15 mad)))))

          ;; the second way
          (let+ ((farr (fft2d:fwd-magnitude wrk-arr))
                 (img  (make-image-array nrows ncols)))
            (plt:window 'fftm :xsize wdx :ysize htx)
            (plt:tvscl 'fftm (vm:shifth farr)
                       :clear t
                       :zlog t))
          )))))
(tst)
|#
(defun tst ()
  (let+ ((arr (img-arr *sub*))
         (med (img-med *sub*))
         ( (ht wd) (array-dimensions arr))
         (wdx (um:ceiling-pwr2 wd))
         (htx (um:ceiling-pwr2 ht))
         (wrk-arr (make-image-array htx wdx :initial-element med))
         (fake (copy-img-array (car (img-fake-star *sub*)) (um:rcurry #'coerce 'single-float)))
         ( (htf wdf) (array-dimensions fake))
         (mnf  (vm:mean fake))
         (gsq  (array*a fake fake))
         (norm (- (vm:total gsq)
                  (/ (sqr (vm:total fake)) (* htf wdf))))
         (krnl (map-array (lambda (x)
                            (/ (- x mnf) norm))
                          fake))
         (krnl-arr (make-image-array htx wdx :initial-element 0f0)))
    (implant-subarray wrk-arr arr 
                      (truncate (- htx ht) 2)
                      (truncate (- wdx wd) 2))
    (implant-subarray krnl-arr krnl
                      (truncate (- htx htf) 2)
                      (truncate (- wdx wdf) 2))
    (let+ ((fimg  (fft2d:fwd wrk-arr)) ;; (vm:shifth wrk-arr)))
           (fkrnl (fft2d:fwd (vm:shifth krnl-arr)))
           (ffilt (map-array #'* fkrnl fimg))
           (chimg (fft2d:inv ffilt))
           (himg  (make-image-array htx wdx)))
      (map-array-into himg #'realpart chimg)
      (let* ((med  (vm:median himg))
             (mad  (vm:mad himg med)))
        (plt:window 'himg :xsize wdx :ysize htx)
        (plt:tvscl 'himg himg ;; (vm:shifth himg)
                   :clear t
                   :flipv t
                   :zrange `(,med ,(+ med (* 30 mad))))
        ))))
#|
(tst)

 |#
