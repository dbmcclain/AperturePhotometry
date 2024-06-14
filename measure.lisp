
(in-package :com.ral.photometry)

;; -----------------------------------------------------------------------------
;; What if we work directly on the star image?

(defun magn (img flux)
  (+ (img-mag-off img)
     (* -2.5f0 (log flux 10f0))))

(defun inv-magn (img mag)
  (expt 10f0 (* -0.4f0 (- mag (img-mag-off img)))))

;; ---------------------------------------------------------------
;; For manual checking
;;
;; Move mouse to visible target, even if not in the detections list.
;; Left-click mouse to reveal its approx position. That position is
;; saved in the clipboard and can be pasted into the following
;; MEASURE-LOCATION expression. See why it failed to show up in the
;; detections.
#|
(measure-location *saved-img* 518.  814.)
|#

(defun measure-location (img x y &key (srch-radius 4))
  (let+ ((himg         (img-himg img))
         (arr          (img-arr img))
         (med          (img-med img))
         (harr         (img-arr himg))
         (s0sq         (img-s0sq img))
         (nsigma       (img-thr img))
         (srch-box     (make-box-of-radius x y srch-radius))
         (:mvb (yc xc) (max-array-pos-in-box harr srch-box))
         (ampl         (aref harr yc xc))
         (pk           (aref arr yc xc)))
    (format t "~%Cursor position ~D, ~D" x y)
    (format t "~%Peak position   ~D, ~D (~@D,~@D)" xc yc (- xc x) (- yc y))
    (when (> pk 62000)
      (format t "~%!! Star likely blown !!"))
    (if (plusp ampl)
        (let+ ((tnoise (sqrt (+ ampl s0sq)))
               (thresh (* nsigma tnoise))
               (snr    (/ ampl tnoise)))
          (format t "~%Peak value ~7,1F (raw ~D)" ampl (round pk))
          (format t "~%Thresh     ~7,1F" thresh)
          (if (>= snr nsigma)
              (print (make-star
                      :x    xc
                      :y    yc
                      :mag  (magn img ampl)
                      :snr  snr
                      :core ampl
                      :bg   med ;; for lack of something better to put here
                      :sd   tnoise))
            (format t "~%Failed: Sum below threshold:~%   Mag ≈ ~4,1F  SNR ≈ ~3,1F"
                    (magn img ampl) snr)))
      (format t "~%Failed: Fitted amplitude not positive,~%   central peak ≈ ~4,1F mag"
              (magn img (- (aref arr yc xc) med)))
      )))

;; -------------------------------------------------------------------
;; Automated Star Detection and Measurement

(defstruct (masked-array
            (:constructor %make-masked-array (arr mask)))
  arr mask)

(defun make-masked-array (arr &optional box)
  (let ((bits (make-array (array-dimensions arr)
                          :element-type 'bit
                          :initial-element (if box 0 1))))
    (when box
      (fill-array-in-box bits box 1))
    (%make-masked-array arr bits)
    ))

(defun bref (arr iy ix)
  (if (and
       (array-in-bounds-p (masked-array-arr arr) iy ix)
       (= 1 (aref (masked-array-mask arr) iy ix)))
      (aref (masked-array-arr  arr) iy ix)
    0f0))

(defun set-mask (arr iy ix &optional (val 1))
  (setf (aref (masked-array-mask arr) iy ix) val))

(defun clr-mask (arr iy ix)
  (set-mask arr iy ix 0))

;; --------------------------------------------

(defun locate-peak (arr y x)
  ;; Starting from X, Y, march to the peak of the star image. Since we
  ;; scan images from Left to right, top to bottom, there is no point
  ;; looking behind us in Y. Anything there would have already been zapped away.
  ;;
  (let+ ((zstart  (bref arr y x))
         (zlf     (bref arr y (1+ x)))
         (zrt     (bref arr y (1- x)))
         (zfwd    (bref arr (1+ y) x)))
    (if (> zlf zrt)
        (if (> zlf zfwd)
            (if (> zlf zstart)
                (locate-peak arr y (1+ x))
              (values y x))
          (if (> zfwd zstart)
              (locate-peak arr (1+ y) x)
            (values y x)))
      (if (> zrt zfwd)
          (if (> zrt zstart)
              (locate-peak arr y (1- x))
            (values y x))
        (if (> zfwd zstart)
            (locate-peak arr (1+ y) x)
          (values y x))))
    ))

(defun zap-peak (arr y x thr)
  ;; Clear the mask over all connected pixels, starting from the peak
  ;; position X, Y.  Clearing mask bits makes the pixels ineligible to
  ;; participate in future operations.
  ;;
  ;; A pixel belongs to the star if it can be reached from a linear
  ;; march in Y, and breadth sweeps in X. Spiral patterns will not be
  ;; taken out this way.
  ;;
  (labels ((zapx (y x)
             (clr-mask arr y x)
             (loop for xx from (1+ x)
                   while (>= (bref arr y xx) thr)
                   do
                     (clr-mask arr y xx))
             (loop for xx from (1- x) by -1
                   while (>= (bref arr y xx) thr)
                   do
                     (clr-mask arr y xx)))
           
           (zap (y x)
             (zapx y x)
             (loop for yy from (1+ y)
                   while (>= (bref arr yy x) thr)
                   do
                     (zapx yy x))
             (loop for yy from (1- y) by -1
                   while (>= (bref arr yy x) thr)
                   do
                     (zapx yy x))
               ))
      (zap y x)))

;; ---------------------------------------------------------
;; Aperture Photometry

#|
(defun cresting-p (arr yc xc)
  (let ((z     (aref arr yc xc)))
    (and
     (>= z (aref arr (1- yc) xc))
     (>= z (aref arr (1+ yc) xc))
     (>= z (aref arr yc (1+ xc)))
     (>= z (aref arr yc (1- xc)))
     
     (> z (aref arr (- yc 2) xc))
     (> z (aref arr (+ yc 2) xc))
     
     (> z (aref arr yc (- xc 2)))
     (> z (aref arr yc (+ xc 2)))
     
     (> z (aref arr (1- yc) (1- xc)))
     (> z (aref arr (1- yc) (1+ xc)))
     
     (> z (aref arr (1+ yc) (1- xc)))
     (> z (aref arr (1+ yc) (1+ xc)))
     )))

(defun ring-med (arr yc xc)
  (let* ((ring-box  (make-box-of-radius xc yc *ring-radius*))
         (moat-box  (make-box-of-radius xc yc *moat-radius*))
         (seqs      nil))
    
    (loop for row from (box-top ring-box) below (box-top moat-box) do
            (push (subseq (array-row arr row) (box-left ring-box) (box-right ring-box)) seqs))

    (loop for row from (box-top moat-box) below (box-bottom moat-box) do
            (let ((vec  (array-row arr row)))
              (push (subseq vec (box-left ring-box) (box-left moat-box)) seqs)
              (push (subseq vec (box-right moat-box) (box-right ring-box)) seqs)))

    (loop for row from (box-bottom moat-box) below (box-bottom ring-box) do
            (push (subseq (array-row arr row) (box-left ring-box) (box-right ring-box)) seqs))

    (let* ((vec (apply #'concatenate 'vector seqs))
           (med (vm:median vec))
           (mad (vm:mad vec med)))
      (values med mad)
      )))

#|
(let ((arr (make-array '(23 23)
                       :element-type 'single-float
                       :initial-element 1f0)))
  (inspect (multiple-value-list (ring-med arr 11 11))))
|#

(defun find-stars (ref-img &optional (thresh 5))
  ;; thresh in sigma units
  (let+ ((ref-arr     (img-arr ref-img))
         (core-radius (img-core ref-img))
         (core-width  (1+ (* 2 core-radius)))
         (med         (img-med ref-img))
         (mad         (img-mad ref-img))
         (nsigma      thresh)
         (thr         (+ med (* (sd-to-mad nsigma) mad)))
         (ref-box     (inset-box (make-box-of-array ref-arr) *ring-radius* *ring-radius*))
         (srch-arr    (make-masked-array ref-arr ref-box)))
    (loop for mult in '(100 50 25 12 6 1) nconc
            ;; peel off from bright to faint, to avoid chasing a
            ;; target that is successively eroded before we reach
            ;; it...
            (loop for y from (box-top ref-box) below (box-bottom ref-box) nconc
                    (loop for x from (box-left ref-box) below (box-right ref-box) nconc
                            (when (>= (aref srch-arr y x) (* mult thr))
                              (let+ ((:mvb (yc xc)   (locate-peak srch-arr ref-box y x))
                                     (_              (zap-peak srch-arr ref-box yc xc thr))
                                     (:mvb (med mad) (ring-med ref-arr yc xc))
                                     (box  (make-box-of-radius xc yc core-radius))
                                     (core (- (sum-array ref-arr box)
                                              (* med (box-area box)))))
                                (when (plusp core)
                                  (let* ((poisson (sqrt core))
                                         (sd      (* +mad/sd+ mad core-width))
                                         (noise   (abs (complex sd poisson)))
                                         (snr     (/ core noise)))
                                    (when (>= snr nsigma)
                                      `(,(make-star
                                          :x    xc
                                          :y    yc
                                          :mag  (magn ref-img core)
                                          :snr  snr
                                          :core core
                                          :sd   sd))
                                      )))))
                          )))))
;; ------------------------------------------------------------------
;; Fitted Photometry - assume a Gaussian core model + Background level
;; As per Stetson & DAOPHOT.
|#
;; --------------------------------------------------------------------
;; Least-squares fitting to Gaussian core + BG
;;
;;   I(i,j) = A*G(i,j; σ) + B
;;
;;   Δ = N * Σ[G(i,j; σ)^2] - (Σ[G(i,j;σ)])^2                             = N G^2_ii - G_ii^2 > 0
;;   A = (N * Σ[I(i,j)*G(i,j;σ)] - Σ[G(i,j;σ)] * Σ[I(i,j)])/Δ             = (N I_ij G_ij - G_ii I_ii)/Δ 
;;   B = (Σ[G(i,h;σ)^2] * Σ[I(i,j)] - Σ[G(i,j;σ)] * Σ[I(i,j)*G(i,j;σ)])/Δ = (G^2_ii I_jj - G_kk G_ij I_ij)/Δ
;;
;;   Noise Variance in star-free region, given SD for image: (= 1.4826 * MAD)
;;
;;      S0^2_meas = N SD^2 / Δ
;;
;;   Total Variance of measured star:
;;
;;     S*^2_meas = A + N SD^2 / Δ
;;
;;   SNR of star meas:
;;
;;     SNR = A / Sqrt(S*^2_meas)
;;

(defun measure-flux (arr y x prof)
  ;; Least squares fit of star core to Gaussian profile
  ;; return estimated amplitude and bg.
  ;;
  ;; Unit Delta function at center measures: A = (N - ksum)/Δ < 1, B = (k2sum-ksum)/Δ < 0.
  ;; Unit Tophat at location measures      : A = 0, B = 1
  ;; Gaussian kernel at location measures  : A = 1, B = 0 => engr mag = 0.
  ;;
  (with-accessors ((krnl   fake-krnl)
                   (box    fake-box)
                   (Δ      fake-Δ)
                   (npix   fake-npix)
                   (ksum   fake-ksum)
                   (k2sum  fake-k2sum)
                   (radius fake-radius)) prof
    (let+ ((box    (move-box box (- x radius) (- y radius)))
           (star   (extract-subarray arr box))
           (ssum   (vm:total star))
           (kssum  (vm:total (map-array #'* krnl star)))
           (ampl   (/ (- (* npix kssum)
                         (* ksum ssum))
                      Δ))
           (bg     (/ (- (* k2sum ssum)
                         (* ksum kssum))
                      Δ))
           (resid  (map-array (lambda (sval gval)
                                ;; residual image, in the sense of (Star - Gaussian)
                                (coerce
                                 (- sval
                                    (+ bg
                                       (* ampl gval)))
                                 'single-float))
                              star krnl)))
      (values ampl bg resid)
      )))

#|
;; Unit Delta Function
(let+ ((arr  (make-image-array 15 15 :initial-element 0f0))
       (prof (img-fake-star *saved-img*)))
  (setf (aref arr 7 7) 1f0)
  (measure-flux arr 7 7 prof))

;; Fake-star profile itself
(let+ ((prof (img-fake-star *saved-img*))
       (arr  (map-array (um:rcurry #'coerce 'single-float) (fake-krnl prof))))
  (measure-flux arr 7 7 prof))
 |#

;; -------------------------------------------------------------
;; Noise behavior of matched filtering

(defun s0sq (ref-img)
  ;; prepare estimated noise from measuring in a star-free region
  (let+ ((prof  (img-fake-star ref-img))
         (Δ     (fake-Δ prof))
         (npix  (fake-npix prof))
         (mad   (img-mad ref-img))
         (sd    (* mad +mad/sd+)))
    (/ (* npix sd sd) Δ)))

;; -------------------------------------------------------------

(defun rtod (x)
  (* x (/ 45f0 (atan 1f0))))

(defun improve-sigma (img &key fit-args (radius 7))
  (let+ ((stars   (remove-if (lambda (star)
                               ;; Improvements will be based on stars
                               ;; with magnitudes: 8.0 <= mag <= 11.5
                               (let ((mag (star-mag star)))
                                 (or (< mag  8.0)
                                     (> mag 11.5))
                                 ))
                             (img-stars img))))
    (format t "~%~D stars selected to guide improvement" (length stars))
    (labels ((quality-sum (vec)
               (let+ ((prof  (make-gaussian-elliptical-fake-star
                              :sigma1 (aref vec 0)
                              :sigma2 (aref vec 1)
                              :theta  (aref vec 2)
                              :radius radius)))
                 (loop for star in stars sum
                         (let+ ((:mvb (_ _ resid-img)
                                    (measure-flux (img-arr img) (star-y star) (star-x star) prof)))
                           (vm:total
                            ;; Squared deviations between core
                            ;; model and star profile,
                            (map-array (lambda (resid-val)
                                         (coerce
                                          (sqr resid-val)
                                          'single-float))
                                       resid-img))
                           ))
                 )))
      (let+ ((:mvb (vfit _ niter)
                 (vm:simplex #'quality-sum (apply #'vector fit-args)))
             (sigma1  (aref vfit 0))
             (sigma2  (aref vfit 1))
             (theta   (aref vfit 2))
             (_  (format t "~%~D iters, σ1 = ~6,4F, σ2 = ~6,4F, Θ = ~6,4F (~6,4F deg)"
                         niter sigma1 sigma2 theta (rtod theta)))
             (prof    (make-gaussian-elliptical-fake-star
                       :sigma1 sigma1
                       :sigma2 sigma2
                       :theta  theta
                       :radius radius))
             (twt     0)
             (resid   nil))
        (loop for star in stars do
                (let+ ((:mvb (ampl _ resid-img)
                           (measure-flux (img-arr img) (star-y star) (star-x star) prof)))
                  (incf twt ampl)
                  (if resid
                      (map-array-into resid (lambda (r s)
                                              (coerce (+ r s) 'single-float))
                                      resid resid-img)
                    (setf resid (map-array (lambda (r)
                                             (coerce r 'single-float))
                                           resid-img)))
                  ))
        (map-array-into resid (um:rcurry #'/ (coerce twt 'single-float)) resid)
        ;; Show weighted mean, amplitude normalized, residual image
        ;; in the sense (Star - Gaussian).
        (plt:window 'resid)
        (plt:tvscl 'resid resid
                   :clear t
                   :title "Mean Normalized Resid"
                   :vflip t
                   :magn 16)
        (let ((xs (map 'vector (um:rcurry #'- radius)
                       (vm:framp (+ radius radius 1)))
                  ))
          (plt:plot 'residx xs (array-row resid radius)
                    :clear t
                    :thick 2
                    :title "Mean Normalized Resid"
                    :xtitle "Position [pix]"
                    :ytitle "Amplitude"
                    :yrange '(-0.11 0.11)
                    :legend "X Cut")
          (plt:plot 'residx xs (reverse (array-col resid radius))
                    :thick 2
                    :color :red
                    :legend "Y Cut"))
        (values (list sigma1 sigma2 theta)
                prof
                resid)
        ))))

#|
(improve-sigma *saved-img*)
(setf *fake-star-130* nil)
(with-seestar
  (setf *saved-img* (photom)))
(with-img *saved-img*
  (measure-stars *saved-img*))
;; Perim = 2*(2r+1) + 2(2r+1-2) = 8r
|#

;; --------------------------------------------------------------------------

#|
(defun prep-kernel (ref-img)
  ;; Precompute the geometric quantities based solely on the detailed
  ;; shape of the Gaussian kernel.
  (let+ ((fake-star   (first (img-fake-star ref-img)))
         (ksum        (vm:total fake-star))
         (k2sum       (vm:total (map-array #'* fake-star fake-star)))
         (npix        (array-total-size fake-star))
         (Δ           (- (* npix k2sum)
                         (sqr ksum)))
         (prof        `(,fake-star ,Δ ,npix ,ksum ,k2sum))
         (mad         (img-mad ref-img))
         (sd          (* mad +mad/sd+))      ;; image noise floor SD from measured image MAD
         (s0sq        (/ (* npix sd sd) Δ))) ;; expected noise from a star-free field.
    (values prof s0sq)))
|#
#|
(defun find-stars (ref-img &optional (thresh 5))
  ;; thresh in sigma units
  ;; Find and measure stars in the image.
  (let+ ((ref-arr     (img-arr ref-img))
         (krnl        (img-fake-star ref-img))
         (s0sq        (img-s0sq ref-img))
         (nsigma      thresh)
         (med         (img-med ref-img))
         (mad         (img-mad ref-img))
         (sd          (* mad +mad/sd+))
         (thr         (+ med (* nsigma sd)))
         (ring-radius (fake-radius krnl))
         (srch-box    (inset-box (make-box-of-array ref-arr) ring-radius ring-radius))
         (srch-arr    (make-masked-array ref-arr srch-box)))
    (loop for mult in '(200 100 50 25 12 6 1)
          for mult-thr = (* mult thr)
          nconc
            ;; peel off from bright to faint, to avoid chasing a
            ;; target that is successively eroded before we reach
            ;; it...
            (loop for y from (box-top srch-box) below (box-bottom srch-box) nconc
                    (loop for x from (box-left srch-box) below (box-right srch-box) nconc
                            (when (>= (bref srch-arr y x) mult-thr)
                              (let+ ((:mvb (yc xc)   (locate-peak srch-arr y x))
                                     ;; (_              (zap-peak srch-arr yc xc thr))
                                     (:mvb (ampl bg) (measure-flux ref-arr yc xc krnl)))
                                (when (plusp ampl)
                                  (let+ ((tnoise (sqrt (+ ampl s0sq)))
                                         (snr    (/ ampl tnoise)))
                                    (when (>= snr nsigma)
                                      (zap-peak srch-arr yc xc thr)
                                      `(,(make-star
                                          :x    xc
                                          :y    yc
                                          :mag  (magn ref-img ampl)
                                          :snr  snr
                                          :core ampl
                                          :bg   (- bg med)
                                          :sd   tnoise))
                                      )))))
                          )))))
|#

(defun find-stars (ref-img &optional (thresh 5))
  ;; thresh in sigma units
  ;; Find and measure stars in the image.
  (let+ ((krnl        (img-fake-star ref-img))
         (s0sq        (img-s0sq ref-img))
         (nsigma      thresh)
         (med         (img-med ref-img))
         (himg        (make-himg ref-img))
         (harr        (img-arr himg))
         (thr         (* nsigma (sqrt s0sq)))
         (ring-radius (fake-radius krnl))
         (srch-box    (inset-box (make-box-of-array harr) ring-radius ring-radius))
         (srch-arr    (make-masked-array harr srch-box)))
    (loop for mult in '(200 100 50 25 12 6 1)
          for mult-thr = (* mult thr)
          nconc
            (loop for y from (box-top srch-box) below (box-bottom srch-box) nconc
                    (loop for x from (box-left srch-box) below (box-right srch-box) nconc
                            (when (>= (bref srch-arr y x) mult-thr)
                              (let+ ((:mvb (yc xc)   (locate-peak srch-arr y x))
                                     (ampl   (aref harr yc xc))
                                     (tnoise (sqrt (+ ampl s0sq)))
                                     (snr    (/ ampl tnoise)))
                                (when (>= snr nsigma)
                                  ;; everything passed, so clear the
                                  ;; mask so that we won't find this
                                  ;; one again.
                                  (zap-peak srch-arr yc xc thr)
                                  `(,(make-star
                                      :x    xc
                                      :y    yc
                                      :mag  (magn ref-img ampl)
                                      :snr  snr
                                      :core ampl
                                      :bg   med ;; just to fill it with something, since we don't have bg here
                                      :sd   tnoise))
                                  )))
                          )))))

(defun conj* (z1 z2)
  (* (conjugate z1) z2))

(defun make-himg (img)
  ;; Cross-correlate the image array with a Gaussian kernel.
  ;; Return the correlation image array.
  ;;
  ;; - Because of the way we construct the kernel, by least-squares
  ;; matching to an average bright star, we need to cross-correlate
  ;; that kernal with the image, *not* convolve. We are constructing a
  ;; matched filter. Hence the use of CONJ* instead of * below.
  (let+ ((arr       (img-arr img))
         ( (ht wd)  (array-dimensions arr))
         (box       (make-box-of-array arr))
         (wdx       (um:ceiling-pwr2 wd))
         (htx       (um:ceiling-pwr2 ht))
         (prof      (img-fake-star img))
         (kradius   (fake-radius prof))
         (Δ         (fake-Δ prof))
         (ksum      (fake-ksum prof))
         (npix      (fake-npix prof))
         (mnf       (/ ksum npix))
         (norm      (/ Δ npix))
         (krnl      (map-array (lambda (x)
                                 (/ (- x mnf) norm))
                               (fake-krnl prof)))
         (wrk-arr   (make-image-array htx wdx :initial-element 0f0))
         (kwrk-arr  (make-image-array htx wdx :initial-element 0f0)))
    (implant-subarray wrk-arr arr 0 0)
    (implant-subarray kwrk-arr krnl 0 0)
    (let+ ((fimg  (fft2d:fwd wrk-arr))
           (fkrnl (fft2d:fwd (vm:shift kwrk-arr `(,(- kradius) ,(- kradius)))))
           (ffilt (map-array #'conj* fkrnl fimg))
           (chimg (fft2d:inv ffilt))
           (harr  (make-image-array ht wd)))
      (map-array-into harr #'realpart (extract-subarray chimg box))
      (let* ((med  (vm:median harr))
             (mad  (vm:mad harr med))
             (himg (copy-img img)))
        (setf (img-himg img)  himg
              (img-arr  himg) harr
              (img-med  himg) med
              (img-mad  himg) mad
              (img-himg himg) himg)
        himg
        ))))

#|
(with-seestar
  (setf *saved-img* (photom)))
(measure-stars *saved-img*)
(show-img 'img *saved-img*)
(phot-limit *saved-img*)

(defvar *qref*)
(setf *qref* (img-slice *saved-img* 610 1075 300))
(measure-stars *qref*)
(show-img 'qref *qref*)
(report-stars *qref*)

(setf *sub* (img-slice *saved-img* 502 691 300))
(measure-stars *sub*)
(show-img 'sub *sub*)
(report-stars *sub*)

(make-himg *saved-img*)
(measure-stars *saved-img*)
(with-seestar
  (setf *saved-img* (photom)))

;; From AAVSO 3c273 Chart
;; Me premised on 3c273 = 12.9 mag
  AAVSO     Me       Me-AAVSO   Me2  M32-AAVSO
  10.2      10.0      -0.2      10.5  +0.3
  12.7      12.4      -0.3      12.8  +0.1
  13.5      13.2      -0.3      13.6  +0.1
  11.9      11.6      -0.3      12.0  +0.1
  13.2      13.0      -0.2      13.4  +0.2
  13.0      12.7      -0.3      13.1  +0.1
  12.5      12.2      -0.3      12.6  +0.1
  13.6      13.2      -0.4      13.7  +0.1
  14.2      13.9      -0.3      14.3  +0.1
  12.1      11.7      -0.4      12.1   0.0

;; Me = original meas with core model sigma = 1.3
;; Me2 = adaptively recomputed sigma

(let ((lst '(-0.2 -0.3 -0.3 -0.3 -0.2 -0.3 -0.3 -0.4 -0.3 -0.4))
      (lst '(0.3 0.1 0.1 0.1 0.2 0.1 0.1 0.1 0.1 0.0)))
  (list (vm:mean lst)     ;; => -0.3 => 0.12
        (vm:stdev lst)))  ;; => 0.07 => 0.079
;; So, 3c273 must really be about 13.2 mag right now.

 |#
;; ---------------------------------------------------------------
;; Live history scanning - how to provide a list of detected and
;; measured stars to the GUI running a mouse cursor, which wants to
;; display the magnitude of the star nearest to the cursor?
;;
;; We need the list of stars to remain a list, so that tools like Lisp
;; FIND, POSITION, SUBSEQ will work on sections of the list. But we
;; want to find the relevant location in the list as quickly as
;; possible.
;;
;; So the solution here uses a hash-table to store the head of the
;; sublist for each section of the sky, but all of the stars remain in
;; one grand list.
;;
;; We start by sorting the stars into Y sequence, then for each 10
;; pixel band, we have one hashtable entry that points to the start of
;; the sublist corresponding to that Y coord value and higher.
;; Thereafter, we turn the process over to SELECT-REGION, which uses
;; POSITION and SUBSEQ, on the sublist of stars.
;;
;; Seems to work pretty well...
;;
;; One thing to note, however, is that this matchup with mouse cursor
;; position will work best when the scale of the display is nearly
;; 1:1, or even more magnified.
;;
;; That typically means working on subimage views instead of the whole
;; image.  Squeezing a whole image into a small window implies the
;; loss of fine positioning within the image array..
;;
;; -------------------------------------------------------------

(defvar *selection-radius*  7)  ;; mouse position must be within 7 1:1 pixels of the star
(defvar *index-granularity* 10) ;; record sublist heads for every 10 Y pixels

(defun find-nearest-star (db xc yc)
  ;; db is the Hash Table database pointing to star sublist sections
  ;; of the overall star list.
  (let* ((index  (truncate (max 0 (- yc *selection-radius*)) *index-granularity*))
         (stars  (gethash index db))
         (sel    (select-region stars yc xc *selection-radius*)))
    (when sel
      (flet ((dist (star)
               (abs (complex (- (star-x star) xc)
                             (- (star-y star) yc)))
               ))
        (car
         (reduce (lambda (ans star)
                   (let ((sdist (dist star))
                         (adist (rest ans)))
                     (if (< sdist adist)
                         `(,star . ,sdist)
                       ans)))
                 (rest sel)
                 :initial-value (let ((star (car sel)))
                                  `(,star . ,(dist star)))
                 ))
        ))))

(defun fast-star-db (img)
  ;; Construct a fast lookup for Y positions in the star list.
  ;; Stars are y-ordered.
  ;; Hash table keeps head of sublist of stars for every 10 units of y.
  ;; Using hash table instead of array avoids array bounds checks.
  (let* ((stars  (copy-seq
                  (sort (img-stars img)
                        #'<
                        :key #'star-y)))
         (star-db (make-hash-table :test #'=)))
    (setf (gethash 0 star-db) stars)
    (um:nlet iter ((stars stars)
                   (ctr   0))
      (if (endp stars)
          star-db
        (let* ((star  (car stars))
               (rest  (cdr stars))
               (y     (truncate (star-y star) *index-granularity*)))
          (cond ((> y ctr)
                 (loop for ix from (1+ ctr) to y do
                         (setf (gethash ix star-db) stars))
                 (go-iter rest y))
                (t
                 (go-iter rest ctr))
                ))
        ))
    ))

(defun star-readout (img)
  ;; Provide a callback function to augment the mouse movement handler
  ;; in Plotter. Handler will be called from within the CAPI thread.
  (let ((star-db (fast-star-db img)))
    (lambda (pane x y xx yy)
      ;; x, y are CAPI mouse coords,
      ;; xx, yy are star image array coords.
      (let ((star (find-nearest-star star-db xx yy)))
        (when star
          (let ((txt (format nil "~4,1F" (star-mag star))))
            (capi:display-tooltip pane
                                  :x  (+ x 10)
                                  :y  (+ y 10)
                                  :text txt)
            ))
        ))))
        
;; ---------------------------------------------------------------

(defun recal (img star)
  ;; Give user a chance to recal the magnitude scale.
  (let ((ans (capi:prompt-for-number "Enter desired magnitude"
                                     :initial-value (star-mag star))
             ))
    (when ans
      (let ((adj  (- ans (star-mag star))))
        (incf (img-mag-off img) adj)
        (dolist (s (img-stars img))
          (incf (star-mag s) adj))
        (phot-limit img)
        ))
    ))

(defun show-crosscuts (img star)
  (let+ ((x    (star-x star))
         (y    (star-y star))
         (arr  (extract-subarray (img-arr img) (make-box-of-radius x y 20)))
         (xs   (vm:bipolar-framp 41))
         (hs   (array-row arr 20))
         (vs   (array-col arr 20))
         (lo   (* 10 (1- (floor (min (reduce #'min hs) (reduce #'min vs)) 10))))
         (hi   (* 10 (1+ (ceiling (max (reduce #'max hs) (reduce #'max vs)) 10)))))
    (plt:plot 'crosscuts xs hs
              :clear t
              :title "Image Cross Cuts"
              :thick 2
              :line-type :histo
              :yrange `(,lo ,hi)
              :alpha 0.5
              :legend "X")
    (plt:plot 'crosscuts xs vs
              :thick 2
              :color :red
              :alpha 0.5
              :line-type :histo
              :legend "Y")
    ))

(defun star-explain (img)
  ;; mouse left-click on star or region puts up a popup containing all the details...
  ;; If Shift-click, then gives us a chance to recal the magnitude scale.
  (lambda (xc yc x y gspec)
    (declare (ignore x y))  ;; these are prior click positions
    (let* ((ans nil)
           (txt (with-output-to-string (s)
                  (let ((*standard-output* s))
                    (setf ans (measure-location img (round xc) (round yc)))))
                ))
      (when (star-p ans)
        (case gspec
          (:SHIFT  (mp:funcall-async #'recal img ans))
          (t
           (mp:funcall-async #'show-crosscuts img ans))
          ))
      txt  ;; needs to return the text to display in a popup
      )))

;; ---------------------------------------------------------------

(defun show-img (pane img &key binarize level)
  ;; Show an img o a pane. Optional binarize.
  ;; Default is linear z scaling from Median to Median + 15*MAD ≈ 10σ
  ;; Binarize is hi contrast at 0 and 5σ levels.
  (let+ ((med    (img-med img))
         (mad    (img-mad img))
         (lo     med)
         (hi     (+ med (* 15 mad)))
         (arr    (img-arr img))
         ( (ht wd) (array-dimensions arr))
         (sf     (/ 1000 (max wd ht)))
         (xsize  (round (* wd sf)))
         (ysize  (round (* ht sf))))
    (print `(:scale ,sf))
    (when binarize
      (let ((mad5 (+ med (* (sd-to-mad (img-thr img)) mad))))
        (setf arr  (map-array (lambda (x)
                                (if (>= x mad5) hi lo))
                              arr))
        ))
    (when level
      (let+ ((flux (inv-magn img level)))
        (setf lo (+ med (truncate flux 2))
              hi (+ med (round flux)))))
    (plt:window pane :xsize xsize :ysize ysize)
    (plt:tvscl pane arr
               :clear  t
               ;; :neg t
               :magn   sf
               :flipv  t
               :zrange `(,lo ,hi))
    (plt:set-move-augmentation pane
                               (when (img-stars img)
                                 (star-readout img)))
    (plt:set-click-augmentation pane
                                (when (img-stars img)
                                  (star-explain img)))
    ))

#|
(defvar *saved-img*)         
(show-img 'img *saved-img*)
|#

(defun hilight-stars (img-pane stars color)
  (plt:with-delayed-update (img-pane)
    (dolist (star stars)
      (with-accessors ((ix  star-x)
                       (iy  star-y)) star
        (plt:draw-rect img-pane ix iy 9 9
                       :raw t
                       :border-color color
                       :filled nil
                       :border-thick 1)
        ))))

;; ---------------------------------------------------------------

(defun improve-stars (img thresh stars0)
  ;; Using initial find, try to improve the Gaussian core model and redo
  (if (>= (count-if (lambda (star)
                      (let ((mag (star-mag star)))
                        (and (> mag 8.0)
                             (< mag 11.5))))
                    stars0)
          5)
      ;; Don't bother doing this unless at least 5 stars to base on
      (progn
        (format t "~%~D initial stars" (length stars0))
        (setf (img-stars img) stars0)
        (let+ ((prof       (img-fake-star img))
               (sigma0     (fake-sigma prof))
               (_          (format t "~%initial sigma ~6,4F" sigma0))
               (radius     (fake-radius prof))
               (:mvb (new-sigma new-prof) (improve-sigma img :fit-args sigma0 :radius radius)))
          (format t "~%Improved sigma: ~6,4F" new-sigma)
          (setf (img-fake-star img) new-prof
                (img-s0sq img)      (s0sq img))
          (find-stars img thresh)
          ))
    ;; else
    stars0))

(defun measure-stars (img &key (thresh 5))
  (let+ ((stars0 (find-stars img thresh))
         (stars  (improve-stars img thresh stars0))
         (himg   (img-himg img)))
    (setf (img-thr   img)  thresh
          (img-stars img)  stars
          (img-thr   himg) thresh
          (img-stars himg) stars)
    (format t "~%~D stars found" (length stars))
    (let* ((snrs  (map 'vector #'star-snr stars))
           (pc25  (vm:percentile snrs 25))
           (pc75  (vm:percentile snrs 75)))
      (print (list :snr_25 pc25 :snr_75 pc75))
      (plt:histogram 'maghist snrs
                     :clear t
                     ;; :ylog  t
                     :xlog t
                     :title "SNR Histo"
                     :xtitle "SNR"
                     :ytitle "Density")
      (plt:with-delayed-update ('himg)
        (show-img 'himg himg :binarize nil)
        (hilight-stars 'himg stars :green))
      (plt:with-delayed-update ('stars)
        (show-img 'stars img)
        (hilight-stars 'stars stars :green))
      (values)
      )))

#|
  *core-radius*
(measure-stars *saved-img* :thresh 20)
|#

