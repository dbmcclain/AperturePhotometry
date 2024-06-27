
(in-package :com.ral.photometry)

;; -----------------------------------------------------------------------------
;; What if we work directly on the star image?

(defun magn (img flux)
  (+ (img-mag-off img)
     (* -2.5f0 (log flux 10f0))))

(defun inv-magn (img mag)
  (expt 10f0 (* -0.4f0 (- mag (img-mag-off img)))))

;; -----------------------------------------------------------------------------
;; Parallelism Support - so easy using Actors...

(defvar *npar*  8)

(defun split-task (farmer-fn start end)
  ;; farm out a long running task
  (flet ((farmed (from to)
           (create
            (lambda (cust)
              (send cust (funcall farmer-fn from to)))
            )))
    (with-recursive-ask
      (ask (apply #'fork
                  (let ((incr (ceiling (- end start) *npar*)))
                    (loop for ix from start below end by incr collect
                            (farmed ix (min end (+ ix incr)))))
                  )))))

(defun split-list-task (farmer-fn lst)
  (let ((qs (make-array *npar*
                        :initial-element nil))
        (qix 0))
    (dolist (item lst)
      (push item (aref qs (mod (incf qix) *npar*))))
    (flet ((farmed (lst)
             (create
              (lambda (cust)
                (send cust (funcall farmer-fn lst)))
              )))
      (with-recursive-ask
        (ask (apply #'fork
                    (loop for q across qs collect (farmed q))
                    )))
      )))

(defmacro merge-splits (expr)
  ;; For use in coalescing the sublists returned from a split task
  `(reduce #'nconc
           (multiple-value-list
            ,expr)))

(defun split-map (mapping-fn lst)
  ;; Apply farmer-fn to every element of lst. Return mapped list.
  ;; Mapping function should return a list item or NIL, since we are
  ;; assembling the result using NCONC.
  (flet ((sublist-handler (lst)
           (loop for item in lst nconc
                   (funcall mapping-fn item))
           ))
  (merge-splits
   (split-list-task #'sublist-handler lst))
  ))

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

(defun format-ra (radeg)
  (let+ ((:mvb (hrs hfrac)  (truncate (/ radeg 15)))
         (:mvb (mins mfrac) (truncate (* 60 hfrac)))
         (:mvb (secs sfrac) (truncate (* 60 mfrac))))
    (format nil "~2,'0D:~2,'0D:~2,'0D.~D"
            hrs mins secs (round sfrac 0.1))))

(defun format-dec (decdeg)
  (let+ ((:mvb (degs dfrac) (truncate (abs decdeg)))
         (:mvb (mins mfrac) (truncate (* 60 dfrac)))
         (secs (round (* 60 mfrac))))
    (format nil "~c~2,'0D:~2,'0D:~2,'0D"
            (if (minusp decdeg) #\- #\+)
            degs mins secs)))

(defun db10 (x)
  (* 10 (log x 10)))

(defun canon-xform (img x y)
  ;; Maybe transform from tilted image to untilted
  (let ((info (img-canon img)))
    (if info
        (let+ ((:mvb (xu yu) (inv-rotxform info x y)))
          (values (round xu)
                  (round yu)
                  (rotxform-arr info)))
      (values (round x)
              (round y)
              (img-arr img))
      )))

(defun measure-location (img x y &key (srch-radius 4))
  (let+ ((:mvb (xx yy arr) (canon-xform img x y)))
    (when (array-in-bounds-p arr yy xx)
      (let+ ((himg         (img-himg img))
             (harr         (img-arr himg))
             (med          (img-med img))
             (s0sq         (img-s0sq img))
             (nsigma       (img-thr img))
             (srch-box     (make-box-of-radius xx yy srch-radius))
             (:mvb (yc xc) (array-max-pos-in-box harr srch-box))
             (:mvb (xcent ycent) (centroid himg xc yc))
             (:mvb (α δ)   (to-radec img xcent ycent))
             (ampl         (aref harr yc xc))
             (pk           (aref arr yc xc)))
        (format t "~%Cursor position ~D, ~D" xx yy)
        (format t "~%Peak position   ~D, ~D (~@D,~@D)" xc yc (- xc xx) (- yc yy))
        (format t "~%ICRS α ~A  δ ~A"
                (format-ra α) (format-dec δ))
        (if (> pk #.(- 65536 256))
            (format t "~%!! Star likely blown !!")
          (when (> pk 32768)
            (format t "~%!! Star possibly in saturation !!")
            ))
        (if (plusp ampl)
            (let+ ((gain   (img-gain img)) ;; e-/ADU
                   (ampl   (* gain ampl))
                   (tnoise (sqrt (+ ampl (* gain gain s0sq))))
                   (thresh (* nsigma tnoise))
                   (snr    (/ ampl tnoise)))
              (format t "~%Peak value ~7,1F (raw ~D)" ampl (round pk))
              (format t "~%Thresh     ~7,1F" thresh)
              (if (>= snr nsigma)
                  (let+ ((star (make-star
                                :x    xcent
                                :y    ycent
                                :pk   pk
                                :ra   α
                                :dec  δ
                                :mag  (magn img ampl)
                                :snr  (db10 snr)
                                :flux ampl
                                :sd   tnoise))
                         (info  (when (img-ncat img)
                                  (find-star-in-cat img star))))
                    (when info
                      (let+ (( (_ dra ddec _ _ cmag) info))
                        (setf (star-catv star) cmag
                              (star-dx star) dra
                              (star-dy star) ddec)
                        ))
                    (print star))
                (format t "~%Failed: Sum below threshold:~%   Mag ≈ ~4,1F  SNR ≈ ~3,1F"
                        (magn img ampl) snr)))
          (format t "~%Failed: Fitted amplitude not positive,~%   central peak ≈ ~4,1F mag"
                  (magn img (- (aref arr yc xc) med)))
          )))))

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
      (array-fill-in-box bits box 1))
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
  (declare (fixnum y x)
           (array single-float arr))
  (let* ((vmax  (bref arr y x))
         (ymax  y)
         (xmax  x))
    (declare (single-float vmax)
             (fixnum ymax xmax))
    (loop for yp fixnum from y to (1+ y) do
            (loop for xp fixnum from (1- x) to (1+ x) do
                    (let ((vval  (bref arr yp xp)))
                      (declare (single-float vval))
                      (when (> vval vmax)
                        (setf vmax vval
                              ymax yp
                              xmax xp))
                      )))
    (if (and (= y ymax)
             (= x xmax))
        (values y x vmax)
      (locate-peak arr ymax xmax))
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

;; ------------------------------------------------------------------
;; Fitted Photometry - Assume an Adaptive Elliptical Bivariate
;; Gaussian Core model + Background level, as per Stetson & DAOPHOT.
;; Major and minor axis sigmas, and angle of major axis of Ellipse,
;; are least-squares fitted values, using brighter stars to guide the
;; fitting. Angle measured with respect to the sensor grid.
#|
;; --------------------------------------------------------------------
;; Least-squares fitting to Gaussian core + BG
;;
;;   I(i,j) = A*G(i,j; σ) + B
;;                                                                          Using Einstein summation conv.
;;   Δ = N * Σ[G(i,j; σ)^2] - (Σ[G(i,j;σ)])^2                             = N G^2_ii - G_ii^2 > 0
;;   A = (N * Σ[I(i,j)*G(i,j;σ)] - Σ[G(i,j;σ)] * Σ[I(i,j)])/Δ             = (N I_ij G_ij - G_ii I_ii)/Δ
;;   B = (Σ[G(i,h;σ)^2] * Σ[I(i,j)] - Σ[G(i,j;σ)] * Σ[I(i,j)*G(i,j;σ)])/Δ = (G^2_ii I_jj - G_kk G_ij I_ij)/Δ
;;
;;   Evan more succinct in Bra-Ket notation:
;;   Σ[G^2] = <G|G>
;;   Σ[G]   = <G|1>, for 1 = a filled array of 1, <1|1> = N.
;;                   Inner prod defined as dot prod between vectors,
;;                   where vectors are row-major access arrays.
;;                   So the entire Gaussian kernel array can be seen as
;;                   a column vector in row-major order.
;;
;;  Treat <X| as a row-vector, and |X> as a column vector.
;;  G is the Gaussian kernel, I is the image of a star.
;;
;;   |I> = A.|G> + B.|1> + |eps>, for noise |eps>
;;
;;   Then:
;;   Δ = N <G|G> - <G|1>^2 > 0,        Proof that Δ > 0? <G|G> and <1|1> are the squared lengths
;;   A = (N <G|I> - <G|1><I|1>)/Δ      of two non-colinear vectors in N-space. Their product is obviously 
;;   B = (<G|G><I|1> - <G|1><G|I>)/Δ   greater than the square of their common projected sub-vector length.
;;
;;   If we define angle Θ from <G|1> = |G|.|1|.cosΘ, where . = ordinary multiplication,
;;   Then Δ = |G|^2.|1|^2.(1 - (cosΘ)^2) = |G|^2.|1|^2.(sinΘ)^2 > 0, unless Θ = 0 or π.
;;   (Also happens to be the squared length of the 2D vector cross prod)
;;
;; Define correlation image <G'| = <G| - <G|1>/N <1|, then <G'|1> = 0, since <1|1> = N.
;; IOW, <G'| is the component of <G| that is orthogonal to <1|. This is Graham-Schmidt orthogonalization.
;; <1| is the everything correlator, a low-pass filter. It causes a sum over all pixels.
;;
;; Then with <G'| we get:
;;
;;   |I> = A'.|G'> + B'.|1> + |eps>, with |eps| << 1 noise, |G'> and |1> orthogonal, <G'|1> = 0.
;;   Ideally, <eps|eps> = |eps|.δ(i,j), <eps|X> ≈ 0 to first order, for any |X>. But this won't really hold
;;   unless the star images are truly Gaussians. 
;;
;;   Δ' = N <G'|G'>
;;   A' = N <G'|I>/Δ' = <G'|I>/<G'|G'>, a measure of similarity or anti-similarity between <G'| and <I|.
;;   B' = <G'|G'><I|1>/Δ' = <I|1>/N = average of image flux spread out over all the pixels.
;;
;;     |G> = |G'> + <G|1>/N.|1>
;;
;;     |I> = A.|G> + B.|1> + |eps>
;;         = A'.|G'> + B'.|1> + |eps>
;;  
;;     <G'|I> = A'.<G'|G'> = A.<G'|G> = A.<G'|G'>
;;                => A = A'
;;     <1|I>  = B'.<1|1> = A.<1|G> + B.<1|1>
;;                => B.N = B'.N - A'.<1|G> => B = B' - A'.<1|G>/N
;;
;;  So, bottom line, we don't give a hoot about B or B'. And A = A'.
;;  So we get the same measured star amplitude using either |G> or
;;  |G'>. But we can do in one convolution with |G'>.
;;
;; --------------------------------------------------------------------------------------------
;;   Noise Variance of measurement in star-free region, given SD for image: (= 1.4826 * MAD)
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
|#

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
    (let+ ((box    (move-box box (- (round x) radius) (- (round y) radius)))
           (star   (extract-subarray arr box))
           (ssum   (vm:total star))
           (kssum  (vm:inner-prod krnl star))
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

(defun acceptable-training-star-p (star)
  ;; Improvements will be based on stars
  ;; with flux: 4000 < adu < 32758
  (let ((z  (star-pk star)))
    (or (< z  4000)
        (> z 32768))
    ))

(defun improve-sigma (img &key fit-args (radius 7))
  (let+ ((stars   (um:take 100
                           (sort 
                            (remove-if #'acceptable-training-star-p
                                       (img-stars img))
                            #'>
                            :key #'star-snr)))
         (nstars  (length stars))
         (stars   (if (> nstars 50)
                      (um:drop (- nstars 50) stars)
                    stars)))
    (format t "~%~D stars selected to guide improvement" (length stars))
    (labels ((quality-sum (vec)
               (let+ ((prof  (make-gaussian-elliptical-fake-star
                              :sigma1 (aref vec 0)
                              :sigma2 (aref vec 1)
                              :theta  (aref vec 2)
                              :radius radius)))
                 (loop for star in stars sum
                         (handler-case
                             (let+ ((:mvb (_ _ resid-img)
                                        (measure-flux (img-arr img)
                                                      (round (star-y star))
                                                      (round (star-x star))
                                                      prof)))
                               (vm:total
                                ;; Squared deviations between core
                                ;; model and star profile,
                                (map-array (lambda (resid-val)
                                             (coerce
                                              (sqr resid-val)
                                              'single-float))
                                           resid-img)))
                           (error ()
                             ;; just ignore - star probably too close to edge of image
                             0
                             ))
                       ))))
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
                (handler-case
                    (let+ ((:mvb (ampl _ resid-img)
                               (measure-flux (img-arr img) (star-y star) (star-x star) prof)))
                      (incf twt ampl)
                      (if resid
                          (map-array-into resid (lambda (r s)
                                                  (coerce (+ r s) 'single-float))
                                          resid resid-img)
                        (setf resid (map-array (lambda (r)
                                                 (coerce r 'single-float))
                                               resid-img))))
                  (error ()
                    ;; just skip the star
                    )))
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
(show-img 'img *saved-img* :level :binary)
(phot-limit *saved-img*)
(report-stars *saved-img*)
(with-seestar
  (setf *saved-img* (photom)))
(with-img *saved-img*
  (measure-stars *saved-img* :thresh 50))
;; Perim = 2*(2r+1) + 2(2r+1-2) = 8r
|#

;; --------------------------------------------------------------------------
;; Parallelized Star Finding

(defun find-stars (ref-img &key (thresh 5) fimg)
  ;; thresh in sigma units
  ;; Find and measure stars in the image.
  (format t "~%Finding stars...")
  (let+ ((krnl        (img-fake-star ref-img))
         (s0sq        (img-s0sq ref-img))
         (nsigma      thresh)
         (:mvb (himg fimg) (make-himg ref-img fimg))
         (harr        (img-arr himg))
         (thr         (* nsigma (sqrt s0sq)))
         (ring-radius (fake-radius krnl))
         (margin      (1+ ring-radius))
         (srch-box    (inset-box (make-box-of-array harr) margin margin))
         ((lf tp rt bt) srch-box)
         (srch-arr    (make-masked-array harr srch-box)))
    (flet ((searcher (start end)
             (let ((srch-box (make-box lf start rt end)))
               (loop for mult in '(200 100 50 25 12 6 1)
                     for mult-thr = (* mult thr)
                     nconc
                       (loop for y from (box-top srch-box) below (box-bottom srch-box)
                             nconc
                               (loop for x from (box-left srch-box) below (box-right srch-box)
                                     nconc
                                       (when (>= (bref srch-arr y x) mult-thr)
                                         (let+ ((:mvb (yc xc pk)   (locate-peak srch-arr y x)))
                                           (when (box-contains-pt-p srch-box xc yc)
                                             (let+ ((gain   (img-gain ref-img)) ;; e-/ADU
                                                    (ampl   (* gain (aref harr yc xc)))
                                                    (tnoise (sqrt (+ ampl (* gain gain s0sq))))
                                                    (snr    (/ ampl tnoise)))
                                               (when (>= snr nsigma)
                                                 ;; everything passed, so clear the
                                                 ;; mask so that we won't find this
                                                 ;; one again.
                                                 (zap-peak srch-arr yc xc thr)
                                                 (let+ ((:mvb (xcent ycent) (centroid himg xc yc)))
                                                   `(,(make-star
                                                       :x    xcent
                                                       :y    ycent
                                                       :pk   pk
                                                       :mag  (magn ref-img ampl)
                                                       :snr  (db10 snr)
                                                       :flux ampl
                                                       :sd   tnoise))
                                                   ))))))
                                     ))))))
      (values
       (merge-splits
        ;; concat farmed results
        (split-task #'searcher tp bt))
       fimg)
      )))

(defun centroid (img xc yc)
  (let+ ((arr   (img-arr img))
         (mass 0)
         (xmass 0)
         (ymass 0))
    (loop for iy from (- yc 3) to (+ yc 3) do
            (loop for ix from (- xc 3) to (+ xc 3) do
                    (let ((v  (max 0 (aref arr iy ix))))
                      (incf mass v)
                      (incf xmass (* (- ix xc) v))
                      (incf ymass (* (- iy yc) v)))))
    (if (plusp mass)
        (values (+ xc (/ xmass mass))
                (+ yc (/ ymass mass)))
      (values xc yc))
    ))
               
(defun conj* (z1 z2)
  (* (conjugate z1) z2))

(defun flexible-fft2d (arr dir)
  (handler-case
      (case dir
        (:fwd (fft2d:fwd arr))
        (t    (fft2d:inv arr)))
    (error () ;; image probably too large to handle in one go...
      (labels ((fft (vec dest dir)
                 (case dir
                   (:fwd  (fft:fwd vec :dest dest))
                   (t     (fft:inv vec :dest dest))))
               (row-handler (dst dir start end)
                 (loop for row from start below end do
                         (let ((rdst (array-row dst row))
                               (rsrc (array-row arr row)))
                           (fft rsrc rdst dir)
                           )))
               (col-handler (src dst dir start end)
                 (let+ ((nrows (array-dimension src 0)))
                   (loop for col from start below end do
                           (let+ ((csrc (array-col src col))
                                  (ans  (fft csrc nil dir)))
                             (loop for row from 0 below nrows do
                                     (setf (aref dst row col) (aref ans row)))
                             )))))
        (let+ (((nrows ncols) (array-dimensions arr))
               (dst           (make-array (array-dimensions arr)
                                          :element-type '(complex single-float))))
         (split-task (um:curry #'row-handler dst dir)     0 nrows)
         (split-task (um:curry #'col-handler dst dst dir) 0 ncols)
         dst
         )))
    ))

(defun make-himg (img fimg)
  ;; Cross-correlate the image array with a Gaussian kernel.
  ;; Return the correlation image.
  ;;
  ;; Convolution? Correlation? With a symmetric kernel there is no
  ;; difference. But after fitting to become a non- cylindrically
  ;; symmetric kernel, it makes a difference.
  ;;
  ;; - Because of the way we construct the kernel, by least-squares
  ;; matching to an average bright star, we need to cross-correlate
  ;; that kernel with the image, *not* convolve. We are constructing a
  ;; matched filter. Hence, correlation and the use of CONJ* instead
  ;; of * below.
  ;;
  ;; Parallelized version - Parallel 2D FFT's of Image and Kernel. If
  ;; FIMG is supplied, it is the Image FFT from a previous run, and
  ;; does not need to be recomputed.
  (let+ ((arr       (img-arr img))
         ( (ht wd)  (array-dimensions arr))
         (box       (make-box-of-array arr))
         (wdx       (um:ceiling-pwr2 wd))
         (htx       (um:ceiling-pwr2 ht))
         (img-proc  (create
                     (lambda (cust)
                       (send cust
                             (or fimg
                                 (let+ ((wrk-arr  (make-image-array htx wdx :initial-element 0f0)))
                                   (implant-subarray wrk-arr arr 0 0)
                                   (flexible-fft2d wrk-arr :fwd))
                                 )))
                     ))
         (krnl-proc (create
                     (lambda (cust)
                       (let+ ((prof      (img-fake-star img))
                              (kradius   (fake-radius prof))
                              (ksum      (fake-ksum prof))
                              (npix      (fake-npix prof))
                              (mnf       (/ ksum npix))
                              (krnl      (map-array (um:rcurry #'- mnf) (fake-krnl prof)))
                              (norm      (vm:inner-prod krnl krnl))
                              (kwrk-arr  (make-image-array htx wdx :initial-element 0f0)))
                         (map-array-into krnl (um:rcurry #'/ norm) krnl)
                         (implant-subarray kwrk-arr krnl 0 0)
                         (send cust (flexible-fft2d (vm:shift kwrk-arr `(,(- kradius) ,(- kradius))) :fwd))
                         ))))
         (:mvb (fimg fkrnl) (with-recursive-ask
                              (ask (fork img-proc krnl-proc))))
         ;; Construct the G' matrix
         ;; FT of a correlation is the conj prod of their FT's
         (ffilt (map-array #'conj* fkrnl fimg))
         (chimg (flexible-fft2d ffilt :inv))
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
      (values himg fimg)
      )))

;; --------------------------------------------------------------------------
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
      (let+ ((:mvb (_xx _yy arr) (canon-xform img xx yy)))
        (when (img-canon img)
          (let+ ((:mvb (ra dec) (to-radec img _xx _yy)))
            (setf (capi:interface-title (capi:element-interface pane))
                  (format nil "Canonical View:  ~A   ~A" (format-ra ra) (format-dec dec)))
            ))
        (when (array-in-bounds-p arr _yy _xx)
          (let+ ((star (find-nearest-star star-db _xx _yy)))
            (when star
              (let ((txt (format nil "~4,1F" (star-mag star))))
                (capi:display-tooltip pane
                                      :x  (+ x 10)
                                      :y  (+ y 10)
                                      :text txt)
                ))
            ))))))
        
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
  (ignore-errors
    (let+ ((x    (round (star-x star)))
           (y    (round (star-y star)))
           (arr  (extract-subarray (img-arr img) (make-box-of-radius (round x) (round y) 20)))
           (xs   (vm:bipolar-framp 41))
           (hs   (array-row arr 20))
           (vs   (array-col arr 20))
           (lo   (* 10 (1- (floor (min (reduce #'min hs) (reduce #'min vs)) 10))))
           (hi   (* 10 (1+ (ceiling (max (reduce #'max hs) (reduce #'max vs)) 10)))))
      (plt:plot 'crosscuts xs hs
                :clear t
                :title  "Image Cross Cuts"
                :xtitle "Pixel Position"
                :ytitle "Image Amplitude [du]"
                :thick  2
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
      )))

(defun star-explain (img)
  ;; mouse left-click on star or region puts up a popup containing all the details...
  ;; If Shift-click, then gives us a chance to recal the magnitude scale.
  (lambda (xc yc x y gspec)
    (declare (ignore x y))  ;; these are prior click positions
    (let* ((ans nil)
           (txt (with-output-to-string (s)
                  (let ((*standard-output* s))
                    (ignore-errors
                      (setf ans (measure-location img (round xc) (round yc))))))
                ))
      (case gspec
        (:SHIFT
         (when ans
           (mp:funcall-async #'recal img ans)))
        (t
         (let+ ((ximg (copy-img img))
                (:mvb (xx yy arr) (canon-xform ximg xc yc))
                (star (or ans
                          (make-star
                           :x  xx
                           :y  yy
                           ))))
           (setf (img-arr ximg) arr)
           (mp:funcall-async #'show-crosscuts ximg star))))
      txt  ;; needs to return the text to display in a popup
      )))

;; ---------------------------------------------------------------

(defun show-img (pane img &key level)
  ;; Show an img o a pane. Optional binarize.
  ;; Default is linear z scaling from Median to Median + 15*MAD ≈ 10σ
  ;; Binarize is hi contrast at 0 and 5σ levels.
  (let+ ((med    (img-med img))
         (mad    (img-mad img))
         (lo     med)
         (hi     (+ med (* 15 mad)))
         (arr    (img-arr img))
         ( (ht wd) (array-dimensions arr))
         (sf     (/ 1250 (max wd ht)))
         (xsize  (round (* wd sf)))
         (ysize  (round (* ht sf))))
    (print `(:scale ,sf))
    (cond
     ((eq level :binary)
      (let ((mad5 (+ med (* (sd-to-mad (img-thr img)) mad))))
        (setf arr  (map-array (lambda (x)
                                (if (>= x mad5) hi lo))
                              arr))
        ))
     ((realp level)
      (let+ ((flux (inv-magn img level)))
        (setf lo (+ med (truncate flux 2))
              hi (+ med (round flux)))))
     )
    (plt:window pane :xsize xsize :ysize ysize :box `(0 0 ,xsize ,ysize))
    (plt:tvscl pane arr
               :clear  t
               ;; :neg t
               :magn   sf
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

(defun improve-stars (img thresh stars0 fimg)
  ;; Using initial find, try to improve the Gaussian core model and redo
  (let ((acceptable (count-if #'acceptable-training-star-p stars0)))
    (cond ((>= acceptable 5)
           ;; Don't bother doing this unless at least 5 stars to base on
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
             (find-stars img :thresh thresh :fimg fimg)
             ))

          (t
           (format t "~%Insufficient number of stars for adapting PSF model")
           stars0)
          )))

(defun measure-stars (img &key (thresh 5))
  (let+ ((:mvb (stars0 fimg) (find-stars img :thresh thresh))
         (stars  (improve-stars img thresh stars0 fimg))
         (himg   (img-himg img)))
    (setf (img-thr   img)  thresh
          (img-stars img)  stars
          (img-thr   himg) thresh
          (img-stars himg) stars)
    (format t "~%~D stars found in image" (length stars))
    #|
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
                     :ytitle "Density"))
    |#
    ;; detections hilighted in green, saturations in red.
    (let* ((arr     (img-arr img))
           (trimmed (remove-if (lambda (star)
                                 (< (aref arr (round (star-y star)) (round (star-x star))) 32768))
                               stars)))
      (plt:with-delayed-update ('himg)
        (show-img 'himg himg)
        (hilight-stars 'himg stars :green)
        (hilight-stars 'himg trimmed :red))
      (plt:with-delayed-update ('stars)
        (show-img 'stars img)
        (hilight-stars 'stars stars :green)
        (hilight-stars 'stars trimmed :red)))
    (values)
    ))

#|
  *core-radius*
(measure-stars *saved-img* :thresh 20)
|#

