
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

(defun #1=measure-location (img x y &key (srch-radius 9))
  (let+ ((krnl        (img-fake-star img))
         (s0sq        (img-s0sq img))
         (arr         (img-arr img))
         (med         (img-med img))
         (nsigma      (img-thr img))
         (thresh      (+ med (* nsigma (sqrt s0sq))))
         (srch-box    (make-box-of-radius x y srch-radius)))
    (format t "~%Starting positon ~D, ~D" x y)
    (loop for y from (box-top srch-box) below (box-bottom srch-box) do
            (loop for x from (box-left srch-box) below (box-right srch-box) do
                    (when (>= (aref arr y x) thresh)
                      (return-from #1#
                        (let+ ((:mvb (yc xc) (locate-peak arr srch-box y x))
                               ( _       (format t "~%Peak position ~D, ~D" xc yc))
                               (:mvb (ampl bg) (measure-flux arr yc xc krnl)))
                          (if (plusp ampl)
                              (let+ ((tnoise (sqrt (+ ampl s0sq)))
                                     (snr    (/ ampl tnoise)))
                                (if (>= snr nsigma)
                                    `(,(make-star
                                        :x    xc
                                        :y    yc
                                        :mag  (magn img ampl)
                                        :snr  snr
                                        :core ampl
                                        :bg   (- bg med)
                                        :sd   tnoise))
                                  (format nil "Failed: Sum below threshold: Mag ≈ ~4,1F  SNR ≈ ~3,1F"
                                          (magn img ampl) snr)))
                            (format nil "Failed: Fitted amplitude not positive, central peak ≈ ~4,1F mag"
                                    (magn img (- (aref arr yc xc) med)))
                            )))
                      )))
    "Failed: Could not find the star"))

;; -------------------------------------------------------------------
;; Automated Star Detection and Measurement

(defun locate-peak (arr box y x)
  ;; Starting from X, Y, march to the peak of the star image. Since we
  ;; scan images from Left to right, top to bottom, there is no point
  ;; looking behind us in Y. Anything there would have already been zapped away.
  ;;
  (let+ ((zstart   (aref arr y x))
         (xp1      (1+ x)))
    (if (and (< xp1 (box-right box))
             (> (aref arr y xp1) zstart))
        (locate-peak arr box y xp1)
      (let ((xm1 (1- x)))
        (if (and (>= xm1 (box-left box))
                 (> (aref arr y xm1) zstart))
            (locate-peak arr box y xm1)
          (let ((yp1 (1+ y)))
            (if (and (< yp1 (box-bottom box))
                     (> (aref arr yp1 x) zstart))
                (locate-peak arr box yp1 x)
              (values y x)
              ))
          ))
      )))

(defun zap-peak (arr box y x thr)
  ;; Set all pixels in the star image to something below detection
  ;; threshold. A pixel belongs to the star if it can be reached from
  ;; a linear march in Y, and breadth sweeps in X. Spiral patterns
  ;; will not be taken out this way.
  ;;
  (let ((repl  (/ thr 2)))
    (labels ((zapx (y x)
               (setf (aref arr y x) repl)
               (loop for xx from (1+ x)
                     while (and (< xx (box-right box))
                                (>= (aref arr y xx) thr))
                     do
                       (setf (aref arr y xx) repl))
               (loop for xx from (1- x) by -1
                     while (and (>= xx (box-left box))
                                (>= (aref arr y xx) thr))
                     do
                       (setf (aref arr y xx) repl)) )
             
             (zap (y x)
               (zapx y x)
               (loop for yy from (1+ y)
                     while (and (< yy (box-bottom box))
                                (>= (aref arr yy x) thr))
                     do
                       (zapx yy x))
               (loop for yy from (1- y) by -1
                     while (and (>= yy (box-top box))
                                (>= (aref arr yy x) thr))
                     do
                       (zapx yy x))
               ))
      (zap y x))))

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
         (srch-arr    (copy-img-array ref-arr))
         (thr         (+ med (* (sd-to-mad nsigma) mad)))
         (ref-box     (inset-box (make-box-of-array ref-arr) *ring-radius* *ring-radius*)))
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
;;   Δ = N * Σ[G(i,j; σ)^2] - (Σ[G(i,j;σ)])^2
;;   A = (N * Σ[I(i,j)*G(i,j;σ)] - Σ[G(i,j;σ)] * Σ[I(i,j)])/Δ
;;   B = (Σ[G(i,h;σ)^2] * Σ[I(i,j)] - Σ[G(i,j;σ)] * Σ[I(i,j)*G(i,j;σ)])/Δ
;; 
(defun measure-flux (arr y x prof)
  ;; Least squares fit of star core to Gaussian profile
  ;; return estimated amplitude and bg.
  ;;
  ;; Unit Delta function at center measures: A = (N - ksum)/Δ < 1, B = (k2sum-ksum)/Δ < 0.
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
                      Δ)))
      (values ampl bg)
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

(defun s0sq (ref-img)
  ;; prepare estimated noise from measuring in a star-free region
  (let+ ((prof  (img-fake-star ref-img))
         (Δ     (fake-Δ prof))
         (npix  (fake-npix prof))
         (mad   (img-mad ref-img))
         (sd    (* mad +mad/sd+)))
    (/ (* npix sd sd) Δ)))

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

(defun find-stars (ref-img &optional (thresh 5))
  ;; thresh in sigma units
  ;; Find and measure stars in the image.
  (let+ ((ref-arr     (img-arr ref-img))
         (krnl        (img-fake-star ref-img))
         (s0sq        (img-s0sq ref-img))
         (med         (img-med ref-img))
         (nsigma      thresh)
         (thr         (coerce (+ med (* nsigma (sqrt s0sq))) 'single-float))
         (srch-arr    (copy-img-array ref-arr))
         (ref-box     (inset-box (make-box-of-array ref-arr) *ring-radius* *ring-radius*)))
    (loop for mult in '(200 100 50 25 12 6 1) nconc
            ;; peel off from bright to faint, to avoid chasing a
            ;; target that is successively eroded before we reach
            ;; it...
            (loop for y from (box-top ref-box) below (box-bottom ref-box) nconc
                    (loop for x from (box-left ref-box) below (box-right ref-box) nconc
                            (when (>= (aref srch-arr y x) (* mult thr))
                              (let+ ((:mvb (yc xc)   (locate-peak srch-arr ref-box y x))
                                     (_              (zap-peak srch-arr ref-box yc xc thr))
                                     (:mvb (ampl bg) (measure-flux ref-arr yc xc krnl)))
                                (when (plusp ampl)
                                  (let+ ((tnoise (sqrt (+ ampl s0sq)))
                                         (snr    (/ ampl tnoise)))
                                    (when (>= snr nsigma)
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

(defun show-img (pane img &key binarize)
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

(defun measure-stars (img &key (thresh 5))
  (let* ((stars (find-stars img thresh)))
    (when stars
      (setf (img-thr   img) thresh
            (img-stars img) stars))
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
      (plt:with-delayed-update ('stars)
        (show-img 'stars img)
        (hilight-stars 'stars stars :green))
      (values)
      )))

#|
  *core-radius*
(measure-stars *saved-img* :thresh 20)
|#

