
(in-package :com.ral.photometry)

;; -----------------------------------------------------------------------------
;; What if we work directly on the star image?

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

(defun magn (img flux)
  (+ (img-mag-off img)
     (* -2.5f0 (log flux 10f0))))

(defun inv-magn (img mag)
  (expt 10f0 (* -0.4f0 (- mag (img-mag-off img)))))

;; ---------------------------------------------------------------
;; For manual checking
#|
(measure-location *saved-img* 676.  985.)

|#

(defun #1=measure-location (img x y &key (srch-radius 9))
  (let+ ((arr         (img-arr img))
         (med         (img-med img))
         (mad         (img-mad img))
         (nsigma      (img-thr img))
         (core-radius (img-core img))
         (thresh      (+ med (* mad (sd-to-mad nsigma))))
         (srch-box    (make-box-of-radius x y srch-radius)))
    (format t "~%Starting positon ~D, ~D" x y)
    (loop for y from (box-top srch-box) below (box-bottom srch-box) do
            (loop for x from (box-left srch-box) below (box-right srch-box) do
                    (when (>= (aref arr y x) thresh)
                      (return-from #1#
                        (let+ ((:mvb (yc xc) (locate-peak arr srch-box y x))
                               ( _       (format t "~%Peak position ~D, ~D" xc yc))
                               (:mvb (med mad) (ring-med arr yc xc))
                               (box  (make-box-of-radius xc yc core-radius))
                               (core (- (sum-array arr box)
                                        (* med (box-area box)))))
                          (if (plusp core)
                              (let* ((poisson (sqrt core))
                                     (sd      (* +mad/sd+ mad (box-width box)))
                                     (noise   (abs (complex sd poisson)))
                                     (snr     (/ core noise)))
                                (if (>= snr nsigma)
                                    `(,(make-star
                                        :x    xc
                                        :y    yc
                                        :mag  (magn img core)
                                        :snr  snr
                                        :core core
                                        :sd   sd))
                                  "Failed sum threshold"))
                            "Failed core sum not positive")))
                      )))
    "Failed to find star"))

;; -------------------------------------------------------------------
;; Automated Star Detection and Measurement

(defun locate-peak (arr box y x)
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

