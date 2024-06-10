
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
(measure-location *sub* 694.  226.)
|#

(defun measure-location (img xc yc)
  (let+ ((arr         (img-arr img))
         (core-radius (img-core img))
         (core-width  (1+ (* 2 core-radius)))
         (sub-arr     (extract-subarray arr (make-box-of-radius xc yc core-radius)))
         (ht          (array-dimension sub-arr 0))
         (row-sum     (let ((acc (copy-seq (array-row sub-arr 0))))
                        (loop for row from 1 below ht do
                                (map-into acc #'+ acc (array-row sub-arr row)))
                        acc))
         (ramp        (vm:framp ht))
         (xcent       (round
                       (+ xc (- core-radius)
                          (/ (vm:total (map 'vector #'* ramp row-sum))
                             (vm:total row-sum)))))
         (col-sum     (let ((acc (array-col sub-arr 0)))
                        (loop for col from 1 below ht do
                                (map-into acc #'+ acc (array-col sub-arr col)))
                        acc))
         (ycent       (round
                       (+ yc (- core-radius)
                         (/ (vm:total (map 'vector #'* ramp col-sum))
                            (vm:total col-sum))))))
    (format t "~%Starting positon ~D,~D" xc yc)
    (format t "~%Centroid position ~D,~D" xcent ycent)

    (let+ ((med         (img-med img))
           (mad         (img-mad img))
           (thresh      (img-thr img)) ;; in SD units
           (thr         (+ med (* (sd-to-mad thresh) mad))))
      (if (< (aref arr ycent xcent) thr)
          "Failed initial threshold"
        ;; else
        (if (not (cresting-p arr ycent xcent))
            "Failed to crest"
          ;; else
          (let+ ((:mvb (med mad) (ring-med arr ycent xcent))
                 (box  (make-box-of-radius xcent ycent core-radius))
                 (core (- (sum-array arr box)
                          (* med (box-area box)))))
            (if (plusp core)
                (let* ((poisson (sqrt core))
                       (sd      (* +mad/sd+ mad core-width))
                       (noise   (abs (complex sd poisson)))
                       (snr     (/ core noise)))
                  (if (>= snr thresh)
                      (make-star
                       :x    xcent
                       :y    ycent
                       :mag  (magn img core)
                       :snr  snr
                       :core core
                       :sd   sd)
                    ;; else
                    "Failed second threshold on summed flux."
                    ))
              ;; else
              "Summed Flux not positive."))
          )))
    ))

;; -------------------------------------------------------------------
;; Automated Star Detection and Measurement

(defun find-stars (ref-img &optional (thresh 5))
  ;; thresh in sigma units
  (let+ ((ref-arr     (img-arr ref-img))
         (core-radius (img-core ref-img))
         (core-width  (1+ (* 2 core-radius)))
         (med         (img-med ref-img))
         (mad         (img-mad ref-img))
         (srch-arr    (copy-array ref-arr))
         (thr         (+ med (* (sd-to-mad thresh) mad)))
         (ref-box     (inset-box (make-box-of-array ref-arr) *ring-radius* *ring-radius*)))
    (loop for yc from (box-top ref-box) below (box-bottom ref-box) nconc
            (loop for xc from (box-left ref-box) below (box-right ref-box) nconc
                    (when (and (>= (aref srch-arr yc xc) thr)
                               (cresting-p srch-arr yc xc))
                      (let+ ((:mvb (med mad) (ring-med ref-arr yc xc))
                             (box  (make-box-of-radius xc yc core-radius))
                             (core (- (sum-array ref-arr box)
                                      (* med (box-area box)))))
                          (when (plusp core)
                            (let* ((poisson (sqrt core))
                                   (sd      (* +mad/sd+ mad core-width))
                                   (noise   (abs (complex sd poisson)))
                                   (snr     (/ core noise)))
                              (when (>= snr thresh)
                                (fill-array srch-arr box med)
                                `(,(make-star
                                    :x    xc
                                    :y    yc
                                    :mag  (magn ref-img core)
                                    :snr  snr
                                    :core core
                                    :sd   sd))
                                )))))
                  ))))

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

