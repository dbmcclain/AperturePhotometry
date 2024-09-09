
(in-package :com.ral.photometry)

;; -------------------------------------------------------

(defun query-header (hdr key)
  (let* ((keylen (length key))
         (npad   (- 8 keylen))
         (keystr (concatenate 'string key (make-string npad :initial-element #\Space) "="))
         (ans    (find keystr hdr
                       :test #'string-equal
                       :key  (um:rcurry #'subseq 0 9))) )
    (when ans
      (string-left-trim '(#\Space #\Tab) (subseq ans 9)))
    ))

(defun nquery-header (hdr key)
  (let ((ans (query-header hdr key)))
    (when ans
      (read-from-string ans))))


#|
(query-header (img-hdr *saved-img*) "NAXIS")
|#

(defun extract-image (fname &optional (chan 0))
  ;; Strip out the first G channel from the Bayer CFA (Seestar S50)
  (let* ((data  (hcl:file-binary-bytes fname)))
    ;; extract the header
    (multiple-value-bind (off hdr)
        (um:nlet iter ((ix   0)
                       (pos  0)
                       (ans  nil))
          (let* ((line (map 'string #'code-char (subseq data pos (+ pos 80)))))
            (if (string-equal "END     " (subseq line 0 8))
                (values (* *fits-segment-length* (ceiling (+ pos 80) *fits-segment-length*))
                        (nreverse ans))
              (go-iter (1+ ix) (+ pos 80) (cons line ans))
              )))

      ;; extract an image
      (let* ((naxis  (nquery-header hdr "naxis"))
             (bitpix (nquery-header hdr "bitpix"))
             (naxis1 (nquery-header hdr "naxis1"))
             (naxis2 (nquery-header hdr "naxis2"))
             (naxis3 (or (nquery-header hdr "naxis3") 1))
             (bzero  (truncate (or (nquery-header hdr "bzero")  0)))
             (bscale (truncate (or (nquery-header hdr "bscale") 1)))
             (bayer  (nquery-header hdr "bayerpat")) ;; should return as a quoted symbol, e.g., 'GRBG
             ;; (gain   (or (nquery-header hdr "gain") 1))
             (instr  (query-header hdr "INSTRUME"))
             (is-seestar (search "Seestar" instr       ;; else Vespara II?
                                 :test #'equalp)) 
             )

        (when (and naxis
                   bitpix
                   naxis1
                   naxis2
                   (integerp naxis)
                   (integerp bitpix)
                   (integerp naxis1)
                   (integerp naxis2)
                   (integerp naxis3))
          (let ((img (case naxis
                       (2
                        (let ((reader-fn (reader-fn data off bitpix bzero bscale)))
                          (cond (bayer
                                 (extract-cfa chan (cadr bayer) bitpix naxis1 naxis2 reader-fn))
                                (t
                                 (extract-mono bitpix naxis1 naxis2 reader-fn)
                                 ))))
                       (3
                        (let* ((plane-size (* naxis1 naxis2))
                               (chan       (mod
                                            (cond
                                             ((integerp chan) chan)
                                             ((eql chan :R)   0)
                                             ((eql chan :G)   1)
                                             (t               2))
                                            naxis3))
                               (pix-wd     (truncate bitpix 8))
                               (plane-off  (+ off (* chan plane-size pix-wd)))
                               (reader-fn  (reader-fn data plane-off bitpix bzero bscale)))
                          (extract-mono bitpix naxis1 naxis2 reader-fn)
                          ))
                       )))

            ;; show image statistics
            (let* ((med  (vm:median img))
                   (mad  (vm:mad img med))
                   (sd   (* +sd/mad+ mad))
                   (max  (reduce #'max (vec img))))
              (print (list :med med :mad mad :sigma sd :max max))

              #|
              (plt:histogram 'img-histo
                             (map 'vector (lambda (x)
                                            (clip
                                             (/ (- x med) mad)
                                             -5 15))
                                  (vm:make-overlay-vector img))
                             :clear t
                             :title "BG Image Statistics"
                             :xtitle "(Value - Median) [MAD]"
                             :ytitle "Density")
              |#
              (plt:histogram 'img-histox
                             (map 'vector (lambda (x)
                                            (clip
                                             (/ (- x med) sd)
                                             -3 5))
                                  (vm:make-overlay-vector img))
                             :clear t
                             :title "BG Image Statistics"
                             :xtitle "(Value - Median) [σ]"
                             :ytitle "Density")

              (make-img
               :arr  img
               :med  med
               :mad  mad
               :hdr  hdr
               :gain (if is-seestar ;; e-/ADU
                         1.1   ;; IMX462 Seestar S50
                       0.6)   ;; IMX585 Vespera II
               :is-see is-seestar)
              ))))
      )))

(defun reader-fn (data off bitpix bzero bscale)
  ;; normalize input data to (0.0,1.0)
  (case bitpix
    (8  (lambda (pos)
          (let* ((dpos  (+ off pos))
                 (val   (aref data dpos)))
            (float (logand (+ (* bscale val) bzero) #xFF)))
          ))
    (16 (lambda (pos)
          (let* ((dpos  (+ pos off))
                 (val   (+ (ash (aref data dpos) 8)
                           (aref data (1+ dpos)))))
            (float (logand (+ (* bscale val) bzero) #xFFFF)))
          ))
    ))

(defun extract-mono (bitpix naxis1 naxis2 reader-fn)
  (let* ((wd        naxis1)
         (ht        naxis2)
         (chan-wd   (truncate bitpix 8))
         (img       (make-image-array ht wd)))
    (loop for ix from 0 below (* wd ht)
          for pos from 0 by chan-wd
          do
            (setf (row-major-aref img ix)
                  (funcall reader-fn pos)))
    img))

(defun extract-cfa (chan bayer bitpix naxis1 naxis2 reader-fn)
  (let* ((cfa-wd naxis1)
         (cfa-ht naxis2)
         (chan   (cond
                  ((integerp chan) (mod chan 4))
                  ((eql chan :R)
                   (position #\R (string bayer) :test #'char-equal))
                  ((eql chan :G)
                   (position #\G (string bayer) :test #'char-equal))
                  ((eql chan :B)
                   (position #\B (string bayer) :test #'char-equal))
                  ))
         (ht        (truncate cfa-ht 2))
         (wd        (truncate cfa-wd 2))
         (chan-wd   (truncate bitpix 8))
         (line-wd   (* chan-wd cfa-wd))
         (pix-step  (* 2 chan-wd))
         (line-step (* 2 line-wd))
         (img       (make-image-array ht wd))
         (base-pos  (+ (* (logand chan 1) chan-wd)     ;; offset to odd number chan
                       (* (ash chan -1)   line-wd))))  ;; offset to chan (0,1) or (2,3)
    (loop for row from 0 below ht
          for pos from base-pos by line-step
          do
            (loop for col from 0 below wd
                  for cfa-pos from pos by pix-step
                  do
                    (setf (aref img row col)
                          (funcall reader-fn cfa-pos))
                  ))
    img))
 
#|
(let ((img
       (extract-image
        "/Volumes/Fornax-T9/Astro/Seestar S50/First Light 240508/3c273/Seestar S50 Subs/Light_Unknown_10.0s_IRCUT_20240508-235415.fit"
        ;; "/Volumes/Fornax-T9/Astro/Seestar S50/240509/Kappa Cygni/Kappa Cygni-sub/Light_Kappa Cygni_10.0s_LP_20240510-014944.fit"
        )))
  (query-header (img-hdr img) "NAXIS")
  (show-img 'img img))
 |#

