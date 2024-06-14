;; array-ops.lisp
;;
;; DM/RAL  2024/06/10 06:02:38 UTC
;; ----------------------------------

(in-package :com.ral.photometry)

;; -------------------------------------------------------

(defun make-image-array (nrows ncols &rest args)
  (apply #'make-array `(,nrows ,ncols)
         :element-type 'single-float
         args))

(defun vec (obj &rest args)
  (apply #'vm:make-overlay-vector obj args))

(defun copy-array-into (dst arr &optional fn)
  (let ((vdst  (vec dst))
        (vsrc  (vec arr)))
    (if fn
        (map-into vdst fn vsrc)
      ;; else
      (replace vdst vsrc))
    dst))

(defun copy-array (arr &optional fn)
  (copy-array-into (make-array (array-dimensions arr)
                               :element-type (array-element-type arr))
                   arr fn))

(defun copy-img-array (arr &optional fn)
  (copy-array-into (apply #'make-image-array (array-dimensions arr))
                   arr fn))

(defun array-fill-in-box (arr box val)
  (let* ((lf  (box-left box))
         (wd  (box-width box)))
    (loop for row from (box-top box) below (box-bottom box) do
          (let* ((start (array-row-major-index arr row lf))
                 (end   (+ start wd))
                 (vec   (vec arr :start start :end end)))
            (fill vec val))
            )))

(defun array-fill (arr val)
  (fill (vec arr) val)
  arr)

(defun reduce-array (fn arr &rest args)
  (apply #'reduce fn (vec arr) args))

(defun reduce-array-in-box (fn arr box &rest args)
  (apply #'reduce-array fn (extract-subarray arr box) args))

(defun pos-in-array (arr fn &rest args)
  (truncate (apply fn (vec arr) args) (array-dimension arr 1)))

(defun array-pos (x arr &rest args)
  (apply #'pos-in-array arr (um:curry #'position x) args))

(defun array-pos-if (fn arr &rest args)
  (apply #'pos-in-array arr (um:curry #'position-if fn) args))

(defun pos-in-box (arr box fn &rest args)
  (let+ ((:mvb (row col) (apply #'pos-in-array (extract-subarray arr box) fn args)))
    (values (+ row (box-top box))
            (+ col (box-left box)))
    ))

(defun array-pos-in-box (x arr box &rest args)
  (apply #'pos-in-box arr box (um:curry #'array-pos x) args))

(defun array-pos-if-in-box (fn arr box &rest args)
  (apply #'pos-in-box arr box (um:curry #'array-pos-if fn) args))

(defun max-pos (vec)
  ;; position of max value in vector
  (position (reduce #'max vec) vec))

(defun array-max-pos (arr)
  ;; 2D index location of max value in array
  (pos-in-array arr #'max-pos))

(defun array-max-pos-in-box (arr box)
  (pos-in-box arr box #'max-pos))


(defun sum-array-in-box (arr box)
  (let* ((lf  (box-left box))
         (wd  (box-width box)))
    (loop for row from (box-top box) below (box-bottom box) sum
          (let* ((start  (array-row-major-index arr row lf))
                 (end    (+ start wd))
                 (vec    (vec arr :start start :end end)))
            (vm:total vec))
            )))

(defun array-row (arr row &key (start 0) end)
  ;; Careful! This returns an overlay on the array, not a copy of the
  ;; row.
  (let* ((offs  (array-row-major-index arr row start))
         (end   (+ offs (if end
                            (- end start)
                          (array-dimension arr 1)))))
    (vec arr :start offs :end   end)
    ))

(defun array-col (arr col &key (start 0) end)
  (let* ((ht   (- (or end (array-dimension arr 0)) start))
         (vec  (make-array ht
                           :element-type (array-element-type arr))))
    (loop for src-row from start below (+ start ht)
          for dst-row from 0
          do
            (setf (aref vec dst-row) (aref arr src-row col)))
    vec))

(defun extract-subarray (arr box)
  (let+ ((dst-arr (make-array `(,(- (box-bottom box) (box-top box))
                                ,(- (box-right box)  (box-left box)))
                              :element-type (array-element-type arr))))
    (loop for src-row from (box-top box) below (box-bottom box)
          for dst-row from 0
          do
            (let ((vsrc (array-row arr src-row))
                  (vdst (array-row dst-arr dst-row)))
              (replace vdst vsrc :start2 (box-left box))
              ))
    dst-arr))

(defun implant-subarray (arr-dst arr-src row-org col-org)
  (let+ ((box-dst  (move-box (make-box-of-array arr-src)
                             col-org row-org)))
    (loop for src-row from 0
          for dst-row from (box-top box-dst) below (box-bottom box-dst)
          do
            (let ((vsrc (array-row arr-src src-row))
                  (vdst (array-row arr-dst dst-row)))
              (replace vdst vsrc :start1 (box-left box-dst) :end1 (box-right box-dst))
              ))
    arr-dst))

(defun make-similar-array (arr &rest args)
  (apply #'make-array (array-dimensions arr)
         :element-type (array-element-type arr)
         args))

(defun map-array (fn arr &rest arrs)
  (let* ((ans  (make-similar-array arr))
         (vs   (mapcar #'vec (cons arr arrs)))
         (vdst (vec ans)))
    (apply #'map-into vdst fn vs)
    ans))

(defun map-array-into (dst fn arr &rest arrs)
  (let* ((vs   (mapcar #'vec (cons arr arrs)))
         (vdst (vec dst)))
    (apply #'map-into vdst fn vs)
    dst))

(defun array-unop (arr fn)
  (map-array fn arr))

(defun array-binop (arr-a arr-b fn)
  (map-array fn arr-a arr-b))

(defun array+a (arr-a arr-b)
  (array-binop arr-a arr-b #'+))

(defun array-a (arr-a arr-b)
  (array-binop arr-a arr-b #'-))

(defun array*a (arr-a arr-b)
  (array-binop arr-a arr-b #'*))

(defun array/a (arr-a arr-b)
  (array-binop arr-a arr-b #'/))

(defun array+k (arr kval)
  (array-unop arr (um:curry #'+ kval)))

(defun array-k (arr kval)
  (array-unop arr (um:rcurry #'- kval)))

(defun k-array (arr kval)
  (array-unop arr (um:curry #'- kval)))

(defun array/k (arr kval)
  (array-unop arr (um:rcurry #'/ kval)))

(defun k/array (arr kval)
  (array-unop arr (um:curry #'/ kval)))

