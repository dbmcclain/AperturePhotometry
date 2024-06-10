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

(defun copy-array (arr &optional fn)
  (let* ((arr2  (apply #'make-image-array (array-dimensions arr)))
         (vec  (vm:make-overlay-vector arr))
         (vec2 (vm:make-overlay-vector arr2)))
    (if fn
        (map-into vec2 fn vec)
      ;; else
      (replace vec2 vec))
    arr2))

(defun fill-array (arr box val)
  (let* ((lf  (box-left box))
         (wd  (box-width box)))
    (loop for row from (box-top box) below (box-bottom box) do
          (let* ((start (array-row-major-index arr row lf))
                 (end   (+ start wd))
                 (vec   (vm:make-overlay-vector arr :start start :end end)))
            (fill vec val))
            )))

(defun sum-array (arr box)
  (let* ((lf  (box-left box))
         (wd  (box-width box)))
    (loop for row from (box-top box) below (box-bottom box) sum
          (let* ((start  (array-row-major-index arr row lf))
                 (end    (+ start wd))
                 (vec    (vm:make-overlay-vector arr :start start :end end)))
            (vm:total vec))
            )))

(defun array-row (arr row &key (start 0) end)
  ;; Careful! This returns an overlay on the array, not a copy of the
  ;; row.
  (let* ((offs  (array-row-major-index arr row start))
         (end   (+ offs (if end
                            (- end start)
                          (array-dimension arr 1)))))
    (vm:make-overlay-vector arr
                            :start offs
                            :end   end)))

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
  (let+ ((dst-arr (make-image-array (- (box-bottom box) (box-top box))
                                    (- (box-right box)  (box-left box)))))
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
         (vs   (mapcar #'vm:make-overlay-vector (cons arr arrs)))
         (vdst (vm:make-overlay-vector ans)))
    (apply #'map-into vdst fn vs)
    ans))

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

  
