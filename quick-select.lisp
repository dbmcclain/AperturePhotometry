;; --------------------------------------

(defun pick-pivot (lst)
  (when lst
    (if (< (length lst) 5)
        (vm:median lst)
      (let* ((chunks (um:group lst 5))
             (full-chunks (remove-if (lambda (chunk)
                                       (< (length chunk) 5))
                                     chunks))
             (sorted-groups (mapcar (um:rcurry #'sort #'<) full-chunks))
             (medians       (mapcar #'third sorted-groups)))
        (quickselect-median medians 'pick-pivot)
        ))))
             
(defun quick-select (lst k &optional (pivot-fn 'random-choice))
  ;; select k'th element in lst
  #F
  (let* ((len  (length lst)))
    (if (= 1 len)
        (progn
          (assert (= k 0))
          0)
      (let* ((pivot (funcall pivot-fn lst))
             lows
             highs
             pivots)
        (map nil (lambda (x)
                   (cond ((< x pivot)
                          (push x lows))
                         ((> x pivot)
                          (push x highs))
                         (t
                          (push x pivots))
                         ))
             lst)
        (let* ((llen (length lows)))
          (if (< k llen)
              (quick-select lows k pivot-fn)
            (let* ((plen (length pivots)))
              (if (< k (+ llen plen))
                  (elt pivots 0)
                (quick-select highs (- k llen plen) pivot-fn)
                ))))
        ))))

(defun random-choice (lst)
  (let ((nel (length lst)))
    (elt lst (random nel))))

(defun quick-select-median (lst &optional (pivot-fn 'random-choice))
  (let* ((len  (length lst))
         (len/2 (truncate len 2)))
    (if (oddp len)
        (quick-select lst len/2 pivot-fn)
      (* 1/2 (+ (quick-select lst (1- len/2) pivot-fn)
                (quick-select lst len/2 pivot-fn)))
      )))

#|
(quick-select-median #(1 2 3 4 5 6 7 8 9 10))
(vm:median #(1 2 3 4 5 6 7 8 9 10))

(time (quick-select-median (vm:make-overlay-vector (img-arr *raw-img*))))
(time (vm:median (img-arr *raw-img*)))
|#

