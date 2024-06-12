
(in-package :com.ral.photometry)

#|
(defun gen-fake-star (img &optional (magsel #'identity))
  (let ((imgstk  (make-image-array 19 19
                                   :initial-element 0.0f0))
        (arr     (img-arr img))
        (sum-snr 0)
        (ct      0))
      (dolist (star (img-stars img))
        (with-accessors ((xc  star-x)
                         (yc  star-y)
                         (snr star-snr)) star
          (when (funcall magsel snr)
            (incf ct)
            (incf sum-snr snr)
            (let* ((med  (ring-med arr yc xc))
                   (box  (make-box-of-radius xc yc 9))
                   (wd   (box-width box))
                   (lf   (box-left box)))
              (loop for src-row from (box-top box) below (box-bottom box)
                    for dst-row from 0
                    do
                      (let* ((start (array-row-major-index arr src-row lf))
                             (end   (+ start wd))
                             (vsrc  (vm:make-overlay-vector arr :start start :end end))
                             (vdst  (array-row imgstk dst-row)))
                        (map-into vdst (lambda (xsrc xdst)
                                         (+ xdst (* snr (- xsrc med))))
                                  vsrc vdst)
                        ))
              ))))
      ;; zap < 1% background
      (let* ((pk  (aref imgstk 9 9))
             (lo  (* 0.01 pk)))
        (loop for ix from 0 below (array-total-size imgstk) do
                (if (< (row-major-aref imgstk ix) lo)
                    (setf (row-major-aref imgstk ix) 0f0)
                  (setf (row-major-aref imgstk ix) (/ (row-major-aref imgstk ix) pk)))))
      (plt:tvscl 'stack imgstk
                 :clear t
                 :magn  16
                 :flipv t)
      (let* ((box  (make-box-of-radius 9 9 5))
             (flux (sum-array imgstk box))
             (mag  (magn img flux)))
        (print (list :stack-mag mag :count ct))
        (setf *fake-star* `(,imgstk 9 ,mag)))
      ))
|#
#|
(gen-fake-star *sub* (lambda (snr) (< 50 snr 150)))
(show-fake-star *fake-star*)
 |#

;; -------------------------------------------------------------------
;; Make Gaussian Fake Star

(defun gaussian (x sigma)
  (let* ((z  (/ x sigma)))
    (exp (* -0.5f0 (sqr z)))
    ))

(defun make-gaussian-fake-star (&key (sigma 0.75f0) (radius *core-radius*))
  (let+ ((xtnt      (1+ (* 2 radius)))
         (xs        (vops:voffset (- radius) (vm:framp xtnt)))
         (exps      (map 'vector (um:rcurry #'gaussian (coerce sigma 'single-float)) xs))
         (arr       (vm:outer-prod exps exps))
         (box       (make-box-of-array arr))
         (ksum      (vm:total arr))
         (k2sum     (vm:total (map-array #'* arr arr)))
         (npix      (* xtnt xtnt))
         (Δ         (- (* npix k2sum)
                       (* ksum ksum))))
    ;; Engineering mag of fake star will be 0.0.
    (make-fake
     :krnl   arr
     :Δ      Δ
     :npix   npix
     :ksum   ksum
     :k2sum  k2sum
     :radius radius
     :sigma  sigma
     :box    box
     )))

(defun show-fake-star (fake-star)
  (let+ ((arr    (fake-krnl fake-star))
         (radius (fake-radius fake-star))
         (dists  (vops:voffset (- radius) (vm:framp (1+ (* 2 radius))))))
    (plt:tvscl 'fake-star arr
               :magn  16
               :zlog  t
               :clear t)
    (plt:spline 'fake-slice dists (array-row arr radius)
                :clear t
                :title "Central Slice of Fake Star"
                :xtitle "Dist from Center [pix]"
                :ytitle "Amplitude"
                :symbol :circle
                :legend "Horizontal")
    (plt:spline 'fake-slice dists (array-col arr radius)
                :color :red
                :symbol :circle
                :legend "Vertical")
    ))

#|
(let ((*core-radius* 3))
  (show-fake-star (make-gaussian-fake-star :sigma 0.75)))
(let ((*core-radius* 7))
  (show-fake-star (make-gaussian-fake-star :sigma 1.3)))
(show-fake-star (img-fake-star *saved-img*))
 |#

