
(in-package #:com.ral.photometry)

(defun select-region (stars iy ix dist)
  ;; Assumes star lists are sorted in ascending Y order
  (when stars
    (let+ (( (xlo ylo xhi yhi) (make-box-of-radius ix iy dist))
           (start (position-if (lambda (star)
                                 (>= (star-y star) ylo))
                               stars))
           (end   (and start
                       (position-if (lambda (star)
                                      (>= (star-y star) yhi))
                                    stars
                                    :start start)))
           (ysel  (and start
                       (subseq stars start end))))
      (remove-if (lambda (star)
                   (let ((x (star-x star)))
                     (or (< x xlo)
                         (>= x xhi))))
                 ysel)
      )))
         
(defun generate-fake-star-positions (img viz)
  (let* ((arr      (img-arr img))
         (wd       (array-dimension arr 1))
         (ht       (array-dimension arr 0))
         (min-dist (* 2 *ring-radius*))
         (margin   *ring-radius*)
         (wd-range (- wd 1 (* 2 margin)))
         (ht-range (- ht 1 (* 2 margin)))
         (stars    (sort (img-stars img)
                         #'<
                         :key #'star-y)))
    (um:nlet iter ((ct     0)
                   (trials 0)
                   (fakes  nil))
      (if (>= ct 100)
          (progn
            (format t "~%100 fake stars after ~d trials" trials)
            fakes)
        ;; else
        (let* ((ix  (+ margin (random wd-range)))
               (iy  (+ margin (random ht-range)))
               (selected-stars (select-region stars iy ix min-dist))
               (selected-fakes (select-region
                                (sort fakes
                                      #'<
                                      :key #'star-y)
                                iy ix min-dist)
                               ))
          (cond ((or selected-stars selected-fakes)
                 (princ #\-)
                 (go-iter ct (1+ trials) fakes))
                (t
                 (princ #\.)
                 (when viz
                   (plt:draw-rect 'img ix iy 9 9
                                  :raw t
                                  :border-color :magenta
                                  :filled nil
                                  :border-thick 1))
                 (go-iter (1+ ct) (1+ trials) (cons (make-star
                                                     :x  ix
                                                     :y  iy)
                                                    fakes)))
                ))))))

(defun plant-fake-stars (img fake-star-positions fake-star mag viz)
  ;; Plant the fake stars and see how well we can recover them.
  (let+ (( (fake-star fake-radius fake-mag) fake-star)
         (fake-img (copy-img img))
         (fake-arr (setf (img-arr fake-img)
                         (copy-img-array (img-arr img))))
         (sf       (expt 10 (* -0.4f0 (- mag fake-mag)))))
    (print (list :sf sf))
    (dolist (star fake-star-positions)
      (with-accessors ((ix star-x)
                       (iy star-y)) star
        (loop for i from (- iy fake-radius) to (+ iy fake-radius)
              for ii from 0
              do
                (loop for j from (- ix fake-radius) to (+ ix fake-radius)
                      for jj from 0
                      do
                        (incf (aref fake-arr i j) (float (round (poidev (* sf (aref fake-star ii jj)))) 1f0))
                      ))))
    (when viz
      (show-img 'img+fakes fake-img))
    fake-img
    ))

(defun recover-fake-stars (img fake-star-positions &optional (thresh 5))
  ;; Show recovery results
  (let* ((found-stars (sort  ;; list of found stars sorted in y index
                       (find-stars img thresh)
                       #'<
                       :key 'star-y))
         (resids      (copy-seq fake-star-positions))
         (mags        nil)
         (snrs        nil))
    #|
    (plt:with-delayed-update ('fake-sky)
      (show-img 'fake-sky img)
      (hilight-stars 'fake-sky found-stars  :green)
      (hilight-stars 'fake-sky fake-star-positions :magenta))
    |#
    (dolist (star fake-star-positions)
      (with-accessors ((ix star-x)
                       (iy star-y)) star
        (let* ((found (find-if (lambda (star)
                                 (with-accessors ((sx  star-x)
                                                  (sy  star-y)) star
                                   (< (abs (complex (- ix sx) (- iy sy))) 2)))
                               found-stars)))
          (when found
            (with-accessors ((mag star-mag)
                             (snr star-snr)) found
              #|
              (print (list
                      :found :Δ (abs (complex (- ix (first found))
                                              (- iy (second found))))
                      :mag (third found)))
              |#
              (push mag mags)
              (push snr snrs)
              (setf resids (remove star resids))))
          )))
    
    (if mags
        (let* ((mn    (vm:mean  mags))
               (sd    (vm:stdev mags mn))
               (mnsnr (vm:mean  snrs)))
          #||#
          (plt:plot 'recovered mags
                    :clear t
                    :title "Recovered Stars"
                    :ytitle "Meas Mag"
                    :xtitle "Star #"
                    :yrange '(16 8)
                    :symbol :circle
                    )
          (plt:plot 'recovered `(0 ,(length mags)) `(,mn ,mn)
                    :color :gray50
                    :thick 2)
          (plt:plot 'recovered `(0 ,(length mags)) `(,(+ mn sd) ,(+ mn sd))
                    :color :gray60
                    :thick 2)
          (plt:plot 'recovered `(0 ,(length mags)) `(,(- mn sd) ,(- mn sd))
                    :color :gray60
                    :thick 2)
          #||#
          (list :recovered-mags :nbr (length mags) :mn mn :sd sd :snr mnsnr))
      ;; else
      "No stars recovered")
    ;; resids
    ;; (values)
    ))

