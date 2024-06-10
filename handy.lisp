;; handy.lisp
;;
;; DM/RAL  2024/06/10 06:03:31 UTC
;; ----------------------------------

(in-package :com.ral.photometry)

(defun clip (x lo hi)
  (if (< x lo)
      lo
    (if (> x hi)
        hi
      x)))

(defun log10 (x)
  (log x 10))

(defun sqr (x)
  (* x x))

(defconstant +mad/sd+  1.4826f0)

(defun sd-to-mad (sigma) ;; convert to sigma units
  (* sigma +mad/sd+)) ;; assumes underlying is Normal Distr.

(defun mad-to-sd (mads)
  (/ mads +mad/sd+))

;; ---------------------------------------------------------------------

(defun poly-eval (x coffs)
  ;; Horner's rule...
  (reduce (lambda (c acc)
            (+ c (* x acc)))
          coffs
          :from-end t
          :initial-value 0))

(defun erfc (x)
  ;; Approx good to 1.2e-7
  #F
  (let* ((z   (abs (float x 1.0d0)))
         (v   (/ (+ 1.0d0 (* 0.5d0 z))))
         (ans (* v (exp (- (poly-eval v
                                      '(-1.26551223d0
                                        1.00002368d0
                                        0.37409196d0
                                        0.09678418d0
                                        -0.18628806d0
                                        0.27886807d0
                                        -1.13520398d0
                                        1.48851587d0
                                        -0.82215223d0
                                        0.17087277d0))
                           (* z z))
                        ))))
    (if (minusp x)
        (- 2.0d0 ans)
      ans)))

(defun erf (x)
  (- 1.0d0 (erfc x)))

(defun nfc (x)
  ;; NFC(x) represents the integral of the Normal Distribution from x to +inf,
  ;; x measured in stdev units
  (* 0.5d0 (erfc (* x #.(/ (sqrt 2.0d0))) )))

(defun nf (x)
  ;; NF(x) represents the integral of the Normal Distribution from -inf to x
  ;; x measured in stdev units
  (- 1.0d0 (nfc x)))

#|
(log (nfc 5) 10) ;; what is 5 σ ? about 1 in 3.5M
(/ (nfc 5))
(let ((x 5))
  (erfc (* x #.(/ (sqrt 2.0d0))) ))
|#

;; -----------------------------------------------
;; Poisson random deviates

(defun gammln (xx)
  (let* ((y   xx)
         (tmp (+ xx 5.5d0))
         (tmp (- tmp
                 (* (+ xx 0.5d0)
                    (log tmp))))
         (ser  1.000000000190015d0)
         (cof  #( 76.18009172947146d0
                 -86.50532032941677d0
                  24.01409824083091d0
                  -1.231739572450155d0
                   0.1208650973866179d-2
                  -0.5395239384953d-5)
               ))
    (loop for j from 0 to 5 do
            (progn
              (incf y 1d0)
              (incf ser (/ (aref cof j) y))
              ))
    (- (log (/ (* 2.506628274631005d0 ser) xx)) tmp)
    ))

(defun poidev (n)
  (let ((xm  (float (round n) 1d0)))
    (if (< xm 12d0)
        (let* ((lim  (exp (- xm))))
          (um:nlet iter ((prod (lw:mt-random 1d0))
                         (ct   0))
            (if (<= prod lim)
                ct
              (go-iter (* prod (lw:mt-random 1d0)) (1+ ct)))
            ))
      ;; else
      (let* ((sq   (sqrt (* 2.0 xm)))
             (alxm (log xm))
             (g    (- (* xm alxm) (gammln (+ xm 1d0)))))
        (um:nlet iter ()
          (let* ((y   (tan (* pi (lw:mt-random 1d0))))
                 (em  (let* ((em  (+ xm (* sq y))))
                        (if (minusp em)
                            (go-iter)
                          (floor em))
                        ))
                 (tt  (* 0.9d0
                         (+ 1d0 (* y y))
                         (exp (- (* em alxm)
                                 (gammln (+ em 1d0))
                                 g))
                         )))
            (if (> (lw:mt-random 1d0) tt)
                (go-iter)
              em)))
        ))))

#|
(plt:histogram 'plt
               (loop repeat 10000 collect
                       (poidev 11))
           :clear t
           )
(plt:Fplot 'plt '(0 100) #'gammln :clear t)

(let ((med 11))
  (plt:histogram 'plt
                 (loop repeat 10000 collect
                         (loop repeat 49 sum (- (poidev 11) med)))
                 :clear t
                 ))
 |#
         
