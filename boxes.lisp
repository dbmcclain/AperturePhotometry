;; boxes.lisp
;;
;; DM/RAL  2024/06/09 18:50:35 UTC
;; ----------------------------------

(in-package #:com.ral.photometry)

;; ------------------------------------------------------------------------------
;; Boxes

(declaim (inline box-left box-top box-right box-bottom make-box))

(defun box-left (box)
  (first box))

(defun box-right (box)
  (third box))

(defun box-top (box)
  (second box))

(defun box-bottom (box)
  (fourth box))

(defun box-width (box)
  (- (box-right box)
     (box-left box)))

(defun box-height (box)
  (- (box-bottom box)
     (box-top box)))

(defun box-area (box)
  (* (- (box-right box)  (box-left box))
     (- (box-bottom box) (box-top box))))

(defun make-box (&rest ltrb)
  ltrb)

(defun make-box-ltwh (lf tp wd ht)
  `(,lf ,tp ,(+ lf wd) ,(+ tp ht)))

(defun make-box-centered (xctr yctr wd ht)
  (make-box-ltwh (- xctr (truncate wd 2))
                 (- yctr (truncate ht 2))
                 wd ht))

(defun make-box-of-radius (xc yc radius)
  (make-box (- xc radius) (- yc radius)
            (+ xc radius 1) (+ yc radius 1)))

(defun make-box-of-array (arr)
  (let+ (( (ht wd) (array-dimensions arr)))
    (make-box 0 0 wd ht)))

(defun make-box-of-img (img)
  (make-box-of-array (img-arr img)))

(defun box-intersection (box1 box2)
  (make-box (max (box-left box1)   (box-left box2))
            (max (box-top box1)    (box-top box2))
            (min (box-right box1)  (box-right box2))
            (min (box-bottom box1) (box-bottom box2))))

(defun box-union (box1 box2)
  (make-box (min (box-left box1)   (box-left box2))
            (min (box-top box1)    (box-top box2))
            (max (box-right box1)  (box-right box2))
            (max (box-bottom box1) (box-bottom box2))))

(defun box-center (box)
  (values (truncate (+ (box-left box) (box-right box)) 2)
          (truncate (+ (box-top box)  (box-bottom box)) 2)))

(defun move-box (box dx dy)
  (let+ (( (lf tp rt bt) box))
    `(,(+ lf dx) ,(+ tp dy) ,(+ rt dx) ,(+ bt dy))
    ))

(defun inset-box (box dx dy)
  (let+ (( (lf tp rt bt) box))
    `(,(+ lf dx) ,(+ tp dy) ,(- rt dx) ,(- bt dy))
    ))

(defun box-contains-pt-p (box x y)
  (let+ (( (lf tp rt bt) box))
    (and (<= lf x (1- rt))
         (<= tp y (1- bt)))
    ))

(defun box-contains-box-p (box-outer box-inner)
  (let+ (( (lf_o tp_o rt_o bt_o) box-outer)
         ( (lf_i tp_i rt_i bt_i) box-inner))
    (and (<= lf_o lf_i)
         (<= tp_o tp_i)
k         (>= rt_o rt_i)
         (>= bt_o bt_i))
    ))

(defun canonicalize-box (box)
  (let+ (( (lf tp rt bt) box))
    (make-box (min lf rt) (min tp bt)
              (max lf rt) (max tp bt))))

(defun copy-box (box)
  (copy-seq box))
