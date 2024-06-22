;; Compare PixInsight Aperture Photometry against SIMBA Catalogs
;; DM/RAL  2024/06/17 12:57:40 UTC
;; ---------------------------------------------------------------------

(in-package :com.ral.photometry)

;; For the Gaia DR2 down to mag 15:
;; Run PixInsight AperturePhotometry script against a G-plane image with Aperture = 3 pix radius and save CSV files.
(defun show-pi-photom ()
  (labels ((remove-nth (vec pos)
             (if (null pos)
                 vec
               (concatenate 'vector (subseq vec 0 pos) (subseq vec (1+ pos)))))
           
           (filter-bad (vec &rest vecs)
             (um:nlet iter ((ix 0)
                            (vec vec)
                            (vecs vecs))
               (if (null ix)
                   (values-list (cons vec vecs))
                 (let ((pos (position-if (um:rcurry #'<= 0) vec)))
                   (go-iter pos
                            (remove-nth vec pos)
                            (mapcar (um:rcurry #'remove-nth pos) vecs))
                   )))))
    
  (um:with-remembered-filename (fname "Choose a file" :photom nil
                                    :filter "*.csv")
  (let+ ((data   (csv:read-file :fname fname :skip-lines 4 :hdr-lines 1 :type :semiv))
         (grp    (csv:get-group nil data))
         (gmag   (map 'vector #'read-from-string (csv:get-column "Catalog_Gmag" grp)))
         (flux3  (map 'vector #'read-from-string (csv:get-column "Flux_Ap3_1" grp)))
         (gain   80)
         (npix   1) ;; (* pi 3 3)))
         (:mvb (flux3 gmag) (filter-bad flux3 gmag))
         (lflux  (map 'vector (um:rcurry #'log 10) flux3)))
    (multiple-value-bind (y0 sigma)
        (linfit:regress-fixed-slope lflux gmag
                                    1     ;; weights
                                    -2.5) ;; slope
      (plt:plot 'flux flux3 gmag
                :clear t
                :xlog  t
                :yrange '(16 5)
                :symbol :dot ;; cross
                :alpha 0.5
                :title "Cat GMag vs Flux_Ap3"
                :xtitle "Flux [ADU]"
                :ytitle "GMag [mag]"
                )
      (print (list y0 sigma))
      (flet ((fit (flux y0)
               (let ((lf (log flux 10)))
                 (+ y0 (* -2.5 lf)))))
        (plt:fplot 'flux '(0.1e3 100e6) (um:rcurry #'fit y0)
                 :color :red
                 :legend "Best Fit Offs")
        (plt:fplot 'flux '(0.1e3 100e6) (um:rcurry #'fit (+ y0 sigma))
                   :color :blue
                   :legend (format nil "Offs±~4,2F mag" sigma))
        (plt:fplot 'flux '(0.1e3 100e6) (um:rcurry #'fit (- y0 sigma))
                   :color :blue))
      (plt:draw-text 'flux
                     (format nil "Mag Offs = ~5,2F" y0)
                     '((:data 11e3) (:data 6.3)))
      (plt:draw-text 'flux
                     (format nil "Sigma = ~5,2F" sigma)
                     '((:data 11e3) (:data 6.8)))
      )))))

#|
(show-pi-photom)

;; This is the string that PixInsight claims it used...
https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/345/gaia2&-c=240.005064 25.967865&-c.r=1.088996&-c.u=deg&-out.form=|&-out.add=RA_ICRS,DE_ICRS&-out=pmRA&-out=pmDE&-out=Source&-out=Gmag&-out=RPmag&-out=BPmag&-out=Plx&-out=RV&-out=Rad&-out=Lum&Gmag=0..15

;; This, from Martin Simmons, works in the Lisp Webkit Demo
https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source%3DI%2F345%2Fgaia2&-c%3D240.005064%2025.967865&-c.r%3D1.088996&-c.u%3Ddeg&-out.form%3D%7C&-out.add%3DRA_ICRS%2CDE_ICRS&-out%3DpmRA&-out%3DpmDE&-out%3DSource&-out%3DGmag&-out%3DRPmag&-out%3DBPmag&-out%3DPlx&-out%3DRV&-out%3DRad&-out%3DLum&Gmag%3D0..15

;; This works from a command line, and from Lisp with SYS:OPEN-PIPE
;; Just need to protect the #\|, a shell symbol for pipes, with single-quotes
;; And the #\space between RA and Dec needs to be %20
vizquery -mime=csv -source=I/345/gaia2 -c=240.005064%2025.967865 -c.r=1.088996 -c.u=deg -out.form='|' -out.add=RA_ICRS,DE_ICRS -out=pmRA -out=pmDE -out=Source -out=Gmag -out=RPmag -out=BPmag -out=Plx -out=RV -out=Rad -out=Lum Gmag=0..15

;; see if Lisp Webkit Demo really needs those #\= turned into %3D? No, and #\/ to %2D not needed either.
;; Needed: the %7c for the #\|
;;         the %20 for the #\space between RA & Dec
https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/345/gaia2&-c=240.005064%2025.967865&-c.r=1.088996&-c.u=deg&-out.form=%7c&-out.add=RA_ICRS,CDE_ICRS&-out=pmRA&-out=pmDE&-out=Source&-out=Gmag&-out=RPmag&-out=BPmag&-out=Plx&-out=RV&-out=Rad&-out=Lum&Gmag=0..15
|#
#|
(defun tst-viz ()
  (let ((s1
"vizquery -mime=csv -source=I/345/gaia2 -c=240.005064%2025.967865 -c.r=1.088996 -c.u=deg -out.form='|' -out.add=RA_ICRS,DE_ICRS -out=pmRA -out=pmDE -out=Source -out=Gmag -out=RPmag -out=BPmag -out=Plx -out=RV -out=Rad -out=Lum Gmag=0..15"))
    (with-open-stream (sin (sys:open-pipe s1 :direction :input))
      (loop for line = (read-line sin nil sin nil)
            until (eql line sin)
            collect line))
    ))
  
(defvar *data*)

(setf *data* (tst-viz))
(sys:call-system-showing-output "echo $SHELL; echo $PATH")
(nth 11 *data*)
|#
;; --------------------------------------------------------------------

(defun dtor (x)
  (/ (* pi x) 180))

#|
(defun rtod (x)
  (/ (* 180 x) pi))

(defun rtop (x y)
  (values (abs (complex x y))
          (rtod (atan y x))))

(defun ptor (rho phi)
  (let ((z  (* rho (cis (dtor phi)))))
    (values (realpart z)
            (imagpart z))
    ))
|#

(defun rtopd (x y)
  (values (abs (complex x y))
          (rtod (atan y x))))

(defun ptord (r th)
  (let ((cs  (cis (dtor th))))
    (values (* r (realpart cs))
            (* r (imagpart cs)))
    ))

(defun vcross (v1 v2)
  (let+ (((x1 y1 z1) v1)
         ((x2 y2 z2) v2))
    `(,(- (* y1 z2) (* z1 y2))
      ,(- (* z1 x2) (* x1 z2))
      ,(- (* x1 y2) (* x2 y1)))
    ))

(defun v*k (k v)
  (mapcar (um:curry #'* k) v))

(defun v+v (v1 v2)
  (mapcar #'+ v1 v2))

(defun vdot (v1 v2)
  (let ((sum 0))
    (map nil (lambda (x y)
               (incf sum (* x y)))
         v1 v2)
    sum
    ))

(defun radec2xyz (ra dec)
  (let+ ((:mvb (rho z) (ptord 1 dec))
         (:mvb (x y)   (ptord rho ra)))
    `(,x ,y ,z)
    ))

(defun xyz2radec (v)
  (let+ (( (x y z) v)
         (:mvb (rho thet) (rtopd x y))
         (:mvb (_ phi)   (rtopd rho z)))
    (values (mod thet 360)
            phi)))

(defun parallactic-angle (img)
  (let+ ((hdr    (img-hdr img))
         ;; transform from pixels to IWS
         (cd1_1  (nquery-header hdr "CD1_1"))
         (cd1_2  (nquery-header hdr "CD1_2"))
         (cd2_1  (nquery-header hdr "CD2_1"))
         (cd2_2  (nquery-header hdr "CD2_2"))
         (:mvb (_ th)  (rtopd cd1_1 cd1_2)))
    th
    ))

#|
(parallactic-angle *saved-img*)
(canon-view *saved-img*)
(show-img 'img *saved-img*)
 |#

;; ---------------------------------------------------

(defstruct rotxform
  theta   ;; parallactic angle + 180 deg
  xf1_1   ;; Cos theta
  xf1_2   ;; Sin theta
  xc      ;; (xc, yc) center of untilted image
  yc
  cxc     ;; (cxc, cyx) center of tilted (canonical) image
  cyc
  arr)    ;; original img array of untilted img

(defun rotxform (xform x y)
  ;; Transform from untilted image to tilted
  (with-accessors ((xf1_1   rotxform-xf1_1)
                   (xf1_2   rotxform-xf1_2)
                   (xc      rotxform-xc)
                   (yc      rotxform-yc)
                   (cxc     rotxform-cxc)
                   (cyc     rotxform-cyc)) xform
    (let+ ((xmc  (- x xc))
           (ymc  (- y yc))
           (cxmc (+ (* xmc xf1_1) (* ymc xf1_2)))
           (cymc (- (* ymc xf1_1) (* xmc xf1_2))))
      (values
       (+ cxmc cxc)
       (+ cymc cyc))
      )))

(defun inv-rotxform (xform x y)
  ;; Transform from tilted image to untilted
  (with-accessors ((xf1_1   rotxform-xf1_1)
                   (xf1_2   rotxform-xf1_2)
                   (xc      rotxform-xc)
                   (yc      rotxform-yc)
                   (cxc     rotxform-cxc)
                   (cyc     rotxform-cyc)) xform
    (let+ ((cxmc  (- x cxc))
           (cymc  (- y cyc))
           (xmc   (- (* cxmc xf1_1) (* cymc xf1_2)))
           (ymc   (+ (* cymc xf1_1) (* cxmc xf1_2))))
      (values
       (+ xmc xc)
       (+ ymc yc))
      )))

;; ---------------------------------------------------

(defun canon-view (img)
  (let+ ((th    (+ 180 (parallactic-angle img)))
         (cisth (cis (dtor th)))
         (xf1_1 (realpart cisth))
         (xf1_2 (imagpart cisth)))
    (labels (#||#
             (proj (x y)
               (values
                (+ (* x xf1_1) (* y xf1_2))
                (- (* y xf1_1) (* x xf1_2))))
             (inv-proj (x y)
               (values
                (- (* x xf1_1) (* y xf1_2))
                (+ (* x xf1_2) (* y xf1_1)))))
      (let+ ((arr (img-arr img))
             ((ht wd) (array-dimensions arr))
             (:mvb (xtl ytl) (proj 0  0))
             (:mvb (xtr ytr) (proj wd 0))
             (:mvb (xbr ybr) (proj wd ht))
             (:mvb (xbl ybl) (proj 0  ht))
             (lf   (min xtl xbl xtr xbr))
             (rt   (max xtr xbr xtl xbl))
             (tp   (min ytl ytr ybl ybr))
             (bt   (max ybl ybr ytl ytr))
             (xc   (/ wd 2))
             (yc   (/ ht 2))
             (cwd  (ceiling (- rt lf)))
             (cht  (ceiling (- bt tp)))
             (cxc  (/ cwd 2))
             (cyc  (/ cht 2))
             (carr (make-image-array cht cwd :initial-element 0f0))
             (cimg (copy-img img)))
        (setf (img-arr cimg) carr
              (img-canon cimg) (make-rotxform
                                :arr   arr
                                :theta th
                                :xf1_1 xf1_1
                                :xf1_2 xf1_2
                                :xc    xc
                                :yc    yc
                                :cxc   cxc
                                :cyc   cyc))
        (loop for row from 0 below cht do
                (loop for col from 0 below cwd do
                        (let+ ((:mvb (sx sy) (inv-proj (- col cxc) (- row cyc)))
                               (srow (round (+ sy yc)))
                               (scol (round (+ sx xc))))
                          (when (array-in-bounds-p arr srow scol)
                            (setf (aref carr row col) (aref arr srow scol)))
                          )))
        (show-img 'canon cimg)
        ))))
         
(defun to-radec (img xp yp)
  (let+ ((hdr    (img-hdr img))
         ;; Center coords
         (cr_ra  (nquery-header hdr "CRVAL1"))
         (cr_dec (nquery-header hdr "CRVAL2"))
         (cr_x   (nquery-header hdr "CRPIX1"))
         (cr_y   (nquery-header hdr "CRPIX2"))

         ;; transform from pixels to IWS
         (cd1_1  (nquery-header hdr "CD1_1"))
         (cd1_2  (nquery-header hdr "CD1_2"))
         (cd2_1  (nquery-header hdr "CD2_1"))
         (cd2_2  (nquery-header hdr "CD2_2"))

         ;; X projections to go from IWS to RA, Dec
         (a_0_0  (nquery-header hdr "A_0_0"))
         (a_0_1  (nquery-header hdr "A_0_1"))
         (a_0_2  (nquery-header hdr "A_0_2"))
         (a_1_0  (nquery-header hdr "A_1_0"))
         (a_1_1  (nquery-header hdr "A_1_1"))
         (a_2_0  (nquery-header hdr "A_2_0"))

         ;; Y projections
         (b_0_0  (nquery-header hdr "B_0_0"))
         (b_0_1  (nquery-header hdr "B_0_1"))
         (b_0_2  (nquery-header hdr "B_0_2"))
         (b_1_0  (nquery-header hdr "B_1_0"))
         (b_1_1  (nquery-header hdr "B_1_1"))
         (b_2_0  (nquery-header hdr "B_2_0")))
    (labels ((to-undistorted (xp yp)
               ;; distortion mapping
               ;; Image coords have origin at top-left corner, middle of pixel 0 is (0.5, 0.5)
               ;; FITS coords have origin at bottom-left corner, pixel 0 is (1, 1)
               ;; Pixel relative (to image center) coords have warp removed:
               ;;   uw = u + F(u,v), F(u,v) = Sum(a_i_j*u^i*v^j)
               ;;   vw = v + G(u,v), G(u,v) = Sum(b_i_j*u^i*v^j)
               ;;
               (let+ ((u   (- (1+ xp) cr_x))
                      (v   (- (1+ yp) cr_y)))
                 (values
                  (+ u
                     a_0_0
                     (* a_0_1 v)
                     (* a_0_2 v v)
                     (* a_1_0 u)
                     (* a_1_1 u v)
                     (* a_2_0 u u))
                  (+ v
                     b_0_0
                     (* b_0_1 v)
                     (* b_0_2 v v)
                     (* b_1_0 u)
                     (* b_1_1 u v)
                     (* b_2_0 u u))
                  )))
             (to-IWC (xw yw)
               ;; Take out skew and rotation of tangent plane
               ;; to become WCS, in degrees
               (values
                (+ (* cd1_1 xw) (* cd1_2 yw))
                (+ (* cd2_1 xw) (* cd2_2 yw))
                ))
             (norm3 (v)
               (let+ (((x y z) v)
                      (norm (sqrt (vdot v v))))
                 `(,(/ x norm)
                   ,(/ y norm)
                   ,(/ z norm))
                 ))
             (iwc-to-xyz (dx dy)
               (let+ ((x (- (dtor dx)))
                      (y (dtor dy))
                      ;; vector r is plate center
                      ((rx ry rz) (radec2xyz cr_ra cr_dec))
                      ;; vector i is r ✕ NCP
                      (norm (abs (complex rx ry)))
                      (ix  (/ ry norm))
                      (iy  (/ (- rx) norm))
                      ;; vector j is i ✕ r
                      (jx  (* iy rz))
                      (jy  (- (* ix rz)))
                      (jz  (- (* ix ry) (* iy rx))))
                 (norm3 `(,(+ (* ix x) (* jx y) rx)
                          ,(+ (* iy x) (* jy y) ry)
                          ,(+ (* jz y) rz)))
                 ))
             (to-thphi (pos)
               (let+ (( (x y z) pos)
                      (:mvb (rho th) (rtopd x y))
                      (:mvb (un phi) (rtopd z rho)))
                 (values (mod th 360) phi un)
                 )))
      (let+ ((:mvb (uw vw)   (to-undistorted xp yp))
             (:mvb (dxw dyw) (to-IWC uw vw))
             (v    (iwc-to-xyz dxw dyw)))
        (xyz2radec v)
        ))))

#|
(to-radec *saved-img* 0 0)              ;; top-left
(to-radec *saved-img* 1080 0)           ;; top-right
(to-radec *saved-img* 1080 1920)        ;; bot-right      
(to-radec *saved-img* 0 1920)           ;; bot-left
  |#
  
(defun get-star-positions (img)
  (dolist (star (img-stars img))
    (let+ ((:mvb (α δ) (to-radec img (star-x star) (star-y star))))
      (setf (star-ra star)  (mod α 360)
            (star-dec star) δ)
      )))

(defun get-catalog (img)
  (let+ ((hdr    (img-hdr img))
         (cr_ra  (nquery-header hdr "CRVAL1"))
         (cr_dec (nquery-header hdr "CRVAL2"))
         ( (nypix nxpix)  (array-dimensions (img-arr img)))
         (xpxsiz (nquery-header hdr "XPIXSZ"))
         (ypxsiz (nquery-header hdr "YPIXSZ"))
         (foclen (nquery-header hdr "FOCALLEN"))
         (radius (* 1.1
                    (rtod
                     (atan (/ (max (* nxpix xpxsiz) (* nypix ypxsiz))
                              foclen 2000)))
                    ))
         (cmdstr (format nil
"vizquery -mime=csv -source=I/345/gaia2 -c=~F%20~F -c.r=~F -c.u=deg -out.form='|' -out.add=RA_ICRS,DE_ICRS -out=pmRA -out=pmDE -out=Source -out=Gmag -out=RPmag -out=BPmag -out=Plx -out=RV -out=Rad -out=Lum Gmag=0..15"
cr_ra cr_dec radius)))
    (setf (img-cat img)
          (with-open-stream (sin (sys:open-pipe cmdstr :direction :input))
            (loop for line = (read-line sin nil sin nil)
                  until (eql line sin)
                  collect line)))
    (prep-cat img)
    (gen-fast-cat-index img)
    ))

(defun gen-fast-cat-index (img)
  (let ((tbl  (make-hash-table :test #'=)))
    (um:nlet iter ((ix  -1)
                   (cat (img-ncat img)))
      (if (endp cat)
          (setf (img-ncat img) tbl)
        (let ((ix-new (floor (car (car cat)))))
          (if (= ix-new ix)
              (go-iter ix (cdr cat))
            (progn
              (loop for jx from (1+ ix) to ix-new do
                      (setf (gethash jx tbl) cat))
              (go-iter ix-new (cdr cat)))
            ))
        ))
    ))

(defun prep-cat (img)
  (um:nlet iter ((cat (img-cat img)))
    (when cat
      (let ((str (car cat)))
        (if (or (zerop (length str))
                (char= (char str 0) #\#))
            (go-iter (cdr cat))
          ;; else - start of catalog proper
          (let ((hdr1  str)
                (hdr2  (cadr cat)))
            (um:nlet inner ((cat (nthcdr 3 cat))
                            (ans nil))
              (if (endp cat)
                  (setf (img-ncat img) (stable-sort
                                        (sort ans #'< :key #'cadr)
                                        #'< :key #'car))
                (let+ ((str   (car cat)))
                  (cond ((plusp (length str))
                         (let+ ((strs  (um:split-string str :delims '(#\|)))
                                ((rastr decstr _ _ _ gmagstr . _) strs)
                                (ra  (read-from-string rastr))
                                (dec (read-from-string decstr))
                                (gmag (read-from-string gmagstr)))
                           (go-inner (cdr cat) (cons `(,ra ,dec ,gmag) ans))
                           ))
                        (t
                         (go-inner (cdr cat) ans))
                        ))))
            ))))))
                         
(defun find-star-in-cat (img star)
  (let+ ((ra     (star-ra star))
         (dec    (star-dec star))
         (cdec   (cos (dtor dec)))
         (tol    1.5) ;; arcsec
         (cat    (gethash (floor (- ra (/ tol 3600 cdec))) (img-ncat img)))
         (start  (um:nlet iter ((cat cat))
                     (unless (endp cat)
                       (let+ (( (catra _ _) (car cat))
                              (dra  (* 3600 cdec (- ra catra))))
                         (if (< dra tol)
                             cat
                           (go-iter (cdr cat))
                           ))))))
    (um:nlet iter ((cat start)
                   (ans nil))
      (if (endp cat)
          ans
        (let+ ((cent (car cat))
               ((catra catdec _) cent)
               (dra (* 3600 cdec (- catra ra))))
          (if (>= dra tol)
              ans
            (let+ ((ddec (* 3600 (- catdec dec)))
                   (dr   (abs (complex dra ddec))))
              (go-iter (cdr cat)
                       (if (and (< dr tol)
                                (or (null ans)
                                    (< dr (car ans))))
                           `(,dr ,dra ,ddec ,@cent)
                         ans))
              )))
        ))))

(defun find-some-star (img)
  (dolist (star (img-stars img))
    (let ((ans (find-star-in-cat img star)))
      (when ans
        (break))
      )))


(defun find-stars-in-cat (img)
  (dolist (star (img-stars img))
    (let ((ans (find-star-in-cat img star)))
      (if ans
          (let+ (( (_ dra ddec _ _ cmag) ans))
            (setf (star-catv star) cmag
                  (star-dx star) dra
                  (star-dy star) ddec))
        (setf (star-catv star) nil
              (star-dx star) nil
              (star-dy star) nil))
      )))
                    
#|
(get-star-positions *saved-img*)
(get-catalog *saved-img*)
(find-stars-in-cat *saved-img*)
(report-stars *saved-img*)
(show-img 'img *saved-img*)
|#

;; -----------------------------------------------------------------------------
#|
(defun tst-viz ()
  (let ((s #>.end
vizquery -mime=csv <<====End
### Example of a list of Targets for vizquery usage
#######################
-source=HIP/hip_main
#-out.add=_1 asks to insert my Hipparcos number as leftmost column.
-out.add=_1
# Be as concise as possible -- don't create a new table for each star
-out.form=mini
# Get all parameters of the Hipparcos main table
-out.all
# Now comes the list of my Hipparcos numbers:
HIP=<<====myHIPsample
# In fact, I want just Hipparcos stars 1 to 10, 
# 101 to 110, 1001 to 1010, 10001 to 10010, and 100001 to 100010
1..10
101..110
1001..1010
10001..10010
100001..100010
====myHIPsample
====End
.end))
    (with-open-stream (sin (sys:open-pipe s1 :direction :input))
      (loop for line = (read-line sin nil sin nil)
            until (eql line sin)
            collect line))))

(defun tst-web ()
  ;; hash  - numeric, needed for op :create-stamp and :get-stamp
  ;; title - optional for op :create-stamp
  ;; auth code obtined from www.originstamp.org/developer
  ;;  and needs to be renewed after every 1M creates
  ;;
  ;; Invocation:
  ;;   (timestamps :get-all-stamps)
  ;;   (timestamps :get-stamp :hash xxxxxx)
  ;;   (timestamps :get-all-transactions)
  ;;   (timestamps :create-stamp :hash xxxx [:title "your title"])
  ;;
  (with-open-stream (http (comm:open-tcp-stream ;; Lispworks specific
                                                ;; "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/345/gaia2&-c=240.005064 25.967865&-c.r=1.088996&-c.u=deg&-out.form=|&-out.add=RA_ICRS,DE_ICRS&-out=pmRA&-out=pmDE&-out=Source&-out=Gmag&-out=RPmag&-out=BPmag&-out=Plx&-out=RV&-out=Rad&-out=Lum&Gmag=0..15"
                                                "vizier.cds.unistra.fr"
                                                ;; 1660
                                                ;; 443
                                                80
                                                )) ;; leave off the "http://" prefix
    (assert http)
    (flet ((crlf ()
             (princ #\return http)
             (princ #\newline http)))

      ;; get one specific TS
      #||#
      (dolist (line `(,(format nil "GET /viz-bin/asu-tsv?-source=I/345/gaia2&-c=240.005064 25.967865&-c.r=1.088996&-c.u=deg&-out.form=|&-out.add=RA_ICRS,DE_ICRS&-out=pmRA&-out=pmDE&-out=Source&-out=Gmag&-out=RPmag&-out=BPmag&-out=Plx&-out=RV&-out=Rad&-out=Lum&Gmag=0..15 HTTP/1.1")
                      "User-Agent: Lispworks Common Lisp/8.0.1"
                      "Host: vizier.cds.unistra.fr"
                      "Accept-Language: en, mi"
                      "Connection: keep-alive"))
        (princ line http)
        (crlf))
      (crlf)
      (force-output http)
      #||#
      (terpri)
      (write-string "Waiting for reply...")
      (loop for ch = (read-char-no-hang http nil :eof)
            until ch
            do (progn
                 (write-char #\.)
                 (force-output))
               (sleep 0.25)
            finally (unless (eq ch :eof)
                      (unread-char ch http))) 
      (terpri)
      (let* ((okay-line "HTTP/1.1 200 OK")
             (nok (length okay-line))
             (lines (loop for line = (read-line http nil nil)
                          while (and line
                                     (> (length line) 1))
                          collect (write-line line))))

        #|
        (unless (and (>= (length (first lines)) nok)
                     (string-equal (subseq (first lines) 0 nok)
                                   okay-line))
          (error "Can't get timestamp data"))
        |#
        
        (let* ((key "Content-Length:")
               (nkey (length key))
               (line (member key lines
                             :key (lambda (s)
                                  (subseq s 0 (min (length s) nkey)))
                             :test #'string-equal)))
          (if line
              (let* ((ndata (parse-integer (car line) :start nkey))
                     (data  (make-string ndata)))
                (read-sequence data http)
                data)
            lines))
        ))))

#|
(tst-web)
 |#
;; ----------------------------------------------------------------------------------
#|
;; full fitting, including slope - not so good, and only buys us 0.01mag improvement
;; Run PixInsight AperturePhotometry script against a G-plane image with Aperture = 3 pix radius and save CSV files.
(um:with-remembered-filename (fname "Choose a file" :photom nil
                                    :filter "*.csv")
  (let* ((data   (csv:read-file :fname fname :skip-lines 4 :hdr-lines 1 :type :semiv))
         (grp    (csv:get-group nil data))
         (gmag   (map 'vector #'read-from-string (csv:get-column "Catalog_Gmag" grp)))
         (flux3  (map 'vector #'read-from-string (csv:get-column "Flux_Ap3_1" grp)))
         (gain   80)
         (npix   1)) ;; (* pi 3 3)))
    (multiple-value-bind (xc yc slope sigma)
        (linfit:regression (map 'vector (um:rcurry #'log 10)
                                (map 'vector (um:rcurry #'/ (* npix gain)) flux3))
                           (coerce gmag 'vector)
                           1)
      (plt:plot 'flux (map 'vector (um:rcurry #'/ (* npix gain)) flux3) gmag
                :clear t
                :xlog  t
                :yrange '(16 5)
                :symbol :dot
                :alpha  0.5
                :title "Cat GMag vs Flux_Ap3"
                :xtitle "Flux [ADU]"
                :ytitle "GMag [mag]"
                )
      (print (list xc yc slope sigma))
      (plt:fplot 'flux '(0.1e3 100e6)
                 (lambda (f)
                   (let ((lf  (log f 10)))
                     (+ yc (* slope (- lf xc)))))
                 :color :red)
      (plt:draw-text 'flux
                     (format nil "Mag Offs = ~5,2F" (- yc (* slope xc)))
                     '((:data 11e3) (:data 6.3)))
      (plt:draw-text 'flux
                     (format nil "Sigma = ~5,2F" sigma)
                     '((:data 11e3) (:data 6.8)))
      (plt:draw-text 'flux
                     (format nil "Slope = ~5,2F" slope)
                     '((:data 11e3) (:data 7.3)))
      )))
    
|#
|#

(defun write-coords-file (img)
  ;; Write a coords file for submission to Astrometry.net File should
  ;; have X Y positions, one pair per line, sorted brightest to
  ;; faintest.
  ;;
  ;; Our images have pixel (0,0) at the upper left corner, to match
  ;; the view coming from the telescope. But FITS conventions have
  ;; pixel (1,1) as the first pixel (Fortran) and located at the lower
  ;; left. RA axis in FITS moves right to left, and Dec moves from
  ;; bottom to top.
  (with-open-file (*standard-output*  "~/astrometry.net.xycoords.txt"
                                      :direction :output
                                      :if-exists :supersede
                                      :if-does-not-exist :create
                                      :element-type 'character)
    (dolist (star (sort (img-stars img) #'< :Key #'star-mag))
      (princ (1+ (star-x star)))
      (princ #\space)
      (princ (1+ (star-y star)))
      (terpri))
    ))

#|
(write-coords-file *saved-img*)

2024-06-20 12:09:46,110 Testing logger
2024-06-20 12:09:46,110 Testing log.msg()
2024-06-20 12:09:46,111 Starting Job processing for Job 10774026
2024-06-20 12:09:46,118 Creating directory /home/nova/astrometry/net/data/jobs/1077/10774026
2024-06-20 12:09:46,118 submission id 10054364
2024-06-20 12:09:46,553 running: augment-xylist --out /home/nova/astrometry/net/data/jobs/1077/10774026/job.axy --scale-low 0.1 --scale-high 180.0 --scale-units degwidth --wcs wcs.fits --corr corr.fits --rdls rdls.fits --downsample 2 --objs 1000 --xylist /home/nova/astrometry/net/data/files/cached/c77/c77028ea24d5cd89a0d9791b11b49a1eca67a9dc --width 831 --height 1164 --tweak-order 2 
2024-06-20 12:09:48,653 created axy file /home/nova/astrometry/net/data/jobs/1077/10774026/job.axy
2024-06-20 12:09:48,653 command: cd /home/nova/astrometry/net/data/jobs/1077/10774026 && /home/nova/astrometry/net/solvescript.sh job-nova-10774026 job.axy >> /home/nova/astrometry/net/data/jobs/1077/10774026/log
2024-06-20 12:10:01,943 Solver completed successfully.
2024-06-20 12:10:01,943 Checking for WCS file /home/nova/astrometry/net/data/jobs/1077/10774026/wcs.fits
2024-06-20 12:10:01,944 WCS file exists
2024-06-20 12:10:02,814 Created TanWCS: <TanWCS: CRVAL (305.880158, 38.367503) CRPIX (645.000000, 636.250000) CD (-0.000240, -0.000621; -0.000620 0.000240) Image size (831.000000, 1164.000000)>
2024-06-20 12:10:02,826 SkyLocation: <SkyLocation: nside(128) healpix(55112)>
2024-06-20 12:10:02,955 Created Calibration Calibration 8341635
2024-06-20 12:10:03,394 Finished job 10774026
|#
