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
               (let+ ((u   (- (+ xp 0.5) cr_x))
                      (v   (- (+ yp 0.5) cr_y)))
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
                  (setf (img-ncat img) ans)
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
  (let+ ((ra  (star-ra star))
         (dec (star-dec star))
         (cdec (cos (dtor dec)))
         (tol  1/720)) ;; 5 arcsec
    (um:nlet iter ((cat (img-ncat img))
                   (ans nil))
      (if (endp cat)
          ans
        (let+ (( (ca cd gmag) (car cat))
               (ddec (- dec cd))
               (dra  (* cdec (- ra ca))))
          (go-iter (cdr cat)
                   (if (and (< (abs ddec) tol)
                            (< (abs dra) tol))
                       (cons `(,ddec ,dra ,(car cat)) ans)
                     ans))
          ))
      )))

(defun find-some-star (img)
  (dolist (star (img-stars img))
    (let ((ans (find-star-in-cat img star)))
      (when ans
        (break))
      )))


(defun find-stars-in-cat (img)
  (dolist (star (img-stars img))
    (um:nlet iter ((ans   (find-star-in-cat img star))
                   (rmin  nil)
                   (which nil))
      (if (endp ans)
          (if which
              (let+ (( (dx dy (_ _ gmag)) which))
                (setf (star-catv star) gmag
                      (star-dx   star) dx
                      (star-dy   star) dy))
            (setf (star-catv star) nil
                  (star-dx   star) nil
                  (star-dy   star) nil))
        ;; else
        (let+ (( (dx dy _) (car ans))
               (r (abs (complex dx dy))))
          (if (or (null rmin)
                  (< r rmin))
              (go-iter (cdr ans) r (car ans))
            (go-iter (cdr ans) rmin which))
          )))
    ))
                    
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

