(in-package #:user)

(asdf:defsystem "com.ral.photometry"
  :description "photometry: Automated Aperture Photometry of Star Frames"
  :version     "1.0"
  :author      "D.McClain <dbm@refined-audiometrics.com>"
  :license     "Copyright (c) 2024 by Refined Audiometrics Laboratory, LLC. All rights reserved."
  :components  ((:file "packages")
                (:file "boxes")
                (:file "handy")
                (:file "photom")
                (:file "fits-reader")
                (:file "measure")
                (:file "fakes")
                (:file "plant+harvest"))
  :serial       t
  :depends-on   ("com.ral.useful-macros"
                 ))

