;; SPDX-License-Identifier: MPL-2.0
;; Guix development environment template.
;; Usage: guix shell -D -f build/guix.scm

(use-modules (guix packages)
             (guix build-system gnu)
             (guix licenses)
             (gnu packages base)
             (gnu packages bash))

(package
  (name "knot-knot")
  (version "0.1.0")
  (source #f)
  (build-system gnu-build-system)
  (inputs (list coreutils bash))
  (synopsis "knot-knot")
  (description "knot-knot — part of the hyperpolymath ecosystem.")
  (home-page "https://github.com/metadatastician/knot-knot")
  (license (@ (guix licenses) mpl2.0)))
