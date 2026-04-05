;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "02-c-causal-compression-reduced"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("tikz" "")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "tikz"))
 :latex)

