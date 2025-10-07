#!/usr/bin/env hy

(import os)
(import re)
(import shutil)
(import sys)
(import tomllib)

(setv config-file (get sys.argv 1))
(with [f (open config-file "rb")]
  (setv config (tomllib.load f)))

(setv config-data (get config "data"))
(setv config-genotypes (get config "genotypes"))
(setv config-fish (get config "fish_used"))
(setv data-files
  (lfor
    suffix
    (get config-data "files")
    (+
      (get config-data "files_prefix")
      suffix)))
(setv genotype-files
  (lfor
    suffix
    (get config-genotypes "files")
    (+
      (get config-genotypes "files_prefix")
      suffix)))
(setv fish-files
  (lfor
    suffix
    (get config-fish "files")
    (+
      (get config-fish "files_prefix")
      suffix)))

(setv date-re (re.compile r"(\d{8})T\d{6}"))
(setv group-re (re.compile r"/([a-d]+)/"))
(setv destination-directory (get sys.argv 2))
(if (>= (len sys.argv) 4)
  (setv genotype (get sys.argv 3))
  (setv genotype "unknown"))
(for [#(data-file genotype-file fish-file) (zip data-files genotype-files fish-files)]
  (print data-file)
  ; make prefix
  (setv date (.group (re.search date-re data-file) 1))
  (setv group (.group (re.search group-re data-file) 1))
  (setv prefix f"{destination-directory}{genotype}_{(cut date None 4)}-{(cut date 4 6)}-{(cut date 6 None)}_{(.upper group)}")
  ; copy/symlink files
  (shutil.copy data-file (+ prefix ".csv"))
  (shutil.copy genotype-file (+ prefix "_genotypes.csv"))
  (setv fish-destination (+ prefix "_fish.txt"))
  (if (.startswith fish-file "templates/fish_used/")
    (os.symlink (os.path.abspath fish_file) (os.path.abspath fish-destination))
    (shutil.copy fish-file fish-destination)))
