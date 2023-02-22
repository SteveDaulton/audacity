;; spec-plot.lsp -- spectral plot function
;;
;; Roger B. Dannenberg, May 2016
;;

(setf *spec-plot-bw* 8000.0) ;; highest frequency to plot (default)
(setf *spec-plot-res* 20.0) ;; bin size (default)
(setf *spec-plot-db* nil) ;; plot dB? (default)

;; We want to allow round-number bin-sizes so plot will be more readable
;; Assuming 20Hz as an example, the FFT size would have to be
;; 44100/20 = 2205, but that's not a power of 2, so we should resample
;; the signal down so that the FFT size is 2048 (or up to 4096). This
;; would result in sample rates of 2048*20 = 40960 or 81120. We should
;; pick the smaller one if it is at least 2x *spec-plot-bw*.

(setf *spec-plot-verbose* nil)  ;; help with debugging/testing
(defmacro sp-display (&rest p)
  (if *spec-plot-verbose* `(display ,@p)))

(defun spec-plot (sound &key (offset 0)
                             (dur nil)
                             (res *spec-plot-res*)
                             (bw *spec-plot-bw*)
                             (db *spec-plot-db*))
  (ny:typecheck (and (not (soundp sound)) (not (stringp sound)))
    (ny:error "SPEC-PLOT" 1 '((SOUND STRING) nil) sound))
  (ny:typecheck (not (or (null offset) (numberp offset)))
    (ny:error "SPEC-PLOT" 0 '((NUMBER NULL) nil) offset))
  (ny:typecheck (not (or (null dur) (numberp dur)))
    (ny:error "SPEC-PLOT" 0 '((NUMBER NULL) nil) dur))
  (ny:typecheck (not (or (null res) (numberp res)))
    (ny:error "SPEC-PLOT" 0 '((NUMBER) nil) res))
  (ny:typecheck (not (or (null bw) (numberp bw)))
    (ny:error "SPEC-PLOT" 0 '((NUMBER) nil) bw))
  (let (sa fft-size power2 filename srate samps)
    (cond ((stringp sound)
           (setf filename sound)
           (setf sound (s-read filename))
           (cond ((arrayp sound)
                  (format t "WARNING: ~A contains a ~A-channel sound. ~A~%"
                          filename (length sound)
                          "Taking FFT of first channel.")))))
    (setf srate (snd-srate sound))
    (setf fft-size (/ srate res))
    ;; limit fft size to between 4 and 2^16
    (cond ((> fft-size 65536)
           (format t "Warning: fft-size reduced from ~A to 65536~%" len)
           (setf fft-size 65536))
          ((< fft-size 8)
           (format t "Warning: fft-size increased from ~A to 8~%" len)
           (setf fft-size 8))
          (t
           (setf fft-size (round fft-size))))
    (setf power2 8) ;; find integer size for FFT
    (while (< power2 fft-size)
      (setf power2 (* 2 power2)))
    (sp-display "spec-plot" fft-size srate dur (/ fft-size srate))
    ;; now power2 >= fft-size. Example: srate=1000, res=10, we need 100 pt FFT
    ;;     power2 = 128, so resample to srate 128 * 10 = 1280
    (cond ((> power2 fft-size) ;; not equal, must resample
           (setf srate (float (* power2 res)))
           (setf sound (snd-resample sound srate))
           (sp-display "resample to" srate)
           (setf fft-size power2)))
    ;; we only need fft-dur samples, but allow an extra second just to
    ;; avoid any rounding errors. We need to be careful if sound is shorter
    ;; than FFT frame -- we want to window the existing sound samples.
    (cond (dur
           (setf samps (round (* dur srate)))
           (setf samps (min samps fft-size)))
          (t
           (setf samps fft-size)))
    (setf dur (/ samps srate))
    (setf sound (extract offset (+ offset dur) sound))
    ;; now set samps to minimum of samps and length of sound
    (setf samps (snd-length sound samps))
    (setf dur (/ samps srate))
    (sp-display "before stretch" dur samps fft-size srate)
    ;(s-plot sound (* 2 dur)) (read)
    (setf sound (mult (stretch dur (sound (hamming-window samps)))
                      sound))
    (sp-display "after stretch")
    ;(s-plot sound (* 2 dur)) (read)
    ;; finally do the FFT
    (setf sa (sa-init :resolution res :input sound :window :none))
    (format t "spec-plot: ") (sa-info sa)
    (setf mag (sa-magnitude (sa-next sa)))
    (setf mag (snd-from-array 0 (/ 1.0 res) mag))
    (if db (setf mag (linear-to-db mag)))
    (s-plot mag bw (round (/ (float bw) res)))))


(defun spec-print (file sound &key (offset 0)
                              (dur nil)
                              (res *spec-plot-res*)
                              (bw *spec-plot-bw*)
                              (threshold nil))
  (ny:typecheck (and (not (filep file)) (not (eq file t)))
    (ny:error "SPEC-PRINT" 1 '((FILE) nil) file))
  (ny:typecheck (and (not (soundp sound)) (not (stringp sound)))
    (ny:error "SPEC-PRINT" 2 '((SOUND STRING) nil) sound))
  (ny:typecheck (not (or (null offset) (numberp offset)))
    (ny:error "SPEC-PRINT" 0 '((NUMBER NULL) nil) offset))
  (ny:typecheck (not (or (null dur) (numberp dur)))
    (ny:error "SPEC-PRINT" 0 '((NUMBER NULL) nil) dur))
  (ny:typecheck (not (or (null res) (numberp res)))
    (ny:error "SPEC-PRINT" 0 '((NUMBER) nil) res))
  (ny:typecheck (not (or (null bw) (numberp bw)))
    (ny:error "SPEC-PRINT" 0 '((NUMBER) nil) bw))
  (ny:typecheck (not (or (null threshold) (numberp threshold)))
    (ny:error "SPEC-PRINT" 0 '((NUMBER) nil) threshold))
  (let (sa fft-size power2 srate samps)
    (cond ((stringp sound)
           (setf filename sound)
           (setf sound (s-read filename))
           (cond ((arrayp sound)
                  (format t "WARNING: ~A contains a ~A-channel sound. ~A~%"
                          filename (length sound)
                          "Taking FFT of first channel.")))))
    ;; this implementation is similar to spec-print, but we omit the
    ;; fancy resampling and just use the "natural" sample rate, and
    ;; let sa-init pick an FFT size to get at least res resolution.
    ;;
    ;; Nevertheless, we still handle offset and windowing so that
    ;; you can get a high resolution plot on a short audio selection
    ;; by zero-padding.
    ;;
    (setf srate (snd-srate sound))
    (setf fft-size (/ srate res))
    ;; limit fft size to between 4 and 2^16
    (cond ((> fft-size 65536)
           (format t "Warning: fft-size reduced from ~A to 65536~%" len)
           (setf fft-size 65536))
          ((< fft-size 8)
           (format t "Warning: fft-size increased from ~A to 8~%" len)
           (setf fft-size 8)))

    (setf power2 8) ;; find integer size for FFT
    (while (< power2 fft-size)
      (setf power2 (* 2 power2)))
    (setf fft-size power2)
    (sp-display "spec-print" fft-size srate dur (/ fft-size srate))
    ;; we only need fft-dur samples, but allow an extra second just to
    ;; avoid any rounding errors. We need to be careful if sound is shorter
    ;; than FFT frame -- we want to window the existing sound samples.
    (cond (dur
           (setf samps (round (* dur srate)))
           (setf samps (min samps fft-size)))
          (t
           (setf samps fft-size)))
    (setf dur (/ samps srate))
    (setf sound (extract offset (+ offset dur) sound))
    ;; now set samps to minimum of samps and length of sound
    (setf samps (snd-length sound samps))
    (setf dur (/ samps srate))
    (sp-display "before stretch" dur samps fft-size srate)
    ;(s-plot sound (* 2 dur)) (read)
    (setf sound (mult (stretch dur (sound (hamming-window samps)))
                      sound))
    (sp-display "after stretch")
    ;(s-plot sound (* 2 dur)) (read)
    ;; finally do the FFT
    (setf sa (sa-init :resolution res :input sound :window :none))
    (format t "spec-plot: ") (sa-info sa)
    (sa-print file sa (sa-normalize (sa-magnitude (sa-next sa)))
              :cutoff bw :threshold threshold)))
