#lang racket/base

;; DISCLAIMER: This program honestly is disconnected from reality in a lot
;; of most drastic ways. If you even consider the thought of using this kind of
;; software to drive decisions, then you are totally misinformed, to say
;; the least.

(require
  racket/string
  racket/match
  racket/flonum
  racket/function
  racket/list
  racket/format
  math/flonum
  math/matrix
  plot)

(struct day
  (ill recovered deceased quarantined observed)
  #:prefab)

(define data
  (call-with-input-file* "data-april.txt"
    (lambda (port)
      (for/list ((l (in-lines port)))
        (apply day (map string->number (string-split l)))))))

(define N 5398000) ;; total population

(struct State
  (S ;; suspectible
   E ;; exposed
   I ;; infected
   R ;; recovered
   D ;; dead
   ) #:prefab)

(struct Params
  (beta  ;; rate of exposure
   delta ;; rate of infection
   rho   ;; rate of dying (NOT death rate)
   alpha ;; fatality rate
   gamma ;; rate of recovery
   E0    ;; number of exposed people at the beginning (not observable)
   ) #:prefab #:mutable)

(define (advance state params)
  (match-let
    (((State S E I R D) state)
     ((Params beta delta rho alpha gamma _) params))
    (State
     (- S (/ (* beta I S) N))
     (+ E (/ (* beta I S) N) (- (* delta E)))
     (+ I (* delta E) (* (- alpha 1) gamma I) (- (* alpha rho I)))
     (+ R (* (- 1 alpha) gamma I))
     (+ D (* alpha rho I)))))

(define (simulate state0 params T)
  (for/fold
   ((acc '())
    (state state0)
    #:result (reverse acc))
   ((t (in-range T)))
    (values
     (cons state acc)
     (advance state params)
     )))

(define (build-state0 I R D params)
  (State (- N (Params-E0 params) I R D)
         (Params-E0 params) I R D))


(struct Deriv
  (dI ;; derivatives of I w.r.t parameters
   dR ;; derivatives of R w.r.t parameters
   dD ;; derivatives of D w.r.t parameters
   ) #:prefab)


(define (params-increment params i)
  (let* ((v (struct->vector params))
         (val (vector-ref v (add1 i))))
    (vector-set! v (add1 i) (* val 1.05))
    (values
     (- (vector-ref v (add1 i)) val)
     (apply Params (cdr (vector->list v))))))

(define (sqrt-signed a)
  (if (> a 0)
      (sqrt a)
      (- (sqrt (- a)))))

(define transform-I sqrt-signed)
(define transform-R sqrt-signed)
(define transform-D sqrt-signed)


(define (make-derivatives state0 params T)
  (let*
      ((nparams (sub1 (vector-length (struct->vector params))))
       (incremented-params (make-vector nparams))
       (increment-values (make-vector nparams))
       )
    (for ((i (in-range nparams)))
      (let-values (((inc par) (params-increment params i)))
        (vector-set! incremented-params i par)
        (vector-set! increment-values i inc)))
    (define states0 (for/vector ((ip (in-vector incremented-params)))
                  (build-state0 (State-I state0) (State-R state0) (State-D state0) ip))) 
    (for/fold
     ((acc '())
      (states states0)
      (base-state state0)
      #:result (reverse acc))
     ((t (in-range T)))
      (let ((dI (for/vector
                    ((s states)
                     (inc increment-values))
                  (/ (- (transform-I (State-I s)) (transform-I (State-I base-state))) inc)))
            (dR (for/vector
                    ((s states)
                     (inc increment-values))
                  (/ (- (transform-R (State-R s)) (transform-R (State-R base-state))) inc)))
            (dD (for/vector
                    ((s states)
                     (inc increment-values))
                  (/ (- (transform-D (State-D s)) (transform-D (State-D base-state))) inc))))
        (values (cons (Deriv dI dR dD) acc)
                (for/vector
                    ((s states)
                     (p incremented-params))
                  (cadr (simulate s p 2)))
                (cadr (simulate base-state params 2)))))))

(define ndays (length data))

(define (make-plot I R D p ndays)
  (let* ((s (build-state0 I R D p))
         (sim (simulate s p ndays)))
    (plot (list (lines (for/list ((s sim)
                                  (i (in-naturals 0)))
                         (vector i (State-I s)))
                       #:label "заболевших" #:color 3)
                (lines (for/list ((s sim)
                                  (i (in-naturals 0)))
                         (vector i (State-D s)))
                       #:label "умерших" #:color 1)
                #;(lines (for/list ((s sim)
                                  (i (in-naturals 0)))
                         (vector i (State-R s)))
                       #:label "выздоровевших" #:color 2))
          #:x-label "День от 1 апреля"
          #:y-label "Количество людей")))

(define (make-diff-plots I R D p1 p2 ndays)
  (let* ((s1 (build-state0 I R D p1))
         (s2 (build-state0 I R D p2))
         (sim1 (list->vector (simulate s1 p1 ndays)))
         (sim2 (list->vector (simulate s2 p2 ndays))))
    (plot (list
           (lines (for/list ((i (in-range 1 (vector-length sim1))))
                         (vector i (- (State-D (vector-ref sim1 i))
                                      (State-D (vector-ref sim1 (sub1 i))))))
                       #:label (format "(β = ~a)"
                                       (~r (Params-beta p1)
                                           #:precision 3)) #:color 1)
           (lines (for/list ((i (in-range 1 (vector-length sim2))))
                    (vector i (- (State-D (vector-ref sim2 i))
                                 (State-D (vector-ref sim2 (sub1 i))))))
                  #:label (format "(β = ~a)"
                                  (~r (Params-beta p2)
                                      #:precision 3)) #:color 0))
          #:x-label "День от 1 апреля"
          #:y-label "Количество умерших за день")))


(define (make-plot-cmp-I sim)
  (plot (list (lines (for/list ((s sim)
                                (i (in-naturals 0)))
                       (vector i (State-I s)))
                     #:label "модель" #:color 3)
              (lines (for/list
                         ((d data)
                          (i (in-naturals 0)))
                       (vector i (day-ill d)))
                     #:label "данные" #:color 0))
        #:x-label "День от 1 апреля"
        #:y-label "Количество заболевших"))

(define (make-plot-cmp-R sim)
  (plot (list (lines (for/list ((s sim)
                                (i (in-naturals 0)))
                       (vector i (State-R s)))
                     #:label "модель" #:color 2)
              (lines (for/list
                         ((d data)
                          (i (in-naturals 0)))
                       (vector i (day-recovered d)))
                     #:label "данные" #:color 0))
          #:x-label "День от 1 апреля"
          #:y-label "Количество выздоровевших"))

(define (make-plot-cmp-D sim)
  (plot (list (lines (for/list ((s sim)
                                (i (in-naturals 0)))
                       (vector i (State-D s)))
                     #:label "модель" #:color 1)
              (lines (for/list
                         ((d data)
                          (i (in-naturals 0)))
                       (vector i (day-deceased d)))
                     #:label "данные" #:color 0))
        #:x-label "День от 1 апреля"
        #:y-label "Количество умерших"))

;; there is no way to choose these weights properly
(define weight-I .1)
(define weight-R .1)
(define weight-D 1.0)

(define (flvector-add! v i a)
  (flvector-set! v i (fl+ (flvector-ref v i) a)))

(define (build-system sim data derivatives)
  (let* ((nparams (vector-length (Deriv-dI (car derivatives))))
         (aty (make-flvector nparams))
         (ata (make-flvector (* nparams nparams))))
    (for
        ((s sim)
         (dat data)
         (der derivatives))
      (let ((residual-I (- (transform-I (day-ill dat)) (transform-I (State-I s))))
            (residual-R (- (transform-R (day-recovered dat)) (transform-R (State-R s)) ))
            (residual-D (- (transform-D (day-deceased dat)) (transform-D (State-D s)) ))
            (dI (Deriv-dI der))
            (dR (Deriv-dR der))
            (dD (Deriv-dD der)))
        (for ((i (in-range nparams)))
          (flvector-add! aty i (* weight-I residual-I (vector-ref dI i)))
          (flvector-add! aty i (* weight-R residual-R (vector-ref dR i)))
          (flvector-add! aty i (* weight-D residual-D (vector-ref dD i))))
        (for ((i (in-range nparams)))
          (for ((j (in-range nparams)))
            (flvector-add! ata (+ i (* j nparams))
                           (* weight-I (vector-ref dI i) (vector-ref dI j)))
            (flvector-add! ata (+ i (* j nparams))
                           (* weight-R (vector-ref dR i) (vector-ref dR j)))
            (flvector-add! ata (+ i (* j nparams))
                           (* weight-D (vector-ref dD i) (vector-ref dD j)))
        ))))
    (values aty ata)))

(define (solve-system aty ata params-to-determine)
  (let* ((nparams (flvector-length aty))
         (AtA (vector->matrix nparams nparams (flvector->vector ata)))
         (AtY (vector->matrix nparams 1 (flvector->vector aty)))
         (AtA-sub (submatrix AtA params-to-determine params-to-determine))
         (AtY-sub (submatrix AtY params-to-determine '(0)))
         (X-sub (matrix-solve AtA-sub AtY-sub))
         (AtA-sub-inv (matrix-inverse AtA-sub))
         (result (make-flvector nparams))
         (result-err (make-flvector nparams))
         )
    (for ((i params-to-determine)
          (k (in-naturals 0)))
      (flvector-set! result i (matrix-ref X-sub k 0))
      (flvector-set! result-err i (sqrt-signed (matrix-ref AtA-sub-inv k k))))
    (values result result-err)))

(define (adjust-parameters p x)
  (let ((v (struct->vector p)))
    (for ((i (in-range (flvector-length x))))
      (vector-set! v (add1 i) (+ (vector-ref v (add1 i)) (flvector-ref x i))))
    (apply Params (cdr (vector->list v)))))

(define p (Params 0.2418 ;; beta = gamma * R0
                  0.1528 ;; delta -- rate of infection
                  0.200  ;; rho -- rate of dying
                  0.004  ;; fatality rate (FIXED) 
                  0.0195 ;; gamma -- rate of recovery
                  327 ;; E0
                  ))

(define p-err (make-flvector 6))

(define sim (void))

(for ((n (in-range 5)))
  (displayln p)
  ;(displayln p-err)
  #;(displayln
   (make-plot (day-ill (car data)) 
              (day-recovered (car data))
              (day-deceased (car data))
              p (length data)))
  (define s (build-state0 (day-ill (car data)) 
                        (day-recovered (car data))
                        (day-deceased (car data))
                        p))
  (define d (make-derivatives s p ndays))
  ;(displayln d)
  (set! sim (simulate s p ndays))
  (define-values (x x-err)
    (let-values (((aty ata) (build-system sim data d)))
      (solve-system aty ata
                    ;(0 3 4 5)
                    '(0 1 2 4)
                    ;'(0 2 4 5)
                    )))
  ;(displayln x)
  (set! p (adjust-parameters p (flvector-scale x 1.0)))
  (set! p-err x-err))

(plot-width 1000)
(plot-height 600)

(displayln (make-plot-cmp-I sim))

(displayln (make-plot-cmp-R sim))

(displayln (make-plot-cmp-D sim))

(displayln
   (make-plot (day-ill (car data)) 
              (day-recovered (car data))
              (day-deceased (car data))
              p 350))


(define  p2 (struct-copy Params p
                         (beta (* (Params-beta p) 2.0))))


(displayln
   (make-plot (day-ill (car data)) 
              (day-recovered (car data))
              (day-deceased (car data))
              p2 350))


(displayln
   (make-diff-plots (day-ill (car data)) 
              (day-recovered (car data))
              (day-deceased (car data))
              p p2 350))
