{-# LANGUAGE BangPatterns #-}
module Integrals ( twoElectronIntegral
                 , kineticIntegral
                 , nuclearIntegral
                 , overlapIntegral
                 , pType
                 , boys
                 )
                  where
import Point3D
import Gaussian
import qualified Numeric.GSL.Special as GSL

boys :: Double -> Double -> Double
boys 0 x
  | x > 0.0e-5 = 0.5 * sqrt (pi / x) * GSL.erf (sqrt x)
  | otherwise = 1.0

boys m x = k `seq` k / d
  where
    k = GSL.hyperg_1F1 (m + 0.5) (m + 1.5) (-x)
    d = (2 * m) + 1

factorial n = product [1..n]

centerOfProduct Gaussian {alpha = a1, center = c1}
  Gaussian {alpha = a2, center = c2} =
    ddiv (add (dmult a1 c1) (dmult a2 c2)) (a1 + a2)
{-# INLINE centerOfProduct #-}

normalizationFactor Gaussian {alpha = a, xIndex = l, yIndex = m, zIndex = n} =
    ((fromIntegral (f1 * f2 * f3)) * ((pi / (2.0*a))**1.5)
    / (fromIntegral (2^(2*(lmn))) * a^lmn)) ** (-0.5)
    where 
      f = factorial
      lmn = l + m + n
      f1 = f (2*l - 1)
      f2 = f (2*m - 1)
      f3 = f (2*n - 1)
{-# INLINE normalizationFactor #-}

pType g = pAxis g > 0

pAxis :: Gaussian -> Int
pAxis Gaussian {xIndex = l, yIndex = m, zIndex = n} 
  | l > 0 && m == 0 && n == 0 = 1
  | l == 0 && m > 0 && n == 0 = 2
  | l == 0 && m == 0 && n > 0 = 3
  | otherwise = 0
{-# INLINE pAxis #-}

pOrbitalAxis :: Gaussian -> Maybe Axis
pOrbitalAxis Gaussian {xIndex = l, yIndex = m, zIndex = n} 
  | l > 0 && m == 0 && n == 0 = Just X
  | l == 0 && m > 0 && n == 0 = Just Y
  | l == 0 && m == 0 && n > 0 = Just Z
  | otherwise = Nothing
{-# INLINE pOrbitalAxis #-}

sType g = pAxis g == 0
{-# INLINE sType #-}

δ i j
  | i == j = 1.0
  | otherwise = 0.0
{-# INLINE δ #-}

s :: Int -> Int -> Gaussian -> Gaussian -> Double
s i j a b 
  | (i,j) == (0,0) =
    n1 * n2 * ((pi / (a1 + a2)) ** 1.5 ) * exp (((-a1) * a2)/(a1+a2) * euclidean2 c1 c2)
  | j == 0 =
    (-a2)/(a1 + a2) * (c1 @@ i - c2 @@ i) * s00ab
  | i == 0 =
    (-a1)/(a1 + a2) * (c2 @@ j - c1 @@ j) * s00ab
  | otherwise = 
    (( (a1 * a2) * (c1 @@ i - c2 @@ i) * (c2 @@ j - c1 @@ j) / ((a1 + a2)^2) )
    + 0.5 * dij / (a1 + a2) ) * s00ab
  where
    s00ab = s 0 0 a b
    !dij = δ i j
    !n1 = normalizationFactor a
    !n2 = normalizationFactor b
    !a1 = alpha a
    !a2 = alpha b
    !c1 = center a
    !c2 = center b
{-# INLINE s #-}

ell :: Int -> Int -> Point3D -> Gaussian -> Gaussian -> Double
ell i j c g1 g2
  | i == 0 && j == 0 = 
    f 0 t
  | j == 0 =
    (c @@ i - p @@ i ) * f 1 t
  | i == 0 =
    (c @@ j - p @@ j) * f 1 t
  | otherwise =
    (p @@ i - c @@ i) * (p @@ j - c @@ j) * (f 2 t) 
    - ((δ i j) * (f 1 t) )/ (2.0 * (a1 + a2))
  where
    !f = boys
    t = (a1 + a2) * euclidean2 p c
    p = centerOfProduct g1 g2
    a1 = alpha g1
    a2 = alpha g2
{-# INLINE ell #-}

k :: Int -> Int -> Gaussian -> Gaussian -> Double
k i j g1 g2
  | i == 0 && j == 0 =
    3.0 * a1 * a2 / (a1 + a2) - (euclidean2 c1 c2) * (2.0 * a1^2 * a2 ^2)/ (a1 + a2)^2 
  | j == 0 =
    (-2.0) * a1 * a2^2 / (a1 + a2)^2 *(c1 @@ i - c2 @@ i)
  | i == 0 =
    (-2.0) * a2 * a1^2 / (a1 + a2)^2 *(c2 @@ j - c1 @@ j)
  | otherwise =
    (δ i j) * a1 * a2 / (a1 + a2)^2
  where
    n1 = normalizationFactor g1
    n2 = normalizationFactor g2
    !a1 = alpha g1
    !a2 = alpha g2
    !c1 = center g1
    !c2 = center g2
{-# INLINE k #-}

-- ABOVE HAS BEEN CHECKED
type FourIndex = (Int, Int, Int, Int)
type FourGaussian = (Gaussian, Gaussian, Gaussian, Gaussian)

g :: FourIndex -> FourGaussian -> Double
-- G0000
g (0,0,0,0) (g1,g2,g3,g4) = boys 0 t
  where
    !t = (s1*s2/s4) * euclidean2 p q
    p = centerOfProduct g1 g2
    q = centerOfProduct g3 g4
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4

-- Gi000 == G0i00
g (0,i,0,0) gs = g (i,0,0,0) gs
g (i,0,0,0) (g1,g2,g3,g4) = ((-s2)/s4) * (pMinQ i) * boys 1 t
  where
    !t = (s1*s2/s4) * euclidean2 p q
    p = centerOfProduct g1 g2
    q = centerOfProduct g3 g4
    pMinQ x = p @@ x - q @@ x
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4

-- G00i0 == G000i
g (0,0,i,0) gs = g (0,0,0,i) gs
g (0,0,0,i) (g1,g2,g3,g4) = (-s1)/s4 * (qMinP i) * boys 1 t
  where
    !t = (s1)*(s2)/(s4) * euclidean2 p q
    p = centerOfProduct g1 g2
    q = centerOfProduct g3 g4
    qMinP x = q @@ x - p @@ x
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4

-- Gij00 
g (i,j,0,0) (g1,g2,g3,g4) = 
  (s2/s4) * ((s2/s4)*(pMinQ i)*(pMinQ j)*(boys 2 t) - (δ i j)*(boys 1 t)/(2.0*s1))
  where
    !t = (s1*s2/s4) * euclidean2 p q
    p = centerOfProduct g1 g2
    q = centerOfProduct g3 g4
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4

-- G00ij
g (0,0,i,j) (g1,g2,g3,g4) = 
  (s1/s4) * ((s1/s4)*(qMinP i)*(qMinP j)*(boys 2 t) - (δ i j)*(boys 1 t)/(2.0*s2))
  where 
    !t = (s1*s2/s4) * euclidean2 p q
    p = centerOfProduct g1 g2
    q = centerOfProduct g3 g4
    qMinP x = q @@ x - p @@ x
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4

-- Gi0j0 == G0ij0 == G0i0j == Gi00j
g (i,0,j,0) gs = g (i,0,0,j) gs
g (0,i,j,0) gs = g (i,0,0,j) gs
g (0,i,0,j) gs = g (i,0,0,j) gs
g (i,0,0,j) (g1,g2,g3,g4) =
  (((s1*s2)/s4) * (pMinQ i) * (qMinP j)*(boys 2 t) + 0.5*(δ i j)*(boys 1 t))/s4
  where 
    !t = (s1*s2/s4) * euclidean2 p q
    p = centerOfProduct g1 g2
    q = centerOfProduct g3 g4
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4

-- Gijk0 == Gij0k
g (i,j,0,k) gs = g (i,j,k,0) gs
g (i,j,k,0) (g1,g2,g3,g4) =
    (-s2)/(s4^2)*(s1*s2/s4 * (pMinQ i)*(pMinQ j)*(qMinP k)
    * (boys 3 t) + 0.5 *((δ i j)*(pMinQ k) + (δ i k)*(pMinQ j) 
      + (δ j k) * (pMinQ i))*(boys 2 t))
  where 
    !t = (s1*s2/s4) * euclidean2 p q
    p = centerOfProduct g1 g2
    q = centerOfProduct g3 g4
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4

-- G0ijk == Gi0jk
g (0,i,j,k) gs = g (i,0,j,k) gs
g (i,0,j,k) (g1,g2,g3,g4) = 
  (-s1)/(s4^2)*(s1*s2/s4 * (pMinQ i)*(qMinP j)*(qMinP k) * (boys 3 t)
  - 0.5 *(boys 2 t) *((δ i j)*(pMinQ k) + (δ i k)*(pMinQ j) 
  + (δ j k) * (pMinQ i)) )
  where 
    !t = (s1*s2/s4) * euclidean2 p q
    p = centerOfProduct g1 g2
    q = centerOfProduct g3 g4
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4


g (i,j,k,l) (g1,g2,g3,g4) = 
    1.0/s4^2 * ((s1^2 * s2^2 / s4^2)*(pMinQ i)*(pMinQ j)*(qMinP k)
    *(qMinP l)*(boys 4 t) - 0.5*s1*s2/s4 *(boys 3 t) *
    (
     (δ i j)*(pMinQ k)*(pMinQ l)
    +(δ i k)*(pMinQ j)*(pMinQ l)
    +(δ i l)*(pMinQ j)*(pMinQ k)
    +(δ j k)*(pMinQ i)*(pMinQ l)
    +(δ j l)*(pMinQ i)*(pMinQ k)
    +(δ k l)*(pMinQ i)*(pMinQ j)
    )
    + 0.25*((δ i j)*(δ k l) 
            + (δ i k)*(δ j l)
            + (δ i l)*(δ j k))*(boys 2 t))
  where
    !t = (s1*s2/s4) * euclidean2 p q
    p = centerOfProduct g1 g2
    q = centerOfProduct g3 g4
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4
{-# INLINE g #-}

overlapIntegral :: Gaussian -> Gaussian -> Double
overlapIntegral g1 g2 
  -- <Sa|Sb>
  | sType g1 && sType g2 = 
    s 0 0 g1 g2
  -- <Pa|Sb>
  | pType g1 && sType g2 = 
    s i 0 g1 g2
  -- <Sa|Pb>
  | sType g1 && pType g2 = 
    s 0 j g1 g2
  -- <Pa|Pb>  
  | otherwise =
    s i j g1 g2
  where
    i = pAxis g1
    j = pAxis g2
{-# INLINE overlapIntegral #-}

kineticIntegral :: Gaussian -> Gaussian -> Double
kineticIntegral g1 g2
  | sType g1 && sType g2 =
    (k 0 0 g1 g2) * (s 0 0 g1 g2)
  | pType g1 && sType g2 =
    (k i 0 g1 g2) * (s 0 0 g1 g2) + (k 0 0 g1 g2) * (s i 0 g1 g2)
  | sType g1 && pType g2 =
    (k 0 j g1 g2) * (s 0 0 g1 g2) + (k 0 0 g1 g2) * (s 0 j g1 g2)
  | otherwise =
    (k i j g1 g2) * (s 0 0 g1 g2) + (k i 0 g1 g2) * (s 0 j g1 g2)
    + (k 0 j g1 g2) * (s i 0 g1 g2) + (k 0 0 g1 g2) * (s i j g1 g2)
  where
    i = pAxis g1
    j = pAxis g2
{-# INLINE kineticIntegral #-}

nuclearIntegral :: [(Double, Point3D)] -> Gaussian -> Gaussian -> Double
nuclearIntegral nuclearCharges g1 g2
  | i == 0 && j == 0 =
    let f c = s00 * l00 c in
      theta * (sum $! zipWith (*) charges (f <$> centers))
  | j == 0 =
    let f c = si0 * (l00 c) + s00 * (li0 c) in
      theta * (sum $ zipWith (*) charges (f <$> centers))
  | i == 0 =
    let f c = s0j * (l00 c) + s00 * (l0j c) in
      theta * (sum $ zipWith (*) charges (f <$> centers))
  | otherwise = 
    let f c = sij * (l00 c) + si0 * (l0j c) + s0j * (li0 c) + s00 * (lij c) in
      theta * (sum $ zipWith (*) charges (f <$> centers))
  where
    -- input of charges/centers subject to change
    (charges, centers) = unzip nuclearCharges
    !theta = (-2.0)*(sqrt (a1 + a2))/ (sqrt pi)
    l00 c = ell 0 0 c g1 g2
    li0 c = ell i 0 c g1 g2
    l0j c = ell 0 j c g1 g2
    lij c = ell i j c g1 g2
    s00 = s 0 0 g1 g2
    si0 = s i 0 g1 g2
    s0j = s 0 j g1 g2
    sij = s i j g1 g2
    a1 = alpha g1
    a2 = alpha g2
    i = pAxis g1
    j = pAxis g2
{-# INLINE nuclearIntegral #-}

-- 0 2 2 5 seems to be wrong!
twoElectronIntegral :: Gaussian -> Gaussian -> Gaussian -> Gaussian -> Double
twoElectronIntegral a b c d 
  -- <SaSb|ScSd>
  | ssss =
    factor * (s00ab) * (s00cd) * g0000
  -- <PaSb|ScSd>
  | psss =
    factor * ((si0 * s00cd * g0000) + (s00ab * s00cd * gi000))
  -- <PaPb|ScSd>
  | ppss =
    factor * ((sij * s00cd * g0000) + (si0 * s00cd * g0j00)
             + (s0j * s00cd * gi000) + (s00ab * s00cd * gij00))
  -- <PaSb|PcSd>
  | psps =
    factor * (si0 * ((s00cd * g00k0) + (sk0 * g0000))
    + s00ab * ((s00cd * gi0k0) + (sk0 * gi000)))
  -- <PaPb|PcSd>
  | ppps = 
    factor * ((sij * (sk0 * g0000 + s00cd * g00k0))
           + (si0 * (sk0 * g0j00 + s00cd * g0jk0))
           + (s0j * (sk0 * gi000 + s00cd * gi0k0))
           + (s00ab * (sk0 * gij00 + s00cd * gijk0)))
  -- <PaPb|PcPd>
  | pppp =
    factor * (sij * (skl * g0000 + sk0 * g000l + s0l * g00k0 + s00cd * g00kl) 
    + si0 * (skl * g0j00 + sk0 * g0j0l + s0l * g0jk0 + s00cd * g0jkl)
    + s0j * (skl * gi000 + sk0 * gi00l + s0l * gi0k0 + s00cd * gi0kl)
    + s00ab * (skl * gij00 + sk0 * gij0l + s0l * gijk0 + s00cd * gijkl))
  where
    ssss = sType a && sType b && sType c && sType d
    psss = pType a && sType b && sType c && sType d
    ppss = pType a && pType b && sType c && sType d
    psps = pType a && sType b && pType c && sType d
    ppps = pType a && pType b && pType c && sType d 
    pppp = pType a && pType b && pType c && pType d 
    !s00ab = s 0 0 a b
    !s00cd = s 0 0 c d
    sij = s i j a b
    skl = s k l c d
    si0 = s i 0 a b
    sk0 = s k 0 c d
    s0j = s 0 j a b
    s0l = s 0 l c d
    !factor = 2.0 * (s1 * s2 / (pi * s4))**0.5
    !g0000 = g (0,0,0,0) (a,b,c,d)
    gi000 = g (i,0,0,0) (a,b,c,d)
    g0j00 = g (0,j,0,0) (a,b,c,d)
    g00k0 = g (0,0,k,0) (a,b,c,d)
    g000l = g (0,0,0,l) (a,b,c,d)
    gij00 = g (i,j,0,0) (a,b,c,d)
    g0jk0 = g (0,j,k,0) (a,b,c,d)
    g00kl = g (0,0,k,l) (a,b,c,d)
    gijk0 = g (i,j,k,0) (a,b,c,d)
    g0jkl = g (0,j,k,l) (a,b,c,d)
    gi0k0 = g (i,0,k,0) (a,b,c,d)
    gi00l = g (i,0,0,l) (a,b,c,d)
    gi0kl = g (i,0,k,l) (a,b,c,d)
    gijkl = g (i,j,k,l) (a,b,c,d)
    g0j0l = g (0,j,0,l) (a,b,c,d)
    gij0l = g (i,j,0,l) (a,b,c,d)
    !i = pAxis a
    !j = pAxis b
    !k = pAxis c
    !l = pAxis d
    !s1 = a1 + a2
    !s2 = a3 + a4
    !s4 = s1 + s2
    !a1 = alpha a
    !a2 = alpha b
    !a3 = alpha c
    !a4 = alpha d
{-# INLINE twoElectronIntegral #-}
