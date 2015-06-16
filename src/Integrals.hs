module Integrals ( twoElectronIntegral
                 , kineticIntegral
                 , nuclearIntegral
                 , overlapIntegral
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

normalizationFactor Gaussian {alpha = a, xIndex = l, yIndex = m, zIndex = n} =
    ((fromIntegral (f1 * f2 * f3)) * ((pi / (2.0*a))**1.5)
    / (fromIntegral (2^(2*(lmn))) * a^lmn)) ** (-0.5)
    where 
      f = factorial
      lmn = l + m + n
      f1 = f (2*l - 1)
      f2 = f (2*m - 1)
      f3 = f (2*n - 1)

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

delta i j
  | i == j = 1.0
  | otherwise = 0.0
{-# INLINE delta #-}

s :: Int -> Int -> Gaussian -> Gaussian -> Double
s i j g1 g2 
  | i == 0 && j == 0 =
    n1 * n2 * ((pi / (a1 + a2))**1.5 )* exp (((-a1) * a2)/(a1+a2) * euclidean2 c1 c2)
  | j == 0 =
    (-a2)/(a1 + a2) * (c1 @@ i - c2 @@ i) * s 0 0 g1 g2
  | i == 0 =
    (-a1)/(a1 + a2) * (c2 @@ j - c1 @@ j) * s 0 0 g1 g2
  | otherwise = 
    ((delta i j) / (2.0 * (a1 + a2)) + (a1 * a2)/((a1 + a2)^2) 
    * (c1 @@ i - c2 @@ i) *(c2 @@ j - c1 @@ j)) * s 0 0 g1 g2
  where
    n1 = normalizationFactor g1
    n2 = normalizationFactor g2
    a1 = alpha g1
    a2 = alpha g2
    c1 = center g1
    c2 = center g2

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
    - (delta i j) * (f 1 t) / (2.0 * (a1 + a2))
  where
    f = boys
    t = (a1 + a2) * euclidean2 p c
    p = centerOfProduct g1 g2
    a1 = alpha g1
    a2 = alpha g2

k :: Int -> Int -> Gaussian -> Gaussian -> Double
k i j g1 g2
  | i == 0 && j == 0 =
    3.0 * a1 * a2 / (a1 + a2) - (euclidean2 c1 c2) * (2.0 * a1^2 * a2 ^2)/ (a1 + a2)^2 
  | j == 0 =
    (-2.0) * a1 * a2^2 / (a1 + a2)^2 *(c1 @@ i - c2 @@ i)
  | i == 0 =
    (-2.0) * a2 * a1^2 / (a1 + a2)^2 *(c2 @@ i - c1 @@ i)
  | otherwise =
    (delta i j) * a1 * a2 / (a1 + a2)^2
  where
    n1 = normalizationFactor g1
    n2 = normalizationFactor g2
    a1 = alpha g1
    a2 = alpha g2
    c1 = center g1
    c2 = center g2

g :: Int -> Int -> Int -> Int -> Gaussian -> Gaussian -> Gaussian -> Gaussian -> Double
g i j k l g1 g2 g3 g4
  -- G0000
  | (i,j,k,l) == (0,0,0,0) = 
    f 0 t
  -- Gi000 == G0j00
  | (j,k,l) == (0,0,0) =
    (-s2)/s4 * (pMinQ i) * f 1 t
  | (i,k,l) == (0,0,0) =
    (-s2)/s4 * (pMinQ j) * f 1 t
  -- G00k0 == G000l
  | (i,j,l) == (0,0,0) =
    (-s1)/s4 * (qMinP k) * f 1 t
  | (i,j,k) == (0,0,0) =
    (-s1)/s4 * (qMinP l) * f 1 t
  -- Gij00 == G00kl
  | (k,l) == (0,0) =
    s2/s4 * ((s2/s4)*(pMinQ i)*(pMinQ j)*(f 2 t) - (delta i j)*(f 1 t)/(2.0*s1))
  | (i,j) == (0,0) =
    s1/s4 * ((s1/s4)*(qMinP k)*(qMinP l)*(f 2 t) - (delta k l)*(f 1 t)/(2.0*s2))
  -- Gi0k0 == G0jk0 == G0j0l == Gi00l
  | (j,l) == (0,0) =
    1.0/s4 *((s1*s2)/s4 * (pMinQ i) * (pMinQ k)*(f 2 t)-(delta i k)*(f 1 t)/(2.0))
  | (i,l) == (0,0) =
    1.0/s4 *((s1*s2)/s4 * (pMinQ j) * (pMinQ k)*(f 2 t)-(delta i k)*(f 1 t)/(2.0))
  | (i,k) == (0,0) =
    1.0/s4 *((s1*s2)/s4 * (pMinQ j) * (pMinQ l)*(f 2 t)-(delta i k)*(f 1 t)/(2.0))
  | (j,k) == (0,0) = 
    1.0/s4 *((s1*s2)/s4 * (pMinQ i) * (pMinQ l)*(f 2 t)-(delta i k)*(f 1 t)/(2.0))
  -- Gijk0 == Gij0l
  | l == 0 =
    (-s2)/(s4^2)*(s1*s2/s4 * (pMinQ i)*(pMinQ j)*(pMinQ k)
    * (f 3 t) + 0.5 *((delta i j)*(pMinQ k) + (delta i k)*(pMinQ j) 
      + (delta j k) * (pMinQ i))*(f 2 t))
  | k == 0 = 
    (-s1)/(s4^2)*(s1*s2/s4 * (pMinQ i)*(pMinQ j)*(pMinQ l )
    * (f 3 t) + 0.5 *((delta i j)*(pMinQ l) + (delta i l)*(pMinQ j) 
      + (delta j l) * (pMinQ i))*(f 2 t))
  -- G0jkl == Gi0kl
  | j == 0 =
    (-s1)/(s4^2)*(s1*s2/s4 * (pMinQ i)*(pMinQ k)*(pMinQ l )
    * (f 3 t) - 0.5 *((delta i k)*(pMinQ l) + (delta i l)*(pMinQ k) 
      + (delta k l) * (pMinQ i))*(f 2 t))
  | i == 0 = 
    (-s1)/(s4^2)*(s1*s2/s4 * (pMinQ j)*(pMinQ k)*(pMinQ l )
    * (f 3 t) - 0.5 *((delta j k)*(pMinQ l) + (delta j l)*(pMinQ k) 
      + (delta k l) * (pMinQ j))*(f 2 t))
  -- Gijkl
  | otherwise =
    1.0/s4^2 * ((s1^2 * s2^2 / s4^2)*(pMinQ i)*(pMinQ j)*(pMinQ k)
    *(pMinQ l)*(f 4 t) - (s1*s2)/(2.0*s4)*(
     (delta i j)*(pMinQ k)*(pMinQ l)
    +(delta i k)*(pMinQ j)*(pMinQ l)
    +(delta i l)*(pMinQ j)*(pMinQ k)
    +(delta k l)*(pMinQ i)*(pMinQ j))*(f 3 t)
    + 0.25*((delta i j)*(delta k l) + (delta i l)*(delta j k))*(f 2 t))
  where
    f = boys
    t = (a1 + a2)*(a3 + a4)/(a1 + a2 + a3 + a4) * euclidean2 p q
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

nuclearIntegral :: [(Double, Point3D)] -> Gaussian -> Gaussian -> Double
nuclearIntegral nuclearCharges g1 g2
  | i == 0 && j == 0 =
    let f c = s00 * l00 c in
      theta * (sum $ zipWith (*) charges (f <$> centers))
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
    theta = (-2.0)*(sqrt (a1 + a2))/ (sqrt pi)
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


twoElectronIntegral :: Gaussian -> Gaussian -> Gaussian -> Gaussian -> Double
twoElectronIntegral g1 g2 g3 g4 
  -- <SaSb|ScSd>
  | ssss =
    factor * (s00ab) * (s00cd) * g0000
  -- <PaSb|ScSd>
  | psss =
    factor * (si0 * s00cd * g0000 + s00ab * s00cd * gi000)
  -- <PaSb|PcSd>
  | psps =
    factor * (si0 * s00cd * g00k0 + si0 * sk0 * g0000 
    + s00ab * s00cd * gi0k0 + s00ab * sk0 * gi000)
  -- <PaPb|PcSd>
  | ppps = 
    factor * (sij * (sk0 * g0000 + s00cd * g00k0)
    + si0 * (sk0 * g0j00 + s00cd * g0jk0)
    + s0j * (sk0 * gi000 + s00cd * gi0k0)
    + s00ab * (sk0 * gij00 + s00cd * gijk0 ))
  -- <PaPb|PcPd>
  | otherwise =
    factor * (sij * (skl * g0000 + sk0 * g000l + s0l * g00k0 + s00cd * g00kl) 
    + si0 * (skl * g0j00 + sk0 * g0j0l + s0l * g0jk0 + s00cd * g0jkl)
    + s0j * (skl * gi000 + sk0 * gi00l + s0l * gi0k0 + s00cd * gi0kl)
    + s00ab * (skl * gij00 + sk0 * gij0l + s0l * gijk0 + s00cd * gijkl))
  where
    ssss = sType g1 && sType g2 && sType g3 && sType g4
    psss = pType g1 && sType g2 && sType g3 && sType g4
    psps = pType g1 && sType g2 && pType g3 && sType g4
    ppps = pType g1 && pType g2 && pType g3 && sType g4 
    s00ab = s 0 0 g1 g2
    s00cd = s 0 0 g3 g4
    sij = s i j g1 g2
    skl = s k l g3 g4
    si0 = s i 0 g1 g2
    sk0 = s k 0 g3 g4
    s0j = s 0 j g1 g2
    s0l = s 0 l g3 g4
    factor = 2.0 * (s1 * s2 / (pi * s4))**0.5
    g0000 = g 0 0 0 0 g1 g2 g3 g4
    gi000 = g i 0 0 0 g1 g2 g3 g4
    g0j00 = g 0 j 0 0 g1 g2 g3 g4
    g00k0 = g 0 0 k 0 g1 g2 g3 g4
    g000l = g 0 0 0 l g1 g2 g3 g4
    gij00 = g i j 0 0 g1 g2 g3 g4
    g0jk0 = g 0 j k 0 g1 g2 g3 g4
    g00kl = g 0 0 k l g1 g2 g3 g4
    gijk0 = g i j k 0 g1 g2 g3 g4
    g0jkl = g 0 j k l g1 g2 g3 g4
    gi0k0 = g i 0 k 0 g1 g2 g3 g4
    gi00l = g i 0 0 l g1 g2 g3 g4
    gi0kl = g i 0 k l g1 g2 g3 g4
    gijkl = g i j k l g1 g2 g3 g4
    g0j0l = g 0 j 0 l g1 g2 g3 g4
    gij0l = g i j 0 l g1 g2 g3 g4
    i = pAxis g1
    j = pAxis g2
    k = pAxis g3
    l = pAxis g4
    s1 = a1 + a2
    s2 = a3 + a4
    s4 = s1 + s2
    a1 = alpha g1
    a2 = alpha g2
    a3 = alpha g3
    a4 = alpha g4


