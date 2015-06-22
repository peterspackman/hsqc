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

centerOfProduct Gaussian {α = α1, center = c1}
  Gaussian {α = α2, center = c2} =
    ddiv (add (dmult α1 c1) (dmult α2 c2)) (α1 + α2)
{-# INLINE centerOfProduct #-}

normalizationFactor Gaussian {α = α, momenta = (l, m, n) } =
    ((fromIntegral (f1 * f2 * f3)) * ((pi / (2.0*α))**1.5)
    / (fromIntegral (2^(2*(lmn))) * α^lmn)) ** (-0.5)
    where 
      f = factorial
      lmn = l + m + n
      f1 = f (2*l - 1)
      f2 = f (2*m - 1)
      f3 = f (2*n - 1)

{-# INLINE normalizationFactor #-}

pType g = pAxis g > 0

pAxis :: Gaussian -> Int
pAxis Gaussian {momenta = (l,m,n)} 
  | l > 0 && m == 0 && n == 0 = 1
  | l == 0 && m > 0 && n == 0 = 2
  | l == 0 && m == 0 && n > 0 = 3
  | otherwise = 0
{-# INLINE pAxis #-}

pOrbitalAxis :: Gaussian -> Maybe Axis
pOrbitalAxis Gaussian {momenta = (l,m,n)} 
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
    n1 * n2 * ((pi / (α1 + α2)) ** 1.5 ) * exp (((-α1) * α2)/(α1+α2) * euclidean2 c1 c2)
  | j == 0 =
    (-α2)/(α1 + α2) * (c1 @@ i - c2 @@ i) * s00ab
  | i == 0 =
    (-α1)/(α1 + α2) * (c2 @@ j - c1 @@ j) * s00ab
  | otherwise = 
    (( (α1 * α2) * (c1 @@ i - c2 @@ i) * (c2 @@ j - c1 @@ j) / ((α1 + α2)^2) )
    + 0.5 * dij / (α1 + α2) ) * s00ab
  where
    s00ab = s 0 0 a b
    !dij = δ i j
    !n1 = normalizationFactor a
    !n2 = normalizationFactor b
    !α1 = α a
    !α2 = α b
    !c1 = center a
    !c2 = center b
{-# INLINE s #-}

ell :: Int -> Int -> Point3D -> Gaussian -> Gaussian -> Double
ell i j c a b
  | i == 0 && j == 0 = 
    f 0 t
  | j == 0 =
    (c @@ i - p @@ i ) * f 1 t
  | i == 0 =
    (c @@ j - p @@ j) * f 1 t
  | otherwise =
    (p @@ i - c @@ i) * (p @@ j - c @@ j) * (f 2 t) 
    - ((δ i j) * (f 1 t) )/ (2.0 * (α1 + α2))
  where
    !f = boys
    t = (α1 + α2) * euclidean2 p c
    p = centerOfProduct a b
    α1 = α a
    α2 = α b
{-# INLINE ell #-}

k :: Int -> Int -> Gaussian -> Gaussian -> Double
k i j a b
  | i == 0 && j == 0 =
    3.0 * α1 * α2 / (α1 + α2) - (euclidean2 c1 c2) * (2.0 * α1^2 * α2 ^2)/ (α1 + α2)^2 
  | j == 0 =
    (-2.0) * α1 * α2^2 / (α1 + α2)^2 *(c1 @@ i - c2 @@ i)
  | i == 0 =
    (-2.0) * α2 * α1^2 / (α1 + α2)^2 *(c2 @@ j - c1 @@ j)
  | otherwise =
    (δ i j) * α1 * α2 / (α1 + α2)^2
  where
    n1 = normalizationFactor a
    n2 = normalizationFactor b
    !α1 = α a
    !α2 = α b
    !c1 = center a
    !c2 = center b
{-# INLINE k #-}

-- ABOVE HAS BEEN CHECKED
type FourIndex = (Int, Int, Int, Int)
type FourGaussian = (Gaussian, Gaussian, Gaussian, Gaussian)

g :: FourIndex -> FourGaussian -> Double
-- G0000
g (0,0,0,0) (a,b,c,d) = boys 0 t
  where
    !t = (ζ*η/ρ) * euclidean2 p q
    p = centerOfProduct a b
    q = centerOfProduct c d
    ζ = α a + α b
    η = α c + α d
    ρ = ζ + η

-- Gi000 == G0i00
g (0,i,0,0) gs = g (i,0,0,0) gs
g (i,0,0,0) (a,b,c,d) = ((-η)/ρ) * (pMinQ i) * boys 1 t
  where
    !t = ζ * η / ρ * (euclidean2 p q)
    p = centerOfProduct a b
    q = centerOfProduct c d
    pMinQ x = p @@ x - q @@ x
    ζ = α a + α b
    η = α c + α d
    ρ = ζ + η

-- G00i0 == G000i
g (0,0,i,0) gs = g (0,0,0,i) gs
g (0,0,0,i) (a,b,c,d) = (-ζ)/ρ * (qMinP i) * boys 1 t
  where
    !t = ζ * η / ρ * (euclidean2 p q)
    p = centerOfProduct a b
    q = centerOfProduct c d
    qMinP x = q @@ x - p @@ x
    ζ = α a + α b
    η = α c + α d
    ρ = ζ + η

-- Gij00 
g (i,j,0,0) (a,b,c,d) = 
  (η/ρ) * ((η/ρ)*(pMinQ i)*(pMinQ j)*(boys 2 t) - (δ i j)*(boys 1 t)/(2.0*ζ))
  where
    !t = (ζ*η/ρ) * euclidean2 p q
    p = centerOfProduct a b
    q = centerOfProduct c d
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    ζ = α a + α b
    η = α c + α d
    ρ = ζ + η

-- G00ij
g (0,0,i,j) (a,b,c,d) = 
  (ζ/ρ) * ((ζ/ρ)*(qMinP i)*(qMinP j)*(boys 2 t) - (δ i j)*(boys 1 t)/(2.0*η))
  where 
    !t = (ζ*η/ρ) * euclidean2 p q
    p = centerOfProduct a b
    q = centerOfProduct c d
    qMinP x = q @@ x - p @@ x
    ζ = α a + α b
    η = α c + α d
    ρ = ζ + η

-- Gi0j0 == G0ij0 == G0i0j == Gi00j
g (i,0,j,0) gs = g (i,0,0,j) gs
g (0,i,j,0) gs = g (i,0,0,j) gs
g (0,i,0,j) gs = g (i,0,0,j) gs
g (i,0,0,j) (a,b,c,d) =
  (((ζ*η)/ρ) * (pMinQ i) * (qMinP j)*(boys 2 t) + 0.5*(δ i j)*(boys 1 t))/ρ
  where 
    !t = (ζ*η/ρ) * euclidean2 p q
    p = centerOfProduct a b
    q = centerOfProduct c d
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    ζ = α a + α b
    η = α c + α d
    ρ = ζ + η

-- Gijk0 == Gij0k
g (i,j,0,k) gs = g (i,j,k,0) gs
g (i,j,k,0) (a,b,c,d) =
    (-η)/(ρ^2)*(ζ*η/ρ * (pMinQ i)*(pMinQ j)*(qMinP k)
    * (boys 3 t) + 0.5 *((δ i j)*(pMinQ k) + (δ i k)*(pMinQ j) 
      + (δ j k) * (pMinQ i))*(boys 2 t))
  where 
    !t = (ζ*η/ρ) * euclidean2 p q
    p = centerOfProduct a b
    q = centerOfProduct c d
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    ζ = α a + α b
    η = α c + α d
    ρ = ζ + η

-- G0ijk == Gi0jk
g (0,i,j,k) gs = g (i,0,j,k) gs
g (i,0,j,k) (a,b,c,d) = 
  (-ζ)/(ρ^2)*(ζ*η/ρ * (pMinQ i)*(qMinP j)*(qMinP k) * (boys 3 t)
  - 0.5 *(boys 2 t) *((δ i j)*(pMinQ k) + (δ i k)*(pMinQ j) 
  + (δ j k) * (pMinQ i)))
  where 
    !t = (ζ*η/ρ) * euclidean2 p q
    p = centerOfProduct a b
    q = centerOfProduct c d
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    ζ = α a + α b
    η = α c + α d
    ρ = ζ + η


g (i,j,k,l) (a,b,c,d) = 
    1.0/ρ^2 * ((ζ^2 * η^2 / ρ^2)*(pMinQ i)*(pMinQ j)*(qMinP k)
    *(qMinP l)*(boys 4 t) - 0.5*ζ*η/ρ *(boys 3 t) *
    ((δ i j)*(pMinQ k)*(pMinQ l)
    +(δ i k)*(pMinQ j)*(pMinQ l)
    +(δ i l)*(pMinQ j)*(pMinQ k)
    +(δ j k)*(pMinQ i)*(pMinQ l)
    +(δ j l)*(pMinQ i)*(pMinQ k)
    +(δ k l)*(pMinQ i)*(pMinQ j))
    + 0.25*((δ i j)*(δ k l) 
            + (δ i k)*(δ j l)
            + (δ i l)*(δ j k))*(boys 2 t))
  where
    !t = (ζ*η/ρ) * euclidean2 p q
    p = centerOfProduct a b
    q = centerOfProduct c d
    pMinQ x = p @@ x - q @@ x
    qMinP x = q @@ x - p @@ x
    ζ = α a + α b
    η = α c + α d
    ρ = ζ + η
{-# INLINE g #-}

overlapIntegral :: Gaussian -> Gaussian -> Double
overlapIntegral a b 
  -- <Sa|Sb>
  | sType a && sType b = 
    s 0 0 a b
  -- <Pa|Sb>
  | pType a && sType b = 
    s i 0 a b
  -- <Sa|Pb>
  | sType a && pType b = 
    s 0 j a b
  -- <Pa|Pb>  
  | otherwise =
    s i j a b
  where
    i = pAxis a
    j = pAxis b
{-# INLINE overlapIntegral #-}

kineticIntegral :: Gaussian -> Gaussian -> Double
kineticIntegral a b
  | sType a && sType b =
    (k 0 0 a b) * (s 0 0 a b)
  | pType a && sType b =
    (k i 0 a b) * (s 0 0 a b) + (k 0 0 a b) * (s i 0 a b)
  | sType a && pType b =
    (k 0 j a b) * (s 0 0 a b) + (k 0 0 a b) * (s 0 j a b)
  | otherwise =
    (k i j a b) * (s 0 0 a b) + (k i 0 a b) * (s 0 j a b)
    + (k 0 j a b) * (s i 0 a b) + (k 0 0 a b) * (s i j a b)
  where
    i = pAxis a
    j = pAxis b
{-# INLINE kineticIntegral #-}

nuclearIntegral :: [(Double, Point3D)] -> Gaussian -> Gaussian -> Double
nuclearIntegral nuclearCharges a b
  | (i,j) == (0,0) =
    let f c = s00 * l00 c in
      θ * (sum $! zipWith (*) charges (f <$> centers))
  | j == 0 =
    let f c = si0 * (l00 c) + s00 * (li0 c) in
      θ * (sum $ zipWith (*) charges (f <$> centers))
  | i == 0 =
    let f c = s0j * (l00 c) + s00 * (l0j c) in
      θ * (sum $ zipWith (*) charges (f <$> centers))
  | otherwise = 
    let f c = sij * (l00 c) + si0 * (l0j c) + s0j * (li0 c) + s00 * (lij c) in
      θ * (sum $ zipWith (*) charges (f <$> centers))
  where
    -- input of charges/centers subject to change
    (charges, centers) = unzip nuclearCharges
    !θ = (-2.0)*(sqrt (α1 + α2))/ (sqrt pi)
    l00 c = ell 0 0 c a b
    li0 c = ell i 0 c a b
    l0j c = ell 0 j c a b
    lij c = ell i j c a b
    s00 = s 0 0 a b
    si0 = s i 0 a b
    s0j = s 0 j a b
    sij = s i j a b
    α1 = α a
    α2 = α b
    i = pAxis a
    j = pAxis b
{-# INLINE nuclearIntegral #-}

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
    !factor = 2.0 * (ζ * η / (pi * ρ))**0.5
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
    !ζ = α a + α b
    !η = α c + α d
    !ρ = ζ + η
{-# INLINE twoElectronIntegral #-}
