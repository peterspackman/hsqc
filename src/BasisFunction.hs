module BasisFunction ( basisFunction
                     , kineticMatrix
                     , nuclearMatrix
                     , twoElectronMatrix
                     , overlapMatrix
                     ) where

import Numeric.LinearAlgebra
import Data.Array.Repa hiding ((++), sum, map, zipWith, replicate, reshape, toList)
import qualified Data.Array.Repa as Repa (map, reshape, toList)
import qualified Numeric.GSL.Special as S
import Control.Applicative

allPairs a b = (,) <$> a <*> b

class GaussianBasisFunction a where
    overlap :: a -> a -> Double
    nuclearAttraction :: Int -> Double -> a -> a -> Double
    kineticEnergy :: a -> a -> Double
    twoElectron :: a -> a -> a -> a -> Double

data Gaussian1S = 
  Gaussian1S { alpha :: Double
             , center :: Double
             } deriving Show

data STONG =
    STONG { d :: [Double] -- contraction coefficients
          , g :: [Gaussian1S] -- primative gaussian functions
          } deriving Show

instance Eq Gaussian1S where
    (==) a b =
      (alpha a == alpha b) && (center a == center b)

instance Eq STONG where
    (==) a b =
      (g a == g b) && (d a == d b)


sto3GH center = sto3G center 1.24
sto3GHe center = sto3G center 2.0925

sto3G center zeta =
    STONG [0.444635, 0.535328, 0.154329] 
          [ Gaussian1S (zeta**2 * 0.109818) center
          , Gaussian1S (zeta**2 * 0.405771) center
          , Gaussian1S (zeta**2 * 2.22766) center]

basisFunction center atomicNumber
  | atomicNumber == 1 = sto3GH center
  | atomicNumber == 2 = sto3GHe center
  | otherwise = sto3GH center -- Hydrogen by default


normalizationTerm alpha beta = 
    (2 * alpha / pi) ** (0.75) * (2 * beta / pi) ** (0.75)

term2 :: Double -> Double -> Double
term2 alpha beta =
  (pi / (alpha + beta)) ** (1.5)


term3 a b r1 r2 =
  exp ( abs (r1-r2)^2 * ((-a)*b)/(a+ b))


instance GaussianBasisFunction Gaussian1S where
    -- OVERLAP INTEGRAL
    overlap a b = 
      if a == b
        then 1.0
      else i * n * (term3 al be r1 r2)
      where
        al = alpha a
        be = alpha b
        r1 = center a
        r2 = center b
        n = normalizationTerm al be
        i = term2 al be
    -- NUCLEAR ATTRACTION
    nuclearAttraction z r g1 g2 =
      if abs(t) < 0.0000001
        then x
      else y
      where
        al = alpha g1
        be = alpha g2
        r1 = center g1
        r2 = center g2
        rp = (al * r1 + be * r2) / (al + be)
        n = normalizationTerm al be
        x = (n * (-2.0) * pi / (al + be) * (term3 al be r1 r2) * (fromIntegral z))
        t = (al + be)*abs(rp - r)^2
        y = (x * 0.5 * sqrt (pi / t) * S.erf (sqrt t))
    -- KINETIC ENERGY INTEGRAL
    kineticEnergy g1 g2 =
      n * x * y * (3.0 - 2.0 * al * be/((al + be)/((abs (r1 - r2))^2))) * (al * be / (al + be))
      where
        al = alpha g1
        be = alpha g2
        c = al * be
        r1 = center g1
        r2 = center g2
        n = normalizationTerm al be
        x = term2 al be
        y = term3 al be r1 r2

    -- TWO ELECTRON INT
    twoElectron g1 g2 g3 g4 =
      if (abs t) < 0.00000001
        then (x)
      else y
      where
        a = alpha g1
        b = alpha g2
        gamma = alpha g3
        del = alpha g4
        ra = center g1
        rb = center g2
        rc = center g3
        rd = center g4
        rp = (a*ra + b*rb)/(a + b)
        rq = (gamma*rc + del*rd)/(gamma + del)
        t = ((abs (rp - rq))**2) * (a+b) * (gamma+del) / (a+b+gamma+del)

        x = (2 * a/pi)**0.75 * (2 * b/pi)**0.75
           *(2 * gamma/pi)**0.75 * (2 * del/pi)**0.75
           *(2 * pi**2.5)
           /((a + b)*(gamma+del)* sqrt (a + b + gamma + del))
           * exp ((-a) * b/(a + b) * (abs (ra -rb))**2 - gamma*del/(gamma+del)*(abs (rc - rd))**2)

        y = (x * 0.5 * sqrt (pi / t) * S.erf (sqrt t))

tupleSTO b =
    zip (d b) (g b)

overlapList (g1:g2:[]) = if g1 == g2 then 1.0 else overlap g1 g2
kineticList (g1:g2:[]) = kineticEnergy g1 g2
twoElectronList (g1:g2:g3:g4:[]) = twoElectron g1 g2 g3 g4
nuclearList z r (g1:g2:[]) = nuclearAttraction z r g1 g2

contraction bs integral =
    sum $ map (f . unzip) (sequence $ map tupleSTO bs)
    where 
      f (ds, gs) = (product ds) * (integral gs)

allEq xs = and $ map (== head xs) (tail xs)

overlapb b1 b2 = contraction [b1,b2] overlapList
kineticb b1 b2 = contraction [b1, b2] kineticList
twoElectronb b1 b2 b3 b4 = contraction [b1,b2,b3,b4] twoElectronList
kineticIntegral bs = contraction bs kineticList
nuclearIntegral z r b1 b2 = contraction [b1, b2] (nuclearList z r)

overlapMatrix basis =
  reshape (length basis) $ fromList $ overlapb <$> basis <*> basis

kineticMatrix basis =
  reshape (length basis) $ fromList $ kineticb <$> basis <*> basis
 
nuclearMatrix zr basis =
    reshape (length basis) $ fromList $ map (v zr) (allPairs basis basis)
  where 
    integral ((z, r),(b1, b2)) = nuclearIntegral z r b1 b2
    v zr bs = sum $ map integral (allPairs zr [bs])

twoElectronMatrix basis =
    fromListUnboxed (Z:.n:.n:.n:.n) t
    where
      t = twoElectronb <$> basis <*> basis <*> basis <*> basis
      n = length basis :: Int
