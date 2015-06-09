module BasisFunction ( basisFunction
                     , kineticMatrix
                     , nuclearMatrix
                     , twoElectronMatrix
                     , twoElectronMatrix'
                     , overlapMatrix
                     ) where

import Numeric.LinearAlgebra
import Data.Array.Repa hiding ((!), (++), sum, map, zipWith, replicate, reshape, toList)
import qualified Data.Array.Repa as Repa (map, reshape, toList)
import qualified Numeric.GSL.Special as S
import Matrix (force)
import Control.Applicative
import Data.Vector ((!))
import qualified Data.Vector as V
import qualified Point3D as P

allPairs a b = (,) <$> a <*> b

class GaussianBasisFunction a where
    overlap :: a -> a -> Double
    nuclearAttraction :: Int -> P.Point3D -> a -> a -> Double
    kineticEnergy :: a -> a -> Double
    twoElectron :: a -> a -> a -> a -> Double

data Gaussian1S = 
  Gaussian1S { alpha :: Double
             , center :: P.Point3D
             } deriving Show

data STONG =
    STONG { coefficients :: [Double] -- contraction coefficients
          , gaussians :: [Gaussian1S] -- primative gaussian functions
          } deriving Show

instance Eq Gaussian1S where
    (==) a b =
      (alpha a == alpha b) && (center a == center b)

instance Eq STONG where
    (==) a b =
      (gaussians a == gaussians b) && (coefficients a == coefficients b)


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

contraction bs integral =
    sum $ map (f . unzip) (sequence $ map tuplesFromSTO bs)
    where 
      f (ds, gs) = (product ds) * (integral gs)


normalizationTerm alpha beta = 
    (2 * alpha / pi) ** (0.75) * (2 * beta / pi) ** (0.75)

term2 :: Double -> Double -> Double
term2 alpha beta =
  (pi / (alpha + beta)) ** (1.5)


term3 :: Double -> Double -> P.Point3D -> P.Point3D -> Double
term3 a b r1 r2 =
  exp ( (P.euclidean2 r1 r2) * ((-a) * b)/(a + b))


instance GaussianBasisFunction Gaussian1S where
    -- OVERLAP INTEGRAL
    {-# INLINE overlap #-}
    overlap Gaussian1S {alpha = a, center = r1}
            Gaussian1S {alpha = b, center = r2}
            = 
      i * n * (term3 a b r1 r2)
      where
        n = normalizationTerm a b
        i = term2 a b
    -- NUCLEAR ATTRACTION
    {-# INLINE nuclearAttraction #-}
    nuclearAttraction z r 
      Gaussian1S {alpha = a, center = r1}
      Gaussian1S {alpha = b, center = r2}
      =
      if abs(t) < 0.0000001
        then x
      else y
      where
        rp = (P.add (a `P.dmult` r1) (b `P.dmult` r2)) `P.ddiv` (a + b)
        n = normalizationTerm a b
        x = (n * (-2.0) * pi / (a + b) * (term3 a b r1 r2) * (fromIntegral z))
        t = (a + b) * (P.euclidean2 rp r)
        y = (x * 0.5 * sqrt (pi / t) * S.erf (sqrt t))
    -- KINETIC ENERGY INTEGRAL
    {-# INLINE kineticEnergy#-}
    kineticEnergy 
      Gaussian1S {alpha = a, center = r1}
      Gaussian1S {alpha = b, center = r2}
      =
      n * x * y * (3.0 - 2.0 * a * b/((a + b)/(P.euclidean2 r1 r2))) * (a * b / (a + b))
      where
        n = normalizationTerm a b
        x = term2 a b
        y = term3 a b r1 r2

    -- TWO ELECTRON INT
    {-# INLINE twoElectron #-}
    twoElectron 
      Gaussian1S {alpha = a, center = ra}
      Gaussian1S {alpha = b, center = rb}
      Gaussian1S {alpha = c, center = rc}
      Gaussian1S {alpha = d, center = rd}
      =
      if (abs t) < 0.00000001
        then (x)
      else y
      where
        rp = (P.add (a `P.dmult` ra) (b `P.dmult` rb)) `P.ddiv` (a + b)
        rq = (P.add (c `P.dmult` rc) (d `P.dmult` rd)) `P.ddiv` (c + d)
        t = (P.euclidean2 rp rq) * (a+b) * (c +d) / (a + b+ c + d)

        x = (2 * a/pi)**0.75 * (2 * b/pi)**0.75
          * (2 * c/pi)**0.75 * (2 * d/pi)**0.75
          * (2 * pi**2.5)
          / ((a + b)*(c + d) * sqrt (a + b + c + d))
          * exp ((-a) * b/(a + b) * (P.euclidean2 ra rb) - c*d/(c+d)*(P.euclidean2 rc rd))

        y = (x * 0.5 * sqrt (pi / t) * S.erf (sqrt t))

tuplesFromSTO STONG {coefficients = c, gaussians =g } =
    zip c g

overlapMatrix basis =
  reshape (length basis) $ fromList $ overlapIntegral <$> basis <*> basis
  where
    overlapIntegral b1 b2 = contraction [b1,b2] overlapList
    overlapList (g1:g2:[])
      | g1 == g2 = 1.0
      | otherwise = overlap g1 g2


kineticMatrix basis =
  reshape (length basis) $ fromList $ kineticIntegral <$> basis <*> basis
  where
    kineticIntegral b1 b2 = contraction [b1, b2] kineticList
    kineticList (g1:g2:[]) = kineticEnergy g1 g2
 

nuclearMatrix zr basis =
  reshape (length basis) $ fromList $ map (v zr) (allPairs basis basis)
  where 
    integral ((z, r),(b1, b2)) = nuclearIntegral z r b1 b2
    v zr bs = sum $ map integral (allPairs zr [bs])
    nuclearIntegral z r b1 b2 = contraction [b1, b2] (nuclearList z r)
    nuclearList z r (g1:g2:[]) = nuclearAttraction z r g1 g2


twoElectronMatrix basis =
  fromListUnboxed (Z:.n:.n:.n:.n) t
  where
    t = twoElectronIntegral <$> basis <*> basis <*> basis <*> basis
    n = length basis :: Int
    twoElectronIntegral b1 b2 b3 b4 = contraction [b1,b2,b3,b4] twoElectronList
    twoElectronList (g1:g2:g3:g4:[]) = twoElectron g1 g2 g3 g4

twoElectronMatrix' basis =
    force $ fromFunction (Z:.n:.n:.n:.n) valAtIndex
    where
      vbasis = V.fromList basis
      valAtIndex = (\(Z:.i:.j:.k:.l) -> twoElectronIntegral (vbasis ! i) (vbasis ! j) (vbasis ! k) (vbasis ! l))
      n = V.length vbasis :: Int
      twoElectronIntegral b1 b2 b3 b4 = contraction [b1,b2,b3,b4] twoElectronList
      twoElectronList (g1:g2:g3:g4:[]) = twoElectron g1 g2 g3 g4


