module BasisFunction ( basisFunction
                     , kineticMatrix
                     , nuclearMatrix
                     , twoElectronMatrix
                     , overlapMatrix
                     , contraction
                     ) where

import Numeric.LinearAlgebra hiding (Gaussian)
import Data.Array.Repa hiding ((!), (++), sum, map, zipWith, replicate, reshape, toList)
import qualified Data.Array.Repa as Repa (map, reshape, toList)
import qualified Numeric.GSL.Special as S
import Matrix (force)
import Control.Applicative
import Data.Vector ((!))
import qualified Data.Vector as V
import qualified Point3D as P
import qualified Geometry as G
import Integrals
import Gaussian
import Element (atomicNumber)
import STO3G

allPairs a b = (,) <$> a <*> b

sto3GH center = sto3G center 1.24
sto3GHe center = sto3G center 2.0925

sto3G center zeta =
    STONG [0.444635, 0.535328, 0.154329] 
          [ Gaussian center (zeta**2 * 0.109818) 0 0 0
          , Gaussian center (zeta**2 * 0.405771) 0 0 0
          , Gaussian center (zeta**2 * 2.22766) 0 0 0]

basisFunction G.Atom {G.center = r, G.element = e}
  | n == 1 = sto3GH r
  | n == 2 = sto3GHe r
  | otherwise = sto3GH r -- Hydrogen by default
  where
    n = Element.atomicNumber e

contraction bs integral =
    sum $ map (f . unzip) (sequence $ map tuplesFromSTO bs)
    where 
      f (ds, gs) = (product ds) * (integral gs)

tuplesFromSTO STONG {contractionCoefficients = c, gaussians =g } =
    zip c g

overlapMatrix basis =
  reshape (length basis) $ fromList $ overlap <$> basis <*> basis
  where
    overlap b1 b2 = contraction [b1,b2] overlapList
    overlapList (g1:g2:[]) = overlapIntegral g1 g2

kineticMatrix basis =
  reshape (length basis) $ fromList $ kinetic <$> basis <*> basis
  where
    kinetic b1 b2 = contraction [b1, b2] kineticList
    kineticList (g1:g2:[]) = kineticIntegral g1 g2
 
-- NEED TO REWRITE THIS TO USE NEW INTEGRAL METHOD
nuclearMatrix :: G.Geometry -> [STONG] -> Matrix Double
nuclearMatrix atoms basis =
  reshape (length basis) $ fromList $ map (v zr) (allPairs basis basis)
  where 
    integral ((z, r),(b1, b2)) = nuclear z r b1 b2
    v zr bs = sum $ map integral (allPairs zr [bs])
    nuclear z r b1 b2 = contraction [b1, b2] (nuclearList z r)
    nuclearList z r (g1:g2:[]) = nuclearIntegral [(z,r)] g1 g2
    zr = zip (map (fromIntegral . G.atomicNumber) atoms) (map G.center atoms)
    (z,r) = unzip zr

twoElectronMatrix basis =
    force $ fromFunction (Z:.n:.n:.n:.n) valAtIndex
    where
      vbasis = V.fromList basis
      valAtIndex = (\(Z:.i:.j:.k:.l) -> twoElectron (vbasis ! i) (vbasis ! j) (vbasis ! k) (vbasis ! l))
      n = V.length vbasis :: Int
      twoElectron b1 b2 b3 b4 = contraction [b1,b2,b3,b4] twoElectronList
      twoElectronList (g1:g2:g3:g4:[]) = twoElectronIntegral g1 g2 g3 g4


