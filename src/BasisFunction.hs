{-# LANGUAGE BangPatterns #-}
module BasisFunction ( Basis
                     , kineticMatrix
                     , nuclearMatrix
                     , twoElectronMatrix
                     , overlapMatrix
                     , contraction
                     , sortBasisFunctions
                     , soadDiagonal
                     ) where

import Numeric.LinearAlgebra hiding (Gaussian)
import Data.Array.Repa hiding ((!), (++), sum, map, zipWith, replicate, reshape, toList)
import qualified Data.Array.Repa as Repa (map, reshape, toList)
import Matrix (Matrix4D, force)
import Control.Applicative
import Data.Vector ((!))
import qualified Data.Vector as V
import Geometry (Geometry)
import qualified Geometry as G
import Data.List (sort, sortBy, groupBy)
import Integrals
import Shell
import Gaussian
import Element (atomicNumber, electronConfig)


contraction :: Basis -> ([Gaussian] -> Double) -> Double
contraction !bs integral =
    sum $ map (f . unzip) (sequence $ map contractionPairs bs)
    where 
      f (!ds, !gs) = (product ds) * (integral gs) 
      contractionPairs Shell {weights = c, gaussians = g} =
        zip c g

overlapMatrix :: Basis -> Matrix Double
overlapMatrix basis =
  reshape (length basis) $! fromList $ overlap <$> basis <*> basis
  where
    overlap b1 b2 = contraction [b1,b2] overlapList
    overlapList (g1:g2:[]) = overlapIntegral g1 g2

kineticMatrix :: Basis -> Matrix Double
kineticMatrix b =
  reshape (length b) $! fromList $ kinetic <$> b <*> b
  where
    kinetic b1 b2 = contraction [b1, b2] kineticList
    kineticList (g1:g2:[]) = kineticIntegral g1 g2
 
nuclearMatrix :: Geometry -> Basis -> Matrix Double
nuclearMatrix atoms b =
  reshape (length b) $! fromList $ map (v zr) (allPairs b b)
  where 
    integral ((z, r),(b1, b2)) = nuclear z r b1 b2
    v zr bs = sum $ map integral (allPairs zr [bs])
    nuclear z r b1 b2 = contraction [b1, b2] (nuclearList z r)
    nuclearList z r (g1:g2:[]) = nuclearIntegral [(z,r)] g1 g2
    !zr = zip (map (fromIntegral . G.atomicNumber) atoms) (map G.center atoms)
    (!z,!r) = unzip zr
    allPairs a b = (,) <$> a <*> b

twoElectronMatrix :: Basis -> Matrix4D Double
twoElectronMatrix basis =
    force $ fromFunction (Z:.n:.n:.n:.n) valAtIndex
    where
      !vbasis = V.fromList basis
      valAtIndex = (\(Z:.i:.j:.k:.l) -> twoElectron (vbasis ! i) (vbasis ! j) (vbasis ! k) (vbasis ! l))
      !n = V.length vbasis :: Int
      twoElectron b1 b2 b3 b4 = contraction [b1,b2,b3,b4] twoElectronList
      twoElectronList (g1:g2:g3:g4:[]) = 
        let (a:b:c:d:[]) = sortBasisFunctions g1 g2 g3 g4 in
          (twoElectronIntegral a b c d)

sortBasisFunctions :: Gaussian -> Gaussian -> Gaussian -> Gaussian -> [Gaussian]
sortBasisFunctions b1 b2 b3 b4 =
    concat $ sortBy numP (map sort [[b1, b2], [b3,b4]])
  where
    pCount s = length $ filter pType s
    numP a b
      | pCount a > pCount b = LT
      | otherwise = GT


-- smear the electrons across the shells,
-- filling up s orbitals first, then p orbitals.
-- Will need changes for d orbitals!
smearElectrons :: Int -> Basis  -> [Double]
smearElectrons ne (s:[]) = [(fromIntegral ne) * 0.5]
smearElectrons ne (s:ps) 
  | ne > 2 = 1.0:(replicate np (0.5 * remaining/(fromIntegral np)))
  | otherwise = smearElectrons ne [s] 
  where
    remaining = fromIntegral (ne - 2)
    np = length ps

soadDiagonal :: Basis -> [Double]
soadDiagonal basis =
    concat $ concat (zipWith (zipWith smearElectrons) configs grouped)
    where
      groupedAtoms = groupBy (\a b -> (atom a) == (atom b)) basis
      grouped = map (groupBy sameShell) groupedAtoms
      configs = map (electronConfig . elem . head) groupedAtoms
      elem b = G.element (atom b)
