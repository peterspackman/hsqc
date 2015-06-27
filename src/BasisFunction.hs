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
import qualified Data.Array.Repa as Repa ((!), map, reshape, toList)
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
import Orbitals hiding (X, Y, Z)
import Element (atomicNumber, electronConfig)

ε = 1e-8

contraction :: Basis -> ([Gaussian] -> Double) -> Double
contraction !bs integral =
    sum $ map (f . unzip) (sequence $ map cgfs bs)
    where 
      f (!weights, !primitives) = (product weights) * (integral primitives) 
{-# INLINE contraction #-}

overlapMatrix :: Basis -> Matrix Double
overlapMatrix basis =
  reshape (length basis) $! fromList $ overlap <$> shellPairs
  where
    shellPairs = [[a,b] | a <- basis, b <- basis]
    overlap pair = contraction pair overlapList
    overlapList (g1:g2:[]) = overlapIntegral g1 g2

kineticMatrix :: Basis -> Matrix Double
kineticMatrix basis =
  reshape (length basis) $! fromList $ kinetic <$> shellPairs
  where
    shellPairs = [[a,b] | a <- basis, b <- basis]
    kinetic pair = contraction pair kineticList
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
twoElectronMatrix basis = computeS $ backpermute sh sortIndex uniqueInts
    where
      !n = V.length vbasis :: Int
      sh = (Z:.n:.n:.n:.n)
      !vbasis = V.fromList basis
      -- compute the integral at index ix
      eri ix@(Z:.i:.j:.k:.l) = 
         if ix == (sortIndex ix) then
            -- Cauchy-Schwarz inequality
            if (gab (Z:.i:.j)) * (gab (Z:.k:.l) ) > ε then
              twoElectron (vbasis ! i) (vbasis ! j) (vbasis ! k) (vbasis ! l)
            else 0.0
         else 0.0 -- 0.0 is a dummy value to fill unique array
      -- calculate all the unique (under permutation of basis functions)
      -- integrals
      uniqueInts = force $ fromFunction sh eri
      -- helper function to get twoCenter eri for ignoring tiny integrals
      gab ix = (Repa.!) twoCenters ix
      -- cache of upper bounds for cauchy-schwarz
      twoCenters = computeUnboxedS $ fromFunction (Z:.n:.n) (\(Z:.i:.j) -> 
          twoElectron (vbasis ! i) (vbasis ! j) (vbasis ! i) (vbasis !j))
      sortIndex (Z:.i:.j:.k:.l) = 
          let (a:b:c:d:[]) = concat $ sort (sort <$> [[i,j],[k,l]]) in
          (Z:.a:.b:.c:.d)
      {-# INLINE sortIndex #-}
      twoElectron !b1 !b2 !b3 !b4 = (contraction [b1,b2,b3,b4] twoElectronList)
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
{-#INLINE sortBasisFunctions #-}


-- smear the electrons across the shells,
-- filling up s orbitals first, then p orbitals.
-- Will need changes for d orbitals!
smearElectrons :: Int -> Basis -> [Double]
smearElectrons 0 bs = []
smearElectrons ne bs = 
  (replicate nshells val) ++ (smearElectrons (ne - eUsed) right)
  where
    val = (fromIntegral eUsed) / (fromIntegral (2 * nshells))
    shell = orbital $ head bs
    eUsed -- check if we can use all available electrons on this shell
      | maxElectrons shell < ne = maxElectrons shell
      | otherwise = ne
    maxElectrons (S _) = 2
    maxElectrons (P _ _) = 6
    (!left, !right) = span ((sameKind shell) . orbital) bs
    nshells = length left

soadDiagonal :: Basis -> [Double]
soadDiagonal basis =
    concat $ concat (zipWith (zipWith smearElectrons) configs grouped)
    where
      groupedAtoms = groupBy (\a b -> (atom a) == (atom b)) basis
      grouped = map (groupBy sameShell) groupedAtoms
      configs = map (electronConfig . elem . head) groupedAtoms
      elem b = G.element (atom b)
