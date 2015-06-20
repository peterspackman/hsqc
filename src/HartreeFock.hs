{-# LANGUAGE MultiParamTypeClasses, BangPatterns #-}
module HartreeFock where
import Data.List (sortBy)
import Numeric.LinearAlgebra hiding (Element)
import Data.Array.Repa hiding ((++), sum, map, zipWith, replicate, reshape, toList)
import qualified Data.Array.Repa as Repa (map, reshape, toList)
import BasisFunction
import STO3G
import Element hiding (atomicNumber)
import Geometry
import Matrix (force, row, col)
import Point3D (Point3D, euclidean)
import Debug.Trace
maxIterations = 10

data Energies = Energies 
              { e :: Double
              , ee :: Double
              , enn :: Double
              , een :: Double
              , ep :: Double
              , ek :: Double
              } deriving Show

data System = System { atoms :: Geometry
                     , densityMatrix :: Array U DIM2 Double
                     , integrals :: Array U DIM4 Double
                     , coreHamiltonian :: Matrix Double
                     , transformation :: Matrix Double
                     , totalEnergy :: Double
                     , electronicEnergy :: Double
                     }

instance Show System where
    show a = (show $ atoms a) ++ "\nE-total: " ++ (show $ totalEnergy a)
             ++ "\nE-electronic: " ++ (show $ electronicEnergy a)
             ++ "\nE-nuc: " ++ (show $ (totalEnergy a - electronicEnergy a))

-- j == Coulomb, k == exchange
gMatrix :: (Array U DIM4 Double) -> (Array U DIM2 Double) -> (Array U DIM2 Double)
gMatrix twoElectron p =
    force $ fromFunction (extent p) g
    where
      -- (ab|cd) = integrals,
      -- sum ix = sum over cd 
      -- g (a,b) = SUM over c,d p(c,d) * [(ab|cd) - 0.5 (ad|bc)]
      g = (\ix -> 2.0 * (sumAllS $! (p  *^ (((j ix) -^ (Repa.map (*0.5) (k ix)))))))
      -- coulomb
      j ix =  slice twoElectron (Z:.(row ix):.(col ix):.All:.All)
      -- exchange
      k ix =  slice twoElectron (Z:.(row ix):.All:.(col ix):.All)




fockMatrix :: Matrix Double -> Array U DIM2 Double -> Matrix Double
fockMatrix hcore g = 
    (trace (show gmat)) (hcore + gmat)
    where 
    (Z:.rows:._) = extent g
    gmat = reshape rows (fromList (Repa.toList g))

scf :: System -> System
scf System { atoms = a
           , densityMatrix = p
           , integrals = t
           , coreHamiltonian = h
           , transformation = x
           } =
    (System a newDensityMatrix t h x totalEnergy electronicEnergy)
    where
      g = gMatrix t p
      n = row (extent p)
      pmat = reshape n (fromList (Repa.toList p))
      f = fockMatrix h g 
      electronicEnergy = 0.5 * sumElements (pmat * (h + f))
      nuc = nuclearEnergy a
      totalEnergy = nuc + electronicEnergy
      c = coefficientMatrix f x
      newDensityMatrix = pMatrix c
      energies = Energies totalEnergy electronicEnergy  nuc  0 0 0


nuclearEnergy :: [Atom] -> Double
nuclearEnergy atoms =
    0.5 * sum (energy <$> [(a,b) | a <- atoms, b <- atoms, not (a == b)])
    where
      energy (a,b) = (fromIntegral ((charge a) * (charge b))) / (distance a b)
      charge a = atomicNumber a 
      distance Atom {center = r1} Atom {center = r2} = euclidean r1 r2

sortedEigenvectors :: Field t => Matrix t -> Matrix (Complex Double)
sortedEigenvectors m =
    (fromColumns (snd . unzip $ sortedVals))
  where
    (values, vectors) = eig m
    sortedVals = sortBy cmpFirst (zip (toList values) (toColumns vectors)) 
    cmpFirst (a1,b1) (a2,b2) = compare (realPart a1) (realPart a2)


coefficientMatrix :: Matrix Double -> Matrix Double -> Matrix Double
coefficientMatrix f x =
    (x <> realC')
    where
      e = eig f'
      c' = sortedEigenvectors f'
      f' = (ctrans x) <> f <> x
      realC' = mapMatrix realPart (c')

pMatrix :: Matrix Double -> Array U DIM2 Double
pMatrix cMatrix =
    force $ fromFunction (Z:.((rows cMatrix)::Int):.((cols cMatrix)::Int)) vals
    where
      vals = (\(Z:.i:.j) -> 2.0 * (cMatrix @@> (i, 0)) * (cMatrix @@> (j, 0)))


initSystem :: Geometry -> (Maybe BasisSet) -> System
initSystem atoms basisSet = 
  (System atoms p t h x 0.0 0.0)
  where
    basis = case basisSet of
              Just b -> concat (map (getAtomicOrbital b) atoms)
              Nothing -> undefined
    o = overlapMatrix basis
    n = length basis
    p = soad basis
    k = kineticMatrix basis
    v = nuclearMatrix atoms basis
    h = k + v
    (s, u) = eigSH' o
    x = u <> diag(mapVector (** (-0.5)) s ) <> (ctrans u)
    !t = twoElectronMatrix basis

converge :: (a -> a -> Bool) -> [a] -> a
converge p (x:ys@(y:_))
    | p x y     = y
    | otherwise = converge p ys

calculateSCF :: System -> Double -> System
calculateSCF system eps =
    converge (\a b -> (abs ((totalEnergy a) - (totalEnergy b)) < eps)) s
    where
      s = take maxIterations (iterate scf system)

sameAtom :: STONG -> STONG -> Bool
sameAtom a b = (atom a == atom b)

soad basis =
    force $ fromFunction (Z:.((n)::Int):.((n)::Int)) vals
  where
    vals = (\(Z:.i:.j) -> if i == j then diag !! i else 0.0)
    diag = [1,1,2/3,2/3,2/3,0.5,0.5::Double]
    n = length basis

-- e.g. 3 basis functions 2 electrons -> density = 2/3
soadValue nbasis ne =
    ne/nbasis

