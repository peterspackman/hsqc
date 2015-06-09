{-# LANGUAGE MultiParamTypeClasses, BangPatterns #-}
module HartreeFock where
import Control.Monad.Identity (runIdentity)
import Data.List (sortBy)
import Numeric.LinearAlgebra hiding (Element)
import Data.Array.Repa hiding ((++), sum, map, zipWith, replicate, reshape, toList)
import qualified Data.Array.Repa as Repa (map, reshape, toList)
import BasisFunction
import Element
import Matrix (row, col)
import Debug.Trace

maxIterations = 1000

data System = System { atoms :: [Atom] 
                     , densityMatrix :: Array U DIM2 Double
                     , integrals :: Array U DIM4 Double
                     , coreHamiltonian :: Matrix Double
                     , transformation :: Matrix Double
                     , totalEnergy :: Double
                     , electronicEnergy :: Double
                     }

instance Show System where
    show a = "(" ++ (show $ atoms a) ++ " -> Et: " ++ (show $ totalEnergy a)
             ++ ", Ee: " ++ (show $ electronicEnergy a)
             ++ "\ncore hamiltonian:\n" ++ (dispf 4 (coreHamiltonian a))


data Atom = Atom { center :: Double
                 , element :: Element
                 }
instance Show Atom where
    show a = "(" ++ (Element.symbol $ element a) ++ ": " ++ show (center a) ++ ")"

instance Eq Atom where
    (==) a b =
      (center a == center b) &&
      (atomicNumber $ element a) == (atomicNumber $ element b)

gMatrix :: (Array U DIM4 Double) -> (Array U DIM2 Double) -> (Array U DIM2 Double)
gMatrix twoElectron p =
    runIdentity . computeP $ fromFunction (extent p) vals
    where
      vals = (\ix -> sumAllS $ (p *^ (((coulomb ix) -^ (Repa.map (*0.5) (corr ix))))))
      coulomb ix =  slice twoElectron (Z:.((col ix)::Int):.((row ix)::Int):.All:.All)
      corr ix =  slice twoElectron (Z:.((col ix)::Int):.All:.All:.((row ix)::Int))

fockMatrix :: Matrix Double -> Array U DIM2 Double -> Matrix Double
fockMatrix hcore g = 
    hcore + gmat
    where 
    (Z:.rows:.columns) = extent g
    gmat = reshape rows (fromList (Repa.toList g))

scf :: System -> System
scf System { atoms = a
           , densityMatrix = p
           , integrals = t
           , coreHamiltonian = h
           , transformation = x
           } =
    System a newDensityMatrix t h x totalEnergy electronicEnergy 
    where
      g = gMatrix t p
      (Z:.n:._) = extent p
      !pmat = reshape n (fromList (Repa.toList p))
      !f = fockMatrix h g 
      electronicEnergy = 0.5 * sumElements (trans pmat * (h + f))
      nuc = nuclearEnergy a
      totalEnergy = nuc + electronicEnergy
      c = coefficientMatrix f x
      !newDensityMatrix = pMatrix c

nuclearEnergy :: [Atom] -> Double
nuclearEnergy atoms =
    sum (energy <$> [(a,b) | a <- atoms, b <- tail atoms, not (a == b)])
    where
      energy (a,b) = (fromIntegral ((charge a) * (charge b))) / (abs ((center a) - (center b)))
      charge a = atomicNumber $ element a

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
    runIdentity . computeP $ fromFunction (Z:.((rows cMatrix)::Int):.((cols cMatrix)::Int)) vals
    where
      vals = (\(Z:.i:.j) -> 2.0 * (cMatrix @@> (i, 0)) * (cMatrix @@> (j, 0)))

initSystem :: [Double] -> [Int] -> System
initSystem r z = 
    (System atoms p t h x 0.0 0.0)
  where
    basis = zipWith (basisFunction) r z
    o = overlapMatrix basis
    n = length basis
    p = fromListUnboxed (Z:.(n::Int):.(n::Int)) (replicate (n*n) 0 ::[Double])
    k = kineticMatrix basis
    v = nuclearMatrix (zip z r) basis
    h = k + v
    (s, u) = eigSH' o
    x = u <> diag(mapVector (** (-0.5)) s ) <> (ctrans u)
    t = twoElectronMatrix basis
    elements = map elementFromNumber z
    atoms = zipWith Atom r elements 

converge :: (a -> a -> Bool) -> [a] -> a
converge p (x:ys@(y:_))
    | p x y     = y
    | otherwise = converge p ys

calculateSCF :: System -> Double -> System
calculateSCF system eps =
    converge (\a b -> (abs ((totalEnergy a) - (totalEnergy b)) < eps)) s
    where
      s = take maxIterations (iterate scf system)
