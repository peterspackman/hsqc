{-# LANGUAGE MultiParamTypeClasses, BangPatterns #-}
module HartreeFock where
import Data.List (sortBy)
import Numeric.LinearAlgebra hiding (Element)
import qualified Data.Vector.Storable as V
import Data.Array.Repa hiding ((++), sum, map, zipWith, replicate, reshape, toList)
import qualified Data.Array.Repa as Repa (map, reshape, toList, transpose)
import BasisFunction
import Shell
import Element hiding (atomicNumber)
import Geometry
import Matrix (fromDiagonal, force, row, col)
import Point3D (Point3D, euclidean)
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
                     , numElectrons :: Int
                     , densityMatrix :: Array U DIM2 Double
                     , integrals :: Array U DIM4 Double
                     , coreHamiltonian :: Matrix Double
                     , overlap :: Matrix Double
                     , energies :: Energies
                     }

instance Show System where
    show System {atoms = a, energies = x} = 
      (unlines $ map show a) 
      ++ "\nTotal energy: " ++ (show $ e x)
      ++ "\nElectronic energy: " ++ (show $ ee x)
      ++ "\nNuclear repulsion energy: " ++ (show $ enn x)

-- j == Coulomb, k == exchange
gMatrix :: (Array U DIM4 Double) -> (Array U DIM2 Double) -> (Array U DIM2 Double)
gMatrix twoElectron p =
  force $ fromFunction (extent p) g
  where
    -- (ab|cd) = integrals,
    -- sum ix = sum over cd 
    -- g (a,b) = SUM over c,d p(c,d) * [(ab|cd) - 0.5 (ad|bc)]
    
    g = (\ix -> 2.0 * (sumAllS $! (p *^ (j ix)) -^ (p *^ (Repa.map (*0.5) (k ix)))))
    -- coulomb
    j ix =  slice twoElectron (Z:.(col ix):.(row ix):.All:.All)
    -- exchange
    k ix =  slice twoElectron (Z:.All:.(col ix):.All:.(row ix))



fockMatrix :: Matrix Double -> Array U DIM2 Double -> Matrix Double
fockMatrix hcore g = 
    (hcore + gmat)
    where 
    (Z:.rows:._) = extent g
    gmat = reshape rows (fromList (Repa.toList g))

scf :: System -> System
scf System { atoms = a
           , densityMatrix = p
           , integrals = t
           , coreHamiltonian = h
           , overlap = s
           , numElectrons = ne
           , energies = eold
           } =
    (System a ne newDensityMatrix t h s e)
    where
      g = gMatrix t p
      n = row (extent p)
      pmat = reshape n (fromList (Repa.toList p))
      newDensityMatrix = pMatrix ne c
      f = fockMatrix h g 
      !c = (coefficientMatrix f s)
      -- calculate energies
      electronicEnergy = (sumElements (pmat * (h + f)))
      newEnergies Energies { ek = ek, een = een, enn = enn} =
        Energies (enn + electronicEnergy) electronicEnergy enn een 0.0 ek
      e = newEnergies eold


nuclearEnergy :: [Atom] -> Double
nuclearEnergy atoms =
    0.5 * sum (energy <$> [(a,b) | a <- atoms, b <- atoms, not (a == b)])
    where
      energy (a,b) = (fromIntegral ((charge a) * (charge b))) / (distance a b)
      charge a = atomicNumber a 
      distance Atom {center = r1} Atom {center = r2} = euclidean r1 r2

sortedEigenvectors :: Matrix Double -> Matrix Double -> Matrix Double
sortedEigenvectors a b =
    (fromColumns (snd . unzip $ sortedVals))
  where
    (values, vectors) = geigSH' a b
    sortedVals = sortBy cmpFirst (zip (toList values) (toColumns vectors)) 
    cmpFirst (a1,b1) (a2,b2) = compare a1 a2


coefficientMatrix :: Matrix Double -> Matrix Double -> Matrix Double
coefficientMatrix f s =
    c
    where
      c = sortedEigenvectors f s

pMatrix :: Int -> Matrix Double -> Array U DIM2 Double
pMatrix ne cMatrix =
    force $ fromFunction (Z:.(rows cMatrix::Int):.(cols cMatrix::Int)) vals
    where
      c = (fromColumns (take (quot ne 2) (toColumns cMatrix)))
      d = c <> (trans c)
      vals = (\(Z:.i:.j) -> (d @@> (i, j)))


initSystem :: Geometry -> (Maybe BasisSet) -> System
initSystem atoms basisSet = 
  (System atoms ne p t h s energies)
  where
    basis = case basisSet of
              Just b -> concat (map (getAtomicOrbital b) atoms)
              Nothing -> undefined
    s = overlapMatrix basis
    n = length basis
    p = soad basis -- superposition of atomic densities
    ne = sum $ concat $ map (electronConfig . element) atoms -- number of electrons
    k = kineticMatrix basis
    enn = nuclearEnergy atoms
    ek = 2.0 * trace k
    een = trace v
    v = nuclearMatrix atoms basis
    h = k + v -- core hamiltonian
    energies = (Energies 0 0 enn een 0 ek )
    !t = twoElectronMatrix basis -- compute the two body integrals

-- Converge an infinite list
converge :: (a -> a -> Bool) -> [a] -> a
converge p (x:ys@(y:_))
    | p x y     = y
    | otherwise = converge p ys

calculateSCF :: System -> Double -> System
calculateSCF system eps =
    converge (\a b -> (abs (e (energies a) - (e (energies b))) < eps)) s
    where
      s = take maxIterations (iterate scf system)

soad :: Basis -> Array U DIM2 Double
soad basis =
  computeS $ fromDiagonal (soadDiagonal basis)


trace :: Matrix Double -> Double
trace m =
  V.sum (takeDiag m)
