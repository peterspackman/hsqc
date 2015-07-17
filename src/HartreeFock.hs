{-# LANGUAGE MultiParamTypeClasses, BangPatterns #-}
module HartreeFock where
import Control.Applicative ((<$>))
import Data.List (sortBy)
import Numeric.LinearAlgebra hiding (Element)
import qualified Data.Vector.Storable as V
import Data.Array.Repa hiding ((++), sum, map, zipWith, replicate, reshape, toList)
import qualified Data.Array.Repa as Repa (sumAllS, map, reshape, toList, transpose)
import BasisFunction
import Shell
import Element hiding (atomicNumber)
import Geometry
import Data.List (unfoldr)
import Matrix (Matrix2D, Matrix4D, genEigSH, fromDiagonal, force, row, col)
import Point3D (Point3D, euclidean)

maxIterations = 100

data Energies = Energies
              { ehf :: Double
              , ee :: Double
              , enn :: Double
              , een :: Double
              , ek :: Double
              }

data System = System { atoms :: Geometry
                     , nElectrons :: Int
                     , twoBodyIntegrals :: Matrix4D Double
                     , density :: Matrix2D Double
                     , hcore :: Matrix Double
                     , overlap :: Matrix Double
                     , kinetic :: Matrix Double
                     , nuclear :: Matrix Double
                     , energies :: Energies
                     , converged :: Bool
                     }

instance Show System where
    show System {atoms = a, energies = x} =
      (unlines $ map show a) ++ (show x)

instance Show Energies where
  show x =
    "Total energy: " ++ (show $ ehf x)
    ++ "\nElectronic energy: " ++ (show $ ee x)
    ++ "\nNuclear repulsion energy: " ++ (show $ enn x)
    ++ "\nKinetic energy: " ++ (show $ ek x)


-- j == Coulomb, k == exchange
twoBodyFockMatrix :: Matrix4D Double -> Matrix2D Double -> Matrix2D Double
twoBodyFockMatrix twoElectron p =
  force $ fromFunction (extent p) g
  where
    -- (ab|cd) = integrals,
    -- sum ix = sum over cd
    -- g (a,b) =  Σ(c,d) p(c,d) * [(ab|cd) - 0.5 (ad|bc)]
    g = (\ix -> 2.0 * (sumAllS $! (p *^ (j ix)) -^ (p *^ (Repa.map (*0.5) (k ix)))))
    -- coulomb
    j ix =  slice twoElectron (Z:.(col ix):.(row ix):.All:.All)
    -- exchange
    k ix =  slice twoElectron (Z:.All:.(col ix):.All:.(row ix))

fockMatrix :: Matrix Double -> Matrix2D Double -> Matrix Double
fockMatrix hcore g =
    (hcore + gmat)
    where
    (Z:.rows:._) = extent g
    gmat = reshape rows (fromList (Repa.toList g))

scf :: System -> System
scf System { atoms = a
           , density = ρ
           , twoBodyIntegrals = t
           , hcore = h
           , overlap = s
           , nElectrons = ne
           , kinetic = k
           , nuclear = v
           , energies = old
           } =
    (System a ne t d h s k v energies done)
    where
      pmat = reshape n (fromList (Repa.toList ρ ))
      g = twoBodyFockMatrix t ρ
      n = row (extent ρ)
      d = densityMatrix ne c
      f = fockMatrix h g
      !c = (coefficientMatrix f s)
      done = (abs ((ehf energies) - (ehf old))) < 1e-12
      -- calculate energies
      electronicEnergy = (sumElements (pmat * (h + f)))
      updateEnergies Energies { enn = enn } =
        Energies (enn + electronicEnergy) electronicEnergy enn
                 (nuclearAttractionEnergy v pmat) (kineticEnergy k pmat)
      energies = updateEnergies old


nuclearRepulsionEnergy :: Geometry -> Double
nuclearRepulsionEnergy atoms =
    0.5 * sum (energy <$> [(a,b) | a <- atoms, b <- atoms, not (a == b)])
    where
      energy (a,b) = (fromIntegral ((charge a) * (charge b))) / (distance a b)
      charge a = atomicNumber a
      distance Atom {center = r1} Atom {center = r2} = euclidean r1 r2


coefficientMatrix :: Matrix Double -> Matrix Double -> Matrix Double
coefficientMatrix f s = genEigSH f s

densityMatrix :: Int -> Matrix Double -> Matrix2D Double
densityMatrix nElectrons cMatrix =
    force $ fromFunction (Z:.(rows cMatrix::Int):.(cols cMatrix::Int)) vals
    where
      c = (fromColumns (take (quot nElectrons 2) (toColumns cMatrix)))
      d = c <> (trans c)
      vals = (\(Z:.i:.j) -> (d @@> (i, j)))

initSystem :: Geometry -> Basis -> System
initSystem atoms basis =
  (System atoms ne t ρ  h s k v energies False)
  where
    s = overlapMatrix basis
    n = length basis
    ρ = soadGuess basis -- superposition of atomic densities
    ne = sum $ concat $ map (electronConfig . element) atoms -- number of electrons
    k = kineticMatrix basis
    pmat = reshape n (fromList (Repa.toList ρ))
    enn = nuclearRepulsionEnergy atoms
    ek = kineticEnergy k pmat
    een = nuclearAttractionEnergy v pmat
    v = nuclearMatrix atoms basis
    h = k + v -- core hamiltonian
    energies = (Energies 0 0 enn een ek )
    !t = twoElectronMatrix basis -- compute the two body integrals

-- Converge an infinite list
converge :: (a -> a -> Bool) -> [a] -> a
converge p (x:ys@(y:_))
    | p x y     = y
    | otherwise = converge p ys

calculateSCF :: System -> (System, Int)
calculateSCF system =
    (last s, length s)
    where
      s = take maxIterations (series system)
      series = unfoldr (\x -> if (converged x) then Nothing else Just (x, scf x))


soadGuess :: Basis -> Matrix2D Double
soadGuess basis =
  computeS $ fromDiagonal (soadDiagonal basis)

kineticEnergy :: Matrix Double -> Matrix Double -> Double
kineticEnergy m ρ =
     2.0 * (trace (m <> ρ))

-- probably wrong
nuclearAttractionEnergy :: Matrix Double -> Matrix Double -> Double
nuclearAttractionEnergy m ρ =
     trace (m <> ρ)

rmsd :: Matrix2D Double -> Matrix2D Double -> Double
rmsd a b =
    Repa.sumAllS $ Repa.map (^2) (a -^ b)

trace :: Matrix Double -> Double
trace m =
  V.sum (takeDiag m)
