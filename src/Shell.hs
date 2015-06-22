{-# LANGUAGE OverloadedStrings, DeriveGeneric #-}
module Shell where

import Geometry
import Orbitals
import Data.Aeson
import GHC.Generics
import Gaussian hiding (center)
import Data.Map (Map, (!))

data Shell =
    Shell { atom :: Atom
          , orbital :: Orbital
          , weights :: [Double] -- contraction coefficients
          , gaussians :: [Gaussian] -- primative gaussian functions
          } deriving (Show, Eq)


type BasisSet = Map String [ContractedGaussian] 

getAtomicOrbital :: BasisSet -> Atom -> [Shell]
getAtomicOrbital b a =
    concat shells
    where
      atomBasis = b ! (atomicSymbol a)
      shells = map shell atomBasis
      shell (ContractedGaussian i k p c) =
        case k of "p" -> [ng (P i X) p c, ng (P i Y) p c, ng (P i Z) p c]
                  "s" -> [ng (S i) p c]
      ng k p c = Shell a k c (map (gaussianFromKind (center a) k) p)

gaussianFromKind c (P _ ax) alpha 
  | ax == X = Gaussian c alpha 1 0 0
  | ax == Y = Gaussian c alpha 0 1 0
  | ax == Z = Gaussian c alpha 0 0 1
gaussianFromKind c (S _) alpha = Gaussian c alpha 0 0 0 

sameShell Shell {orbital = o1} Shell {orbital = o2} =
  (iFromO o1) == (iFromO o2)
  where 
    iFromO (S i) = i
    iFromO (P i _) = i

data ContractedGaussian = 
  ContractedGaussian { shell :: Int
                     , kind :: String
                     , primitives :: [Double]
                     , coefficients :: [Double]
  } deriving (Show, Generic)

instance FromJSON ContractedGaussian 
instance ToJSON ContractedGaussian


