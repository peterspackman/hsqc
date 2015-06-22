{-# LANGUAGE OverloadedStrings, DeriveGeneric #-}
module BasisSets where
import GHC.Generics (Generic)
import Geometry
import Data.Map (Map, (!))
import Gaussian hiding (center)
import Data.Aeson
import Orbitals
import Point3D hiding (X, Y, Z)
import Shell
import qualified Data.ByteString.Lazy as B

sto3GFile :: FilePath
sto3GFile = "basis/STO-3G.json"

type BasisSet = Map String [ContractedGaussian] 

data ContractedGaussian = 
  ContractedGaussian { shell :: Int
                     , kind :: String
                     , primitives :: [Double]
                     , coefficients :: [Double]
} deriving (Show, Generic)

instance FromJSON ContractedGaussian 
instance ToJSON ContractedGaussian


getAtomicOrbitals :: BasisSet -> Atom -> Basis
getAtomicOrbitals b a =
    concat shells
    where
      atomicBasisFunctions = b ! (atomicSymbol a)
      shells = map shell atomicBasisFunctions
      -- c = center
      -- p = primitives
      -- i = shell i.e. the 2 in 2S
      shell (ContractedGaussian i k p c) =
        case k of "p" -> [ng (P i X) p c, ng (P i Y) p c, ng (P i Z) p c]
                  "s" -> [ng (S i) p c]
      ng k p c = Shell a k c (map (gaussianFromKind (center a) k) p)

gaussianFromKind :: Point3D -> Orbital -> Double -> Gaussian
gaussianFromKind c (P _ ax) alpha 
  | ax == X = Gaussian c alpha 1 0 0
  | ax == Y = Gaussian c alpha 0 1 0
  | ax == Z = Gaussian c alpha 0 0 1
gaussianFromKind c (S _) alpha = Gaussian c alpha 0 0 0 

sto3GBasis :: IO (Maybe BasisSet)
sto3GBasis = 
    decode <$> (B.readFile sto3GFile ):: IO (Maybe BasisSet)

formBasis :: String -> Geometry -> IO Basis
formBasis "STO-3G" g = do
    bset <- sto3GBasis
    case bset of
      Just b -> return $ concat (map (getAtomicOrbitals b) g)
      Nothing -> fail $ "Error reading STO-3G basis set at: " ++ sto3GFile 
