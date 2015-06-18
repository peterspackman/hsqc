{-# LANGUAGE OverloadedStrings, DeriveGeneric #-}
module STO3G where
import GHC.Generics
import Gaussian hiding (center)
import Geometry
import Data.Aeson
import Data.Map (Map, (!))
import qualified Data.ByteString.Lazy as B

type BasisSet = Map String [ContractedGaussian] 

sto3gBasis :: IO (Maybe BasisSet)
sto3gBasis = 
    decode <$> getJSON :: IO (Maybe BasisSet)

getAtomicOrbital :: BasisSet -> Atom -> [STONG]
getAtomicOrbital b a =
    concat stos
    where
      atomBasis = b ! (atomicSymbol a)
      stos = map shell atomBasis
      shell (ContractedGaussian k p c) =
        case k of "p" -> [ng "px" p c, ng "py" p c, ng "pz" p c]
                  "s" -> [ng "s" p c]
      ng k p c = STONG a c (map (gaussianFromKind (center a) k) p)

gaussianFromKind c kind alpha
  | kind == "px" = Gaussian c alpha 1 0 0
  | kind == "py" = Gaussian c alpha 0 1 0
  | kind == "pz" = Gaussian c alpha 0 0 1
  | kind == "s" = Gaussian c alpha 0 0 0 

data STONG =
    STONG { atom :: Atom
          , contractionCoefficients :: [Double] -- contraction coefficients
          , gaussians :: [Gaussian] -- primative gaussian functions
          } deriving (Show, Eq)

data ContractedGaussian = 
  ContractedGaussian { kind :: String
                     , primitives :: [Double]
                     , coefficients :: [Double]
  } deriving (Show, Generic)

instance FromJSON ContractedGaussian 
instance ToJSON ContractedGaussian

jsonFile :: FilePath
jsonFile = "basis/STO-3G.json"

getJSON :: IO B.ByteString
getJSON = B.readFile jsonFile
