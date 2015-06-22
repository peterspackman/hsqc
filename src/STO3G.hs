{-# LANGUAGE OverloadedStrings, DeriveGeneric #-}
module STO3G where
import Gaussian hiding (center)
import Data.Aeson
import Orbitals
import Shell
import qualified Data.ByteString.Lazy as B

sto3gBasis :: IO (Maybe BasisSet)
sto3gBasis = 
    decode <$> getJSON :: IO (Maybe BasisSet)

jsonFile :: FilePath
jsonFile = "basis/STO-3G.json"

getJSON :: IO B.ByteString
getJSON = B.readFile jsonFile

