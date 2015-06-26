{-# LANGUAGE OverloadedStrings, DeriveGeneric, BangPatterns #-}
module Shell where
import Geometry
import Orbitals
import Gaussian hiding (center)

data Shell =
    Shell { atom :: Atom
          , orbital :: Orbital
          , cgfs :: [CGF] -- contraction coefficients
          } deriving (Show, Eq)

type CGF = (Double, Gaussian)

type Basis = [Shell]


sameShell Shell {orbital = o1} Shell {orbital = o2} =
  (iFromO o1) == (iFromO o2)
  where 
    iFromO (S i) = i
    iFromO (P i _) = i
