module Geometry where

import Element hiding (atomicNumber)
import qualified Element as E
import Point3D

data Atom = Atom { center :: Point3D
                 , element :: Element
                 } deriving (Eq)

instance Show Atom where
    show a = 
      "" ++ (show (Element.symbol $ element a)) ++ ": " ++ show (center a)

type Geometry = [Atom]

bohrPerAngstrom = 1.889725989

unitConvert :: Geometry -> Geometry
unitConvert g =
    map c g
    where
      c Atom {center = c, element = e} =
        Atom (dmult bohrPerAngstrom c) e
atomicNumber :: Atom -> Int
atomicNumber Atom {element = e} = E.atomicNumber e

atomicSymbol :: Atom -> AtomicSymbol
atomicSymbol Atom {element = e} = E.symbol e
