module Geometry where

import Element hiding (atomicNumber)
import qualified Element as E
import Point3D

data Atom = Atom { center :: Point3D
                 , element :: Element
                 } deriving (Eq)

instance Show Atom where
    show a = 
      "" ++ (Element.symbol $ element a) ++ ": " ++ show (center a)

type Geometry = [Atom]


atomicNumber :: Atom -> Int
atomicNumber Atom {element = e} = E.atomicNumber e

atomicSymbol :: Atom -> String
atomicSymbol Atom {element = e} = E.symbol e
