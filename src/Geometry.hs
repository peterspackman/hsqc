module Geometry where

import Element hiding (atomicNumber)
import qualified Element as E
import Point3D

data Atom = Atom { center :: Point3D
                 , element :: Element
                 }

instance Show Atom where
    show a = "" ++ (Element.symbol $ element a) ++ "-" ++ show (center a)

instance Eq Atom where
    (==) a b =
      (center a == center b) &&
      (atomicNumber a) == (atomicNumber b)

type Geometry = [Atom]


atomicNumber :: Atom -> Int
atomicNumber Atom {element = e} = E.atomicNumber e
