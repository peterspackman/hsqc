module Orbitals where

data Axis = X | Y | Z deriving (Enum, Show, Eq)
data Orbital = S Int | P Int Axis deriving (Eq, Show)
type Triplet = (Int, Int, Int)

class BasisFunction b where
  overlap :: (Num a) => b -> b -> a

instance Ord Orbital where
    compare (S _) (P _ _) = LT
    compare (P _ _) (S _) = GT
    compare (P i _) (P j _ ) = compare i j
    compare (S i) (S j) = compare i j

sameKind :: Orbital -> Orbital -> Bool
sameKind (S i) (S j) = i == j
sameKind (P i _) (P j _ ) = i == j
sameKind a b = False
