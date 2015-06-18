module Orbitals where

data Axis = X | Y | Z deriving (Enum, Show, Eq)
data Orbital = S Int | P Int Axis deriving (Eq, Show)
type Triplet = (Int, Int, Int)

class BasisFunction b where
  overlap :: (Num a) => b -> b -> a

type Basis a = [a] 

instance Ord Orbital where
    compare (S _) (P _ _) = LT
    compare (P _ _) (S _) = GT
    compare (P i _) (P j _ ) = compare i j
    compare (S i) (S j) = compare i j

data PrimitiveGaussian = 
  PrimitiveGaussian { angularMomentum :: Triplet
                    , center :: (Double, Double, Double) 
                    , alpha :: Double
                    } deriving (Show, Eq)

primitiveGaussian :: PrimitiveGaussian -> (Double, Double, Double) -> Double
primitiveGaussian PrimitiveGaussian { angularMomentum = (a1,a2,a3)
                                    , center = (rx,ry,rz)
                                    , alpha = alpha
                                    } (x,y,z)
  = (x - rx)^a1 * (y - ry)^a2 * (z - rz)^a3 * exp ((-alpha) * (distance))
  where distance = (x - rx)^2 + (y - ry)^2 + (z - rz)^2
