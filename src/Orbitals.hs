module Orbitals where

data OrbitalType = S Int | P Int | D Int | F Int
  deriving (Show, Eq)


data PrimitiveG = PrimitiveG { angularMomentum :: (Int, Int, Int) 
                             , center :: (Double, Double, Double) 
                             , alpha :: Double
                             } deriving (Show, Eq)
primitiveGaussian :: PrimitiveG -> (Double, Double, Double) -> Double
primitiveGaussian PrimitiveG { angularMomentum = (a1,a2,a3)
                             , center = (rx,ry,rz)
                             , alpha = alpha
                             } (x,y,z)
  = (x - rx)^a1 * (y - ry)^a2 * (z - rz)^a3 * exp ((-alpha) * (distance))
  where distance = (x - rx)^2 + (y - ry)^2 + (z - rz)^2
