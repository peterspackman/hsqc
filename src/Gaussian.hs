
module Gaussian where

import Point3D

data Gaussian =
  Gaussian { center :: Point3D
           , α :: Double
           , momenta :: (Int, Int, Int)
           } deriving (Eq)


instance Show Gaussian where
    show Gaussian { α = α, momenta = (l,m,n)}
      | l > 0 = "<Px: " ++ (show α) ++ ">"
      | m > 0 = "<Py: " ++ (show α) ++ ">"
      | n > 0 = "<Pz: " ++ (show α) ++ ">"
      | otherwise = "<S: " ++ (show $ α) ++ ">"


instance Ord Gaussian where
    compare Gaussian {momenta = (x1, y1, z1) }
            Gaussian {momenta = (x2, y2, z2) }
            | o1 > o2 = LT
            | o1 < o2 = GT
            | otherwise = EQ
            where
              o1 = (x1 + y1 + z1)
              o2 = (x2 + y2 + z2)
