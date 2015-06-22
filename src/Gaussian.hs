
module Gaussian where

import Point3D

data Gaussian = 
  Gaussian { center :: Point3D
           , α :: Double
           , xIndex :: Int
           , yIndex :: Int
           , zIndex :: Int
           } deriving (Eq)

instance Show Gaussian where
    show g 
      | xIndex g > 0 = "<Px: " ++ (show $ α g) ++ ">"
      | yIndex g > 0 = "<Py: " ++ (show $ α g) ++ ">"
      | zIndex g > 0 = "<Pz: " ++ (show $ α g) ++ ">"
      | otherwise = "<S: " ++ (show $ α g) ++ ">"


instance Ord Gaussian where
    compare Gaussian {xIndex = x1, yIndex = y1, zIndex = z1}
            Gaussian {xIndex = x2, yIndex = y2, zIndex = z2} 
            | o1 > o2 = LT
            | o1 < o2 = GT
            | otherwise = EQ
            where 
              o1 = (x1 + y1 + z1) 
              o2 = (x2 + y2 + z2)
