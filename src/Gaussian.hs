module Gaussian where

import Point3D

data Gaussian = 
  Gaussian { center :: Point3D
           , alpha :: Double
           , xIndex :: Int
           , yIndex :: Int
           , zIndex :: Int
           } deriving (Show, Eq)


