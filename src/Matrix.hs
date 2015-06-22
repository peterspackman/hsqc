module Matrix where

import Data.Array.Repa as R hiding ((++), map) 
import qualified Data.Array.Repa (map)
import Data.Array.Repa.Unsafe as R
import Data.Array.Repa.Repr.Unboxed (Unbox)
import Data.List (intersperse)
import Data.List.Split (chunksOf)
import Text.Printf (printf, PrintfArg)
import System.Random (getStdRandom, randomR)
import Data.Array.Repa.Algorithms.Randomish (randomishDoubleArray)
import Data.List
import Control.Monad.Identity (runIdentity)
import qualified Data.Text as T
import Numeric.LinearAlgebra hiding (Element)


type Matrix2D a = R.Array R.U R.DIM2 a 
type Matrix3D a = R.Array R.U R.DIM3 a
type Matrix4D a = R.Array R.U R.DIM4 a

force :: (R.Shape sh, Unbox e) => R.Array R.D sh e -> R.Array R.U sh e
force a = runIdentity ( R.computeP a)

kroneckerDelta :: Eq a => a -> a -> Int
kroneckerDelta a b
  | a == b = 1
  | otherwise = 0

δ a b = kroneckerDelta a b
{-# INLINE δ #-}

randomMatrix :: Int -> Int -> Matrix2D Double
randomMatrix rows columns =
    randomishDoubleArray (Z:. (rows::Int):. (columns::Int)) 0.0 100.0 4

random4D w x y z =
    randomishDoubleArray (Z:.(w::Int):.(x::Int):.(y::Int):.(z::Int)) 0.0 100.0 4

row :: DIM2 -> Int
row (Z:.r:._) = r
{-# INLINE row #-}

col :: DIM2 -> Int
col (Z:._:.c) = c
{-# INLINE col #-}

fromDiagonal :: [Double] -> Array D DIM2 Double
fromDiagonal d =
    fromFunction (Z:.n:.n) vals
    where
      vals = (\ix -> if (row ix == col ix) then d !! (row ix) else 0.0)
      n = length d

--
genEigSH :: Matrix Double -> Matrix Double -> Matrix Double
genEigSH a b =
    (fromColumns (snd . unzip $ sortedVals))
  where
    (values, vectors) = geigSH' a b
    sortedVals = sortBy cmpFirst (zip (Numeric.LinearAlgebra.toList values) (toColumns vectors)) 
    cmpFirst (a1,b1) (a2,b2) = compare a1 a2


