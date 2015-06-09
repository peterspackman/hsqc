module Matrix where

import Data.Array.Repa as R hiding ((++), map) 
import qualified Data.Array.Repa (map)
import Data.Array.Repa.Unsafe as R
import Data.List (intersperse)
import Data.List.Split (chunksOf)
import Text.Printf (printf, PrintfArg)
import System.Random (getStdRandom, randomR)
import Data.Array.Repa.Algorithms.Randomish (randomishDoubleArray)
import Data.List
import qualified Data.Text as T


type Matrix2D a = R.Array R.U R.DIM2 a 
type Matrix3D a = R.Array R.U R.DIM3 a
type Matrix4D a = R.Array R.U R.DIM4 a

kroneckerDelta :: Eq a => a -> a-> Int
kroneckerDelta a b
  | a == b = 1
  | otherwise = 0

randomMatrix :: Int -> Int -> Matrix2D Double
randomMatrix rows columns =
    randomishDoubleArray (Z:. (rows::Int):. (columns::Int)) 0.0 100.0 4

random4D w x y z =
    randomishDoubleArray (Z:.(w::Int):.(x::Int):.(y::Int):.(z::Int)) 0.0 100.0 4

dispf d arr = do
    let dims = listOfShape $ extent arr
        lines = (chunksOf (2 * head dims) (formatFixed d arr))
    putStrLn $ show dims
    mapM_ putStrLn $ map concat lines

fmt :: (Text.Printf.PrintfArg a, Show a) => Int -> a -> String
fmt d a = 
    printf formatString $ a
    where
      formatString = "%."++show d++"f"

textRepresentation :: Show a => [Char] -> [a] -> String
textRepresentation sep list = intercalate sep unpackedWords 
  where words = map show list
        places = maximum (map length words) + 1
        textWords = map T.pack words
        justifiedWords = map (T.justifyRight places ' ') textWords
        unpackedWords = map T.unpack justifiedWords

format sep f m = intersperse sep (map f (toList m)) 

formatFixed d x = format " " (printf ("%."++show d++"f")) $ x

mmultP :: Monad m
       => Array U DIM2 Double
       -> Array U DIM2 Double
       -> m (Array U DIM2 Double)
mmultP arr brr = do
    trr <- transpose2P brr
    let (Z:.h1:._) = extent arr
    let (Z:._:.w2) = extent brr
    computeP
      $ fromFunction (Z:.h1:.w2)
      $ \ix -> R.sumAllS
            $ R.zipWith (*)
              (unsafeSlice arr (Any:.(row ix):.All))
              (unsafeSlice trr (Any:.(col ix):.All))
{-# NOINLINE mmultP #-}

transpose2P :: Monad m
            => Array U DIM2 Double
            -> m (Array U DIM2 Double)
transpose2P arr = 
    computeUnboxedP
    $ unsafeBackpermute new_extent swap arr
    where
      swap (Z:.i:.j) = Z:.j:.i
      new_extent = swap (extent arr)
{-# INLINE transpose2P #-}

row :: DIM2 -> Int
row (Z:.r:._) = r
{-# INLINE row #-}

col :: DIM2 -> Int
col (Z:._:.c) = c
{-# INLINE col #-}
