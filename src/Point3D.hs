module Point3D where

import Data.Vector.Storable hiding ((++))
import qualified Data.Vector.Storable as V
import Foreign
import Foreign.C.Types

data Point3D = Point3D {-# UNPACK #-} !Double
                       {-# UNPACK #-} !Double
                       {-# UNPACK #-} !Double

data Axis = X | Y | Z deriving (Eq, Show)

instance Show Point3D where
    show (Point3D x y z) = "{" ++ (show x) ++ " " ++ (show y) ++ " " ++ (show z) ++"}"

instance Storable Point3D where
  sizeOf _ = sizeOf (undefined :: Double) * 3
  alignment _ = alignment (undefined :: Double)

  {-# INLINE peek #-}
  peek p = do
    a <- peekElemOff q 0
    b <- peekElemOff q 1
    c <- peekElemOff q 2
    return (Point3D a b c)
    where
      q = castPtr p

  poke p (Point3D a b c)= do
    pokeElemOff q 0 a
    pokeElemOff q 1 b
    pokeElemOff q 2 c
    where
      q = castPtr p

instance Eq Point3D where
    (==) (Point3D a b c) (Point3D a' b' c') = (a == a') && (b == b') && (c == c')

(@@) :: Point3D -> Int -> Double
(Point3D x y z) @@ i
  | i == 1 = x
  | i == 2 = y
  | i == 3 = z

(#) :: Point3D -> Axis -> Double
(Point3D x y z) # ax
  | ax == X = x
  | ax == Y = y
  | ax == Z = z
{-# INLINE (#) #-}

add :: Point3D -> Point3D -> Point3D
{-# INLINE add #-}
add (Point3D a b c) (Point3D a' b' c') = Point3D (a + a') (b + b') (c + c')

sub :: Point3D -> Point3D -> Point3D
{-# INLINE sub #-}
sub (Point3D a b c) (Point3D a' b' c') = Point3D (a - a') (b - b') (c - c')

mult :: Point3D -> Point3D -> Point3D
{-# INLINE mult #-}
mult (Point3D a b c) (Point3D a' b' c') = Point3D (a * a') (b * b') (c * c')

vsum :: Point3D -> Double
{-# INLINE vsum #-}
vsum (Point3D a b c) = a + b + c

norm :: Point3D -> Double
{-# INLINE norm #-}
norm (Point3D a b c) = sqrt (a*a + b*b + c*c)

norm2 :: Point3D -> Double
{-# INLINE norm2 #-}
norm2 (Point3D a b c) = (a*a + b*b + c*c)

euclidean :: Point3D -> Point3D -> Double
{-# INLINE euclidean #-}
euclidean a b = norm (sub a b)

euclidean2 :: Point3D -> Point3D -> Double
{-# INLINE euclidean2 #-}
euclidean2 a b = norm2 (sub a b)

dmult :: Double -> Point3D -> Point3D
{-# INLINE dmult #-}
dmult d (Point3D a b c) = Point3D (a*d) (b*d) (c*d)

ddiv :: Point3D -> Double -> Point3D
{-# INLINE ddiv #-}
ddiv (Point3D a b c) d = Point3D (a/d) (b/d) (c/d)

