{-# LANGUAGE OverloadedStrings #-}
import XYZParser (xyzParser)
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B
import System.Environment
import Data.List
import HartreeFock


main :: IO ()
main = do
    putStrLn $ "HF (H2): " ++ (show (hartreeFock [0, 1.4] [1,1]))
    putStrLn $ "Szabo energy: " ++ (show (-1.8310))
