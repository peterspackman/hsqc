{-# LANGUAGE OverloadedStrings #-}
import XYZParser (xyzParser)
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B
import System.Environment
import HartreeFock
import System.Environment

main :: IO ()
main = do
    (path:_) <- getArgs
    file <- B.readFile $ path
    case parseOnly xyzParser file of
      Left err -> putStrLn $ "Error while parsing " ++ (show path) ++ ": " ++ err
      Right atoms -> print $ (calculateSCF (initSystem atoms) 0.001)

    
