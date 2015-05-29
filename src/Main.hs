{-# LANGUAGE OverloadedStrings #-}
import XYZParser (xyzParser)
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B
import System.Environment
import Data.List


main :: IO ()
main = do
    (path:_) <- getArgs
    progName <- getProgName
    file <- B.readFile $ path
    case parseOnly xyzParser file of
      Left err -> putStrLn $ "Error while parsing .xyz file: " ++ err
      Right xyz -> mapM_ print xyz
