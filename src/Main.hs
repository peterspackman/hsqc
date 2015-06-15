{-# LANGUAGE OverloadedStrings, QuasiQuotes #-}
import XYZParser (xyzParser)
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B
import Control.Monad (when)
import System.Environment
import HartreeFock
import System.Console.Docopt

patterns :: Docopt
patterns = [docoptFile|USAGE.txt|]

getArgOrExit = getArgOrExitWith patterns

main :: IO ()
main = do
  args <- parseArgsOrExit patterns =<< getArgs

  when (args `isPresent` (command "scf")) $ do
    path <- args `getArgOrExit` (argument "file")
    file <- B.readFile $ path
    case parseOnly xyzParser file of
      Left err -> putStrLn $ "Error while parsing " ++ (show path) ++ ": " ++ err
      Right geom -> print $ (calculateSCF (initSystem geom) 0.001)

