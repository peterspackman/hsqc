{-# LANGUAGE OverloadedStrings, QuasiQuotes #-}
import XYZParser (xyzParser)
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B
import Control.Monad (when)
import System.Environment
import HartreeFock
import System.Console.Docopt
import Shell
import BasisSets
import qualified Test as T

patterns :: Docopt
patterns = [docoptFile|USAGE.txt|]

getArgOrExit = getArgOrExitWith patterns

formGeometry path = do
  file <- B.readFile $ path
  case parseOnly xyzParser file of
      Right geom -> return geom
      Left err -> fail $ "Error while parsing " ++ (show path) ++ ": " ++ err


main :: IO ()
main = do
  args <- parseArgsOrExit patterns =<< getArgs
  when (args `isPresent` (command "scf")) $ do
    path <- args `getArgOrExit` (argument "file")
    geom <- formGeometry path
    basis <- formBasis "STO-3G" geom
    let (final, steps) = (calculateSCF (initSystem geom basis)) -- convergence at 12 decimal points
    putStrLn $ "Convergence after " ++ (show steps) ++ " cycles:"
    print final

