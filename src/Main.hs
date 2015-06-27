{-# LANGUAGE OverloadedStrings, QuasiQuotes #-}
import XYZParser (xyzParser)
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B
import Control.Monad (when)
import System.Environment
import Geometry (Geometry, unitConvert)
import HartreeFock
import System.Console.Docopt
import System.CPUTime
import Shell
import BasisSets
import qualified Test as T

patterns :: Docopt
patterns = [docoptFile|USAGE.txt|]

getArgOrExit :: Arguments -> Option -> IO String
getArgOrExit = getArgOrExitWith patterns

formGeometry :: FilePath -> IO Geometry
formGeometry path = do
  file <- B.readFile $ path
  case parseOnly xyzParser file of
      Right geom -> return geom
      Left err -> fail $ "Error while parsing " ++ (show path) ++ ": " ++ err

basisSet :: String
basisSet = "STO-3G"

main :: IO ()
main = do
  args <- parseArgsOrExit patterns =<< getArgs
  when (args `isPresent` (command "scf")) $ do
    path <- args `getArgOrExit` (argument "file")
    geom <- formGeometry path
    let g = unitConvert geom
    basis <- formBasis basisSet g
    putStrLn $ "Using " ++ (show basisSet) ++ " basis: " ++ (show $ length basis) ++ " basis functions."
    let (final, steps) = (calculateSCF (initSystem g basis)) -- convergence at 12 decimal points
    putStrLn $ "Convergence after " ++ (show steps) ++ " cycles:"
    print final
    time <- getCPUTime
    putStrLn $ "cpu time: " ++ (show $ (fromIntegral time) * 1e-12) ++"s"

