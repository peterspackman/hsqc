{-# LANGUAGE OverloadedStrings, QuasiQuotes #-}
import XYZParser (xyzParser)
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B
import Control.Monad (when)
import System.Environment
import HartreeFock
import System.Console.Docopt
import STO3G

patterns :: Docopt
patterns = [docoptFile|USAGE.txt|]

getArgOrExit = getArgOrExitWith patterns

water = do
  file <- B.readFile $ path
  Just b <- sto3gBasis
  case parseOnly xyzParser file of
      Right geom -> return $ concat (map (getAtomicOrbital b) geom)
  where
    path = "h2o.xyz" :: FilePath


main :: IO ()
main = do
  args <- parseArgsOrExit patterns =<< getArgs

  when (args `isPresent` (command "scf")) $ do
    path <- args `getArgOrExit` (argument "file")
    file <- B.readFile $ path
    basisSet <- sto3gBasis
    case parseOnly xyzParser file of
      Left err -> putStrLn $ "Error while parsing " ++ (show path) ++ ": " ++ err
      Right geom -> print $ (calculateSCF (initSystem geom basisSet) 0.001)

