{-# LANGUAGE OverloadedStrings #-}
import XYZParser (xyzParser)
import System.Environment
import HartreeFock


main :: IO ()
main = do
    putStrLn $ "HF (H2): " ++ (show $ electronicEnergy (calculateSCF (initSystem [0, 1.4] [1,1]) 0.000001) )
    putStrLn $ "Szabo energy: " ++ (show (-1.8310))
    putStrLn $ "HF (HeH): " ++ (show $ electronicEnergy (calculateSCF (initSystem [0, 1.4632] [2,1]) 0.000001) )
    putStrLn $ "Szabo energy: " ++ (show (-4.2275))

