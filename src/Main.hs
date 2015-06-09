{-# LANGUAGE OverloadedStrings #-}
import XYZParser (xyzParser)
import System.Environment
import HartreeFock
import Point3D

hcenters = [(Point3D 0.0 1.4 0.0), (Point3D 0.0 0.0 0.0)]
hehcenters = [(Point3D 0.0 1.4632 0.0), (Point3D 0.0 0.0 0.0)]

main :: IO ()
main = do
    putStrLn $ "Integrals: " ++ (show $ integrals (initSystem hcenters [1,1])) 
    putStrLn $ "\nHF (H2): " ++ (show (calculateSCF (initSystem hcenters [1,1]) 0.000001) )
    putStrLn $ "Szabo energy: " ++ (show (-1.8310))
    putStrLn $ "\nHF (HeH): " ++ (show (calculateSCF (initSystem hehcenters [2,1]) 0.000001) )
    putStrLn $ "Szabo energy: " ++ (show (-4.2275))
