import Numeric.LinearAlgebra
import Numeric.GSL

main :: IO ()
main = do
    let w = 4 |> [2,0,-3,0::Double]
    let m = (3><4) [1..] :: Matrix Double
    putStrLn (show (w <.> w))
    putStrLn (show (m <> w))
    putStrLn (show (sqrt . eigenvalues $ m <> trans m))
