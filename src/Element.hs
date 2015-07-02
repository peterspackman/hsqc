module Element ( Element
               , atomicNumber
               , elements
               , elementFromNumber
               , elementFromSymbol
               , electronConfig
               , AtomicNumber
               , AtomicSymbol(..)
               , metals
               , symbol) where

import qualified Data.Map as Map  
import Data.List (sort)
data Element = Element { atomicNumber :: AtomicNumber
                       , symbol :: AtomicSymbol
                       , _name :: String
                       , mass :: Double
                       } | Unknown deriving (Show, Eq)
type AtomicNumber = Int                       
data AtomicSymbol = H  | He | Li | Be | B  | C  | N  | O  | F  | Ne | 
                    Na | Mg | Al | Si | P  | S  | Cl | Ar | K  | Ca | 
                    Sc | Ti | V  | Cr | Mn | Fe | Co | Ni | Cu | Zn | 
                    Ga | Ge | As | Se | Br | Kr | Rb | Sr | Y  | Zr | 
                    Nb | Mo | Tc | Ru | Rh | Pd | Ag | Cd | In | Sn | 
                    Sb | Te | I  | Xe | Cs | Ba | La | Ce | Pr | Nd | 
                    Pm | Sm | Eu | Gd | Tb | Dy | Ho | Er | Tm | Yb | 
                    Lu | Hf | Ta | W  | Re | Os | Ir | Pt | Au | Hg | 
                    Tl | Pb | Bi | Po | At | Rn | Fr | Ra | Ac | Th | 
                    Pa | U  | Np | Pu | Am | Cm | Bk | Cf | Es | Fm |
                    Md | No | Lr | Rf | Db | Sg | Bh | Hs | 
                    Mt | Ds | Rg | Cn | Uut | Uug | Uup | Uuh | 
                    Uus | Uuo deriving (Show, Enum, Eq, Ord)

elements :: [Element]
elements = [ Element 1 H "Hydrogen"         1.008     
           , Element 2 He "Helium"          4.002602  
           , Element 3 Li "Lithium"         6.941     
           , Element 4 Be "Beryllium"       9.012182  
           , Element 5 B "Boron"            10.811    
           , Element 6 C "Carbon"           12.011    
           , Element 7 N "Nitrogen"         14.007    
           , Element 8 O "Oxygen"           15.999    
           , Element 9 F "Fluorine"         18.9984032
           , Element 10 Ne "Neon"           20.1797   
           , Element 11 Na "Sodium"         22.989768 
           , Element 12 Mg "Magnesium"      24.305    
           , Element 13 Al "Aluminium"      26.981    
           , Element 14 Si "Silicon"        28.085    
           , Element 15 P "Phosphous"       30.973    
           , Element 16 S "Sulphur"         32.066    
           , Element 17 Cl "Chlorine"       35.4527   
           , Element 18 Ar "Argon"          39.948    
           , Element 19 K "Potassium"       39.0983   
           , Element 20 Ca "Calcium"        40.078    
           , Element 21 Sc "Scandium"       44.95591  
           , Element 22 Ti "Titanium"       47.88     
           , Element 23 V "Vanadium"        50.9415   
           , Element 24 Cr "Chromium"       51.9961   
           , Element 25 Mn "Manganese"      54.93805  
           , Element 26 Fe "Iron"           55.845    
           , Element 27 Co "Cobalt"         58.933    
           , Element 28 Ni "Nickel"         58.6934   
           , Element 29 Cu "Copper"         63.546    
           , Element 30 Zn "Zinc"           65.39     
           , Element 31 Ga "Gallium"        69.723    
           , Element 32 Ge "Germanium"      72.61     
           , Element 33 As "Arsenic"        74.92159  
           , Element 34 Se "Selenium"       78.96     
           , Element 35 Br "Bromine"        79.904    
           , Element 36 Kr "Krypton"        83.8      
           , Element 37 Rb "Rubidium"       85.4678   
           , Element 38 Sr "Strontium"      87.62     
           , Element 39 Y "Yttrium"         88.90585  
           , Element 40 Zr "Zirconium"      91.224    
           , Element 41 Nb "Niobium"        92.90638  
           , Element 42 Mo "Molybdenum"     95.94     
           , Element 43 Tc "Technetium"     98.9063   
           , Element 44 Ru "Ruthenium"      101.07    
           , Element 45 Rh "Rhodium"        102.9055  
           , Element 46 Pd "Palladium"      106.42    
           , Element 47 Ag "Silver"         107.8682  
           , Element 48 Cd "Cadmium"        112.411   
           , Element 49 In "Indium"         114.82    
           , Element 50 Sn "Tin"            118.71    
           , Element 51 Sb "Antimony"       121.75    
           , Element 52 Te "Tellurium"      127.6     
           , Element 53 I "Iodine"          126.90447 
           , Element 54 Xe "Xenon"          131.29    
           , Element 55 Cs "Caesium"        132.90543 
           , Element 56 Ba "Barium"         137.327   
           , Element 57 La "Lanthanum"      138.9055  
           , Element 58 Ce "Cerium"         140.115   
           , Element 59 Pr "Praseodymium"   140.90765 
           , Element 60 Nd "Neodymium"      144.24    
           , Element 61 Pm "Promethium"     146.9151  
           , Element 62 Sm "Samarium"       150.36    
           , Element 63 Eu "Europium"       151.965   
           , Element 64 Gd "Gadolinium"     157.25    
           , Element 65 Tb "Terbium"        158.92534 
           , Element 66 Dy "Dysprosium"     162.5     
           , Element 67 Ho "Holmium"        164.93032 
           , Element 68 Er "Erbium"         167.26    
           , Element 69 Tm "Thulium"        168.93421 
           , Element 70 Yb "Ytterbium"      173.04    
           , Element 71 Lu "Lutetium"       174.967   
           , Element 72 Hf "Hafnium"        178.49    
           , Element 73 Ta "Tantalum"       180.9479  
           , Element 74 W "Tungsten"        183.85    
           , Element 75 Re "Rhenium"        186.207   
           , Element 76 Os "Osmium"         190.2     
           , Element 77 Ir "Iridium"        192.22    
           , Element 78 Pt "Platinum"       195.08    
           , Element 79 Au "Gold"           196.96654 
           , Element 80 Hg "Mercury"        200.59    
           , Element 81 Tl "Thallium"       204.3833  
           , Element 82 Pb "Lead"           207.2     
           , Element 83 Bi "Bismuth"        208.98037 
           , Element 84 Po "Polonium"       208.9824  
           , Element 85 At "Astatine"       209.9871  
           , Element 86 Rn "Radon"          222.0176  
           , Element 87 Fr "Francium"       223.0197  
           , Element 88 Ra "Radium"         226.0254  
           , Element 89 Ac "Actinium"       227.0278  
           , Element 90 Th "Thorium"        232.0381  
           , Element 91 Pa "Protactinium"   231.0359  
           , Element 92 U "Uranium"         238.0289  
           , Element 93 Np "Neptunium"      237.0482  
           , Element 94 Pu "Plutonium"      244.0642  
           , Element 95 Am "Americium"      243.0614
           , Element 96 Cm "Curium"         247.0703
           , Element 97 Bk "Berkelium"      247.0703
           , Element 98 Cf "Californium"    251.0796
           , Element 99 Es "Einsteinium"    252.0829
           , Element 100 Fm "Fermium"       257.0951
           , Element 101 Md "Mendelevium"   258.0986
           , Element 102 No "Nobelium"      259.1009
           , Element 103 Lr "Lawrencium"    260.1053  
           , Element 104 Rf "Rutherfordium" 261.1087  
           , Element 105 Db "Dubnium"       262.1138  
           , Element 106 Sg "Seaborgium"    263.1182  
           , Element 107 Bh "Bohrium"       262.1229  
           , Element 108 Hs "Hassium"       265       
           , Element 109 Mt "Meitnerium"    266       
           , Element 110 Ds "Darmstadtium"  269       
           , Element 111 Rg "Roentgenium"   272       
           , Element 112 Cn "Copernicium"   285       
           , Element 113 Uut "Ununtrium"    284       
           , Element 114 Uug "Ununquadium"  289       
           , Element 115 Uup "Ununpentium"  288       
           , Element 116 Uuh "Ununhexium"   293       
           , Element 117 Uus "Ununseptium"  294       
           , Element 118 Uuo "Ununoctium"   294       
           ]

transitionMetals = [Sc .. Zn] ++ [Y .. Cd] ++ [La .. Hg]
alkaliMetals = tail group1
alkalineEarthMetals = group2
chalcogens = [O, S, Se, Te, Po]
group1 = [H, Li, Na, K, Rb, Cs, Fr]
group2 = [Be, Mg, Ca, Sr, Ba, Ra]
halogens = [F, Cl, Br, I, At]
lanthanides = [La .. Lu]
actinides = [Ac .. Lr]
metalloids = [B, Si, Ge, As, Sb, Te]
metals = sort $ alkaliMetals ++ alkalineEarthMetals ++ transitionMetals ++ [Al, Ga, In, Tl, Sn, Pb, Bi, Po]
         ++ lanthanides ++ [Hf .. Hg] ++ actinides ++ [Rf .. Cn]


electronConfig :: Element -> [Int]
electronConfig e = 
    case Map.lookup (atomicNumber e) configExceptions of
      Just val -> val
      _ -> filter (> 0) $ f (fillShells (atomicNumber e)) 
      where f :: [(Int, Int, Int)] -> [Int]
            f ss = [sum (g n ss) | n <- [1..m]]
              where m = length ss
                    g l = map (\(a,_,c) -> if a == l then c else 0)

valenceElectrons :: Element -> Int
valenceElectrons e = last (electronConfig e)

covalentBounds :: Element -> Int
covalentBounds e = min n (8-n)
    where n = valenceElectrons e

subshellMaxElectrons :: [Int]
subshellMaxElectrons = [2, 6, 10, 14, 18]

shellConfigGen :: Int -> [(Int, Int)]
shellConfigGen n = [(i,m-i) | m <- [2..n+n], i <- [((m+1) `div` 2)..m-1] ] 

fillShells :: Int -> [(Int, Int, Int)]
fillShells = f (shellConfigGen 5)
        where f :: [(Int, Int)] -> Int -> [(Int, Int, Int)]
              f [] _ = []
              f ((i,j):xs) m | m == 0 = []
                             | m < l = [(i, j, m)]
                             | otherwise = (i, j, l) : f xs (m - l)
                             where l = subshellMaxElectrons !! (j-1)  


configExceptions :: Map.Map Int [Int] 
configExceptions = Map.fromList [ (24, [2, 8, 13, 1])
                                , (29, [2, 8, 18, 1]) 
                                , (41, [2, 8, 18, 12, 1])
                                , (42, [2, 8, 18, 13, 1])
                                , (44, [2, 8, 18, 15, 1])
                                , (45, [2, 8, 18, 16, 1])
                                , (46, [2, 8, 18, 18])
                                , (47, [2, 8, 18, 18, 1]) 
                                , (57, [2, 8, 18, 18, 9, 2])
                                , (58, [2, 8, 18, 19, 9, 2])
                                , (64, [2, 8, 18, 25, 9, 2])
                                , (79, [2, 8, 18, 32, 18, 1])
                                , (89, [2, 8, 18, 32, 18, 9, 2])
                                , (90, [2, 8, 18, 32, 18, 10, 2])
                                , (91, [2, 8, 18, 32, 20, 9, 2])
                                , (92, [2, 8, 18, 32, 21, 9, 2])
                                , (93, [2, 8, 18, 32, 22, 9, 2]) 
                                , (96, [2, 8, 18, 32, 25, 9, 2]) ]
elementFromNumber :: Int -> Element
elementFromNumber n = f n elements
  where f :: Int -> [Element] -> Element
        f _ [] = Unknown
        f x (e:es) | atomicNumber e == x = e
                   | otherwise = f x es

elementFromSymbol :: String -> Element
elementFromSymbol s = f s elements
  where f :: String -> [Element] -> Element
        f _ [] = Unknown
        f x (e:es) | (show $ symbol e) == x = e
                   | otherwise = f x es
