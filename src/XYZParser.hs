{-# LANGUAGE OverloadedStrings #-}
module XYZParser where

import Element (elementFromSymbol, Element)
import Control.Applicative
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B

data Atom =
  Atom { element :: Element
       , x :: Double
       , y :: Double
       , z :: Double
         } deriving Show


type Geometry = [Atom]


skipSpaceNoNewline :: Parser ()
skipSpaceNoNewline = 
    skipWhile (\x -> isSpace x && not (inClass "\n\r" x))


-- | Parser of geometry entries
atomParser :: Parser Atom
atomParser = do
  t <- many letter_ascii
  skipSpace
  x <- double
  skipSpace
  y <- double
  skipSpace
  z <- double
  skipSpaceNoNewline
  return $ Atom (elementFromSymbol t) x y z


countEntryParser :: Parser ()
countEntryParser = do
    skipSpaceNoNewline
    decimal
    return ()


commentLineParser :: Parser ()
commentLineParser =
    skipWhile (not . inClass "\n\r")


xyzParser :: Parser Geometry
xyzParser = do
    countEntryParser <* endOfLine
    commentLineParser <* endOfLine
    many $ atomParser <* endOfLine

