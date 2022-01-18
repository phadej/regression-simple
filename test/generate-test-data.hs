module Main where

import Control.Monad      (unless)
import System.Environment (getArgs)
import Text.Printf (printf)

import qualified System.Random.SplitMix as SM

main :: IO ()
main = do
    args <- getArgs
    unless (null args) $ do
        let g        = SM.mkSMGen 42
        let (g1, g2) = SM.splitSMGen g

        writeFile "linear.dat" $ unlines
            [ printf "%9.05f %9.05f %9.05f" x y (sqrt var)
            | (n, x, var) <- zip3
                (normals (doubles g1))
                [ 1 .. 20 :: Double ]
                vars 
            , let y = 3 * x + 5 + n * var
            ]

        writeFile "quad.dat" $ unlines
            [ printf "%9.05f %9.05f %9.05f" x y (sqrt var)
            | (n, x, var) <- zip3
                (normals (doubles g2))
                [ 1 .. 20 :: Double ]
                vars 
            , let y = 0.1 * x * x - 3 * x + 5 + n * var
            ]
        
boxMuller :: Double -> Double -> (Double, Double)
boxMuller u v = (x, y) where
    x = k * sin (2 * pi * v)
    y = k * cos (2 * pi * u)
    k = sqrt (negate 2 * log u)

doubles :: SM.SMGen -> [Double]
doubles g = let (d, g') = SM.nextDouble g in d : doubles g'

normals :: [Double] -> [Double]
normals []       = []
normals [_]      = []
normals (u:v:ds) = let (x, y) = boxMuller u v in x : y : normals ds

vars :: [Double]
vars = 0.4 : 0.5 : vars
