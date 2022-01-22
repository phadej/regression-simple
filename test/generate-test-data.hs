module Main where

import Control.Monad      (unless)
import Data.List          (unfoldr, zip5)
import System.Environment (getArgs)
import Text.Printf        (printf)

import qualified System.Random.SplitMix as SM

main :: IO ()
main = do
    args <- getArgs
    unless (null args) $ do
        let (g1 : g2 : g3 : g4 : _) = unfoldr (Just . SM.splitSMGen) (SM.mkSMGen 42)

        writeFile "gnuplot/linear.dat" $ unlines
            [ printf "%9.05f %9.05f %9.05f %9.05f" (x + dx) (y + dy) devX devY
            | (x, devX, devY, px, py) <- zip5
                [ 1 .. 20 :: Double ]
                devsX
                devsY
                (normals (doubles g1))
                (normals (doubles g2))
            , let y = 3 * x + 5
            , let dx = devX * px
            , let dy = devY * py
            ]

        writeFile "gnuplot/quad.dat" $ unlines
            [ printf "%9.05f %9.05f %9.05f %9.05f" (x + dx) (y + dy) devX devY
            | (x, devX, devY, px, py) <- zip5
                [ 1 .. 20 :: Double ]
                devsX
                devsY
                (normals (doubles g3))
                (normals (doubles g4))
            , let y = 0.1 * x * x - 3 * x + 5
            , let dx = devX * px
            , let dy = devY * py
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

devsY :: [Double]
devsY = 1.2 : 1.4 : 1.6 : devsY

devsX :: [Double]
devsX = 0.2 : 0.3 : devsX
