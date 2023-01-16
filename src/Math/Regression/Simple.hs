{-# LANGUAGE CPP                    #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE GADTs                  #-}
{-# LANGUAGE RecordWildCards        #-}
-- |
--
-- @regression-simple@ provides (hopefully) simple regression functions.
--
-- The @'linear' :: Foldable f => (a -> (Double, Double)) -> f a -> 'V2'@
-- is the simplest one.
--
-- There are variants with weights, y-errors, and x and y-errors.
-- In addition, package includes Levenberg–Marquardt algorithm implementation
-- to fit arbitrary functions (with one, two or three parameters),
-- as long as you can give their partial derivatives as well (@ad@ package is handy for that).
--
-- For multiple independent variable ordinary least squares
-- or Levenberg-Marquard with functions with \> 3 parameter you should look elsewhere.
--
-- Package has been tested to return similar results as @fit@ functionality in @gnuplot@
-- (L-M doesn't always converge to exactly the same points in parameter space).
--
module Math.Regression.Simple (
    -- * Linear regression
    linear,
    linearFit,
    linearWithWeights,
    linearWithYerrors,
    linearWithXYerrors,
    -- ** Step-by-step interface
    linearFit',
    LinRegAcc (..),
    zeroLinRegAcc,
    addLinReg,
    addLinRegW,
    -- * Quadratic regression
    quadratic,
    quadraticFit,
    quadraticWithWeights,
    quadraticWithYerrors,
    quadraticWithXYerrors,
    -- ** Step-by-step interface
    quadraticFit',
    QuadRegAcc (..),
    zeroQuadRegAcc,
    addQuadReg,
    addQuadRegW,
    quadRegAccToLin,
    -- * Levenberg–Marquardt algorithm
    -- ** One parameter
    levenbergMarquardt1,
    levenbergMarquardt1WithWeights,
    levenbergMarquardt1WithYerrors,
    levenbergMarquardt1WithXYerrors,
    -- ** Two parameters
    levenbergMarquardt2,
    levenbergMarquardt2WithWeights,
    levenbergMarquardt2WithYerrors,
    levenbergMarquardt2WithXYerrors,
    -- ** Three parameters
    levenbergMarquardt3,
    levenbergMarquardt3WithWeights,
    levenbergMarquardt3WithYerrors,
    levenbergMarquardt3WithXYerrors,
    -- * Auxiliary types
    Fit (..),
    V2 (..),
    V3 (..),
) where

import Control.DeepSeq (NFData (..))

import qualified Data.Foldable      as F
import qualified Data.List.NonEmpty as NE

import Math.Regression.Simple.LinAlg
import Numeric.KBN

-- $setup
-- >>> :set -XDeriveFunctor -XDeriveFoldable -XDeriveTraversable
--
-- >>> import Numeric (showFFloat)
-- >>> import Data.List.NonEmpty (NonEmpty (..))
-- >>> import qualified Data.List.NonEmpty as NE
--
-- Don't show too much decimal digits
--
-- >>> let showDouble x = showFFloat (Just (min 5 (5 - ceiling (logBase 10 x)))) (x :: Double)
-- >>> showDouble 123.456 ""
-- "123.46"
--
-- >>> showDouble 1234567 ""
-- "1234567"
--
-- >>> showDouble 123.4567890123456789 ""
-- "123.46"
--
-- >>> showDouble 0.0000000000000012345 ""
-- "0.00000"
--
-- >>> newtype PP a = PP a
-- >>> class Show' a where showsPrec' :: Int -> a -> ShowS
-- >>> instance Show' a => Show (PP a) where showsPrec d (PP x) = showsPrec' d x
-- >>> instance Show' Double where showsPrec' d x = if x < 0 then showParen (d > 6) (showChar '-' . showDouble (negate x)) else showDouble x
-- >>> instance Show' Int where showsPrec' = showsPrec
-- >>> instance (Show' a, Show' b) => Show' (a, b) where showsPrec' d (x, y) = showParen True $ showsPrec' 0 x . showString ", " . showsPrec' 0 y
-- >>> instance Show' v => Show' (Fit v) where showsPrec' d (Fit p e ndf wssr) = showParen (d > 10) $ showString "Fit " . showsPrec' 11 p . showChar ' ' . showsPrec' 11 e . showChar ' ' . showsPrec' 11 ndf . showChar ' ' . showsPrec' 11 wssr
-- >>> instance Show' V2 where showsPrec' d (V2 x y)   = showParen (d > 10) $ showString "V2 " . showsPrec' 11 x . showChar ' ' . showsPrec' 11 y
-- >>> instance Show' V3 where showsPrec' d (V3 x y z) = showParen (d > 10) $ showString "V3 " . showsPrec' 11 x . showChar ' ' . showsPrec' 11 y . showChar ' ' . showsPrec' 11 z
-- >>> instance Show' a => Show' [a] where showsPrec' _ [] = id; showsPrec' _ [x] = showsPrec' 0 x; showsPrec' _ (x:xs) = showsPrec' 0 x . showChar '\n' . showsPrec' 0 xs
-- >>> instance Show' a => Show' (NonEmpty a) where showsPrec' _ (x :| []) = showsPrec' 0 x; showsPrec' _ (x:|xs) = showsPrec' 0 x . showChar '\n' . showsPrec' 0 xs
-- >>>
--
-- Inputs:
-- >>> let input1 = [(0, 1), (1, 3), (2, 5)]
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> let input3 = [(0, 2), (1, 3), (2, 6), (3, 11)]
--
-- >>> let sq z = z * z
--

-------------------------------------------------------------------------------
-- Linear
-------------------------------------------------------------------------------

-- | Linear regression.
--
-- >>> let input1 = [(0, 1), (1, 3), (2, 5)]
-- >>> PP $ linear id input1
-- V2 2.0000 1.00000
--
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> PP $ linear id input2
-- V2 2.0063 0.88685
--
linear :: F.Foldable f => (a -> (Double, Double)) -> f a -> V2
linear f = fitParams . linearFit f

-- | Like 'linear' but returns complete 'Fit'.
--
-- To get confidence intervals you should multiply the errors
-- by @quantile (studentT (n - 2)) ci'@ from @statistics@ package
-- or similar.
-- For big @n@ using value 1 gives 68% interval and using value 2 gives 95% confidence interval.
-- See https://en.wikipedia.org/wiki/Student%27s_t-distribution#Table_of_selected_values
-- (@quantile@ calculates one-sided values, you need two-sided, thus adjust @ci@ value).
--
-- The first input is perfect fit:
--
-- >>> let fit = linearFit id input1
-- >>> PP fit
-- Fit (V2 2.0000 1.00000) (V2 0.00000 0.00000) 1 0.00000
--
-- The second input is quite good:
--
-- >>> PP $ linearFit id input2
-- Fit (V2 2.0063 0.88685) (V2 0.09550 0.23826) 3 0.25962
--
-- But the third input isn't so much,
-- standard error of a slope parameter is 20%.
--
-- >>> let input3 = [(0, 2), (1, 3), (2, 6), (3, 11)]
-- >>> PP $ linearFit id input3
-- Fit (V2 3.0000 1.00000) (V2 0.63246 1.1832) 2 4.0000
--
linearFit :: F.Foldable f => (a -> (Double, Double)) -> f a -> Fit V2
linearFit f = linearFit' . linRegAcc f

-- | Weighted linear regression.
--
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> PP $ linearFit id input2
-- Fit (V2 2.0063 0.88685) (V2 0.09550 0.23826) 3 0.25962
--
-- >>> let input2w = [(0.1, 1.2, 1), (1.3, 3.1, 1), (1.9, 4.9, 1), (3.0, 7.1, 1/4), (4.1, 9.0, 1/4)]
-- >>> PP $ linearWithWeights id input2w
-- Fit (V2 2.0060 0.86993) (V2 0.12926 0.23696) 3 0.22074
--
linearWithWeights :: F.Foldable f => (a -> (Double, Double, Double)) -> f a -> Fit V2
linearWithWeights f = linearFit' . linRegAccW f

-- | Linear regression with y-errors.
--
-- >>> let input2y = [(0.1, 1.2, 0.12), (1.3, 3.1, 0.31), (1.9, 4.9, 0.49), (3.0, 7.1, 0.71), (4.1, 9.0, 1.9)]
-- >>> let fit = linearWithYerrors id input2y
-- >>> PP fit
-- Fit (V2 1.9104 0.98302) (V2 0.13006 0.10462) 3 2.0930
--
-- When we know actual y-errors, we can calculate the Q-value using @statistics@ package:
--
-- >>> import qualified Statistics.Distribution            as S
-- >>> import qualified Statistics.Distribution.ChiSquared as S
-- >>> S.cumulative (S.chiSquared (fitNDF fit)) (fitWSSR fit)
-- 0.446669639443138
--
-- or using @math-functions@
--
-- >>> import Numeric.SpecFunctions (incompleteGamma)
-- >>> incompleteGamma (fromIntegral (fitNDF fit) / 2) (fitWSSR fit / 2)
-- 0.446669639443138
--
-- It is not uncommon to deem acceptable on equal terms any models with, say, Q > 0.001.
-- If Q is too large, too near to 1 is most likely caused by overestimating
-- the y-errors.
--
--
linearWithYerrors :: F.Foldable f => (a -> (Double, Double, Double)) -> f a -> Fit V2
linearWithYerrors f = linearWithWeights f' where
    f' a = case f a of
        (x, y, dy) -> (x, y, recip (dy * dy))

-- | Iterative linear regression with x and y errors.
--
-- /Orear, J. (1982). Least squares when both variables have uncertainties. American Journal of Physics, 50(10), 912–916. doi:10.1119\/1.12972/
--
-- >>> let input2xy = [(0.1, 1.2, 0.01, 0.12), (1.3, 3.1, 0.13, 0.31), (1.9, 4.9, 0.19, 0.49), (3.0, 7.1, 0.3, 0.71), (4.1, 9.0, 0.41, 1.9)]
-- >>> let fit :| fits = linearWithXYerrors id input2xy
--
-- First fit is done using 'linearWithYerrors':
--
-- >>> PP fit
-- Fit (V2 1.9104 0.98302) (V2 0.13006 0.10462) 3 2.0930
--
-- After that the effective variance is used to refine the fit,
-- just a few iterations is often enough:
--
-- >>> PP $ take 3 fits
-- Fit (V2 1.9092 0.99251) (V2 0.12417 0.08412) 3 1.2992
-- Fit (V2 1.9092 0.99250) (V2 0.12418 0.08414) 3 1.2998
-- Fit (V2 1.9092 0.99250) (V2 0.12418 0.08414) 3 1.2998
--
linearWithXYerrors
    :: F.Foldable f
    => (a -> (Double, Double, Double, Double))  -- ^ \(x_i, y_i, \delta x_i, \delta y_i\)
    -> f a                                      -- ^ data
    -> NE.NonEmpty (Fit V2)
linearWithXYerrors f xs = iterate1 go fit0 where
    fit0   = linearWithYerrors (\a -> case f a of (x,y,_,dy) -> (x,y,dy)) xs
    go fit = linearWithWeights (\a -> case f a of (x,y,dx,dy) -> (x,y,recip $ sq (param1 * dx) + sq dy)) xs where
        V2 param1 _ = fitParams fit

-- >>> import qualified Numeric.AD.Mode.Reverse.Double as AD
-- >>> data H3 a = H3 a a a deriving (Functor, Foldable, Traversable)
-- >>> let linearF (H3 a b x) = a * x + b
-- >>> let lin' (V2 a b) (x, y, dx, dy) = case AD.grad' linearF (H3 a b x) of (f, H3 da db f') -> (y, f, V2 da db, recip $ sq (f' * dx) + sq dy)
--
-- >>> PP $ NE.last $ levenbergMarquardt2WithWeights lin' (V2 1 1) input2xy
-- Fit (V2 1.9092 0.99250) (V2 0.12418 0.08414) 3 1.2998

-- | Calculate linear fit from 'LinRegAcc'.
linearFit' :: LinRegAcc -> Fit V2
linearFit' LinRegAcc {..} = Fit params errors ndf wssr where
    matrix@(SM22 a11 _ a22) = inv (SM22 x2 x w)
    params@(V2 a b)         = mult matrix (V2 xy y)

    errors = V2 sa sb

    -- ensure that error is always non-negative.
    -- Due rounding errors, in perfect fit situations it can be slightly negative.
    wssr = max 0 (y2 - a * xy - b * y)
    ndf  = lra_n - 2
    ndf' = fromIntegral ndf :: Double

    sa   = sqrt (a11 * wssr / ndf')
    sb   = sqrt (a22 * wssr / ndf')

    w  = getKBN lra_w
    x  = getKBN lra_x
    x2 = getKBN lra_x2
    y  = getKBN lra_y
    xy = getKBN lra_xy
    y2 = getKBN lra_y2

    -- is it useful?
    -- r2 = (n * xy - x*y) / (sqrt (n * x2 - x*x) * sqrt (n * y2 - y*y))

linRegAcc :: F.Foldable f => (a -> (Double, Double)) -> f a -> LinRegAcc
linRegAcc f = F.foldl' (\acc a -> case f a of (x,y) -> addLinReg acc x y) zeroLinRegAcc

linRegAccW :: F.Foldable f => (a -> (Double, Double, Double)) -> f a -> LinRegAcc
linRegAccW f = F.foldl' (\acc a -> case f a of (x,y,w) -> addLinRegW acc x y w) zeroLinRegAcc

-------------------------------------------------------------------------------
-- Quadractic
-------------------------------------------------------------------------------

-- | Quadratic regression.
--
-- >>> let input1 = [(0, 1), (1, 3), (2, 5)]
-- >>> quadratic id input1
-- V3 0.0 2.0 1.0
--
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> PP $ quadratic id input2
-- V3 (-0.00589) 2.0313 0.87155
--
-- >>> let input3 = [(0, 2), (1, 3), (2, 6), (3, 11)]
-- >>> PP $ quadratic id input3
-- V3 1.00000 0.00000 2.0000
--
quadratic :: F.Foldable f => (a -> (Double, Double)) -> f a -> V3
quadratic f = fitParams . quadraticFit f

-- | Like 'quadratic' but returns complete 'Fit'.
--
-- >>> PP $ quadraticFit id input2
-- Fit (V3 (-0.00589) 2.0313 0.87155) (V3 0.09281 0.41070 0.37841) 2 0.25910
--
-- >>> PP $ quadraticFit id input3
-- Fit (V3 1.00000 0.00000 2.0000) (V3 0.00000 0.00000 0.00000) 1 0.00000
--
quadraticFit :: F.Foldable f => (a -> (Double, Double)) -> f a -> Fit V3
quadraticFit f = quadraticFit' . quadRegAcc f

-- | Weighted quadratic regression.
--
-- >>> let input2w = [(0.1, 1.2, 1), (1.3, 3.1, 1), (1.9, 4.9, 1), (3.0, 7.1, 1/4), (4.1, 9.0, 1/4)]
-- >>> PP $ quadraticWithWeights id input2w
-- Fit (V3 0.02524 1.9144 0.91792) (V3 0.10775 0.42106 0.35207) 2 0.21484
--
quadraticWithWeights :: F.Foldable f => (a -> (Double, Double, Double)) -> f a -> Fit V3
quadraticWithWeights f = quadraticFit' . quadRegAccW f

-- | Quadratic regression with y-errors.
--
-- >>> let input2y = [(0.1, 1.2, 0.12), (1.3, 3.1, 0.31), (1.9, 4.9, 0.49), (3.0, 7.1, 0.71), (4.1, 9.0, 0.9)]
-- >>> PP $ quadraticWithYerrors id input2y
-- Fit (V3 0.08776 1.6667 1.0228) (V3 0.10131 0.31829 0.11917) 2 1.5398
--
quadraticWithYerrors :: F.Foldable f => (a -> (Double, Double, Double)) -> f a -> Fit V3
quadraticWithYerrors f = quadraticWithWeights f' where
    f' a = case f a of
        (x, y, dy) -> (x, y, recip (dy * dy))

-- | Iterative quadratic regression with x and y errors.
--
-- /Orear, J. (1982). Least squares when both variables have uncertainties. American Journal of Physics, 50(10), 912–916. doi:10.1119\/1.12972/
--
quadraticWithXYerrors
    :: F.Foldable f
    => (a -> (Double, Double, Double, Double))  -- ^ \(x_i, y_i, \delta x_i, \delta y_i\)
    -> f a                                      -- ^ data
    -> NE.NonEmpty (Fit V3)
quadraticWithXYerrors f xs = iterate1 go fit0 where
    fit0   = quadraticWithYerrors (\a -> case f a of (x,y,_,dy) -> (x,y,dy)) xs
    go fit = quadraticWithWeights (\a -> case f a of (x,y,dx,dy) -> (x,y,recip $ sq ((2 * p1 * x + p2) * dx) + sq dy)) xs where
        V3 p1 p2 _ = fitParams fit

-- | Calculate quadratic fit from 'QuadRegAcc'.
quadraticFit' :: QuadRegAcc -> Fit V3
quadraticFit' QuadRegAcc {..} = Fit params errors ndf wssr where
    matrix@(SM33 a11
                 _   a22
                 _   _   a33) = inv (SM33 x4
                                          x3 x2
                                          x2  x w)

    params@(V3 a b c) = mult matrix (V3 x2y xy y)

    errors = V3 sa sb sc

    wssr = max 0 (y2 - a * x2y - b * xy - c * y)
    ndf  = qra_n - 3
    ndf' = fromIntegral ndf :: Double

    sa  = sqrt (a11 * wssr / ndf')
    sb  = sqrt (a22 * wssr / ndf')
    sc  = sqrt (a33 * wssr / ndf')

    w   = getKBN qra_w
    x   = getKBN qra_x
    x2  = getKBN qra_x2
    x3  = getKBN qra_x3
    x4  = getKBN qra_x4
    y   = getKBN qra_y
    xy  = getKBN qra_xy
    x2y = getKBN qra_x2y
    y2  = getKBN qra_y2

    -- is it useful?
    -- r2 = (n * xy - x*y) / (sqrt (n * x2 - x*x) * sqrt (n * y2 - y*y))

quadRegAcc :: F.Foldable f => (a -> (Double, Double)) -> f a -> QuadRegAcc
quadRegAcc f = F.foldl' (\acc a -> case f a of (x,y) -> addQuadReg acc x y) zeroQuadRegAcc

quadRegAccW :: F.Foldable f => (a -> (Double, Double, Double)) -> f a -> QuadRegAcc
quadRegAccW f = F.foldl' (\acc a -> case f a of (x,y,w) -> addQuadRegW acc x y w) zeroQuadRegAcc

-------------------------------------------------------------------------------
-- Levenberg–Marquardt 1
-------------------------------------------------------------------------------

-- | Levenberg–Marquardt for functions with one parameter.
--
-- See 'levenbergMarquardt2' for examples, this is very similar.
--
-- For example we can fit \(f = x \mapsto \beta x + 1\), its derivative is \(\partial_\beta f = x \mapsto x\).
--
-- >>> let scale a (x, y) = (y, a * x + 1, x)
-- >>> PP $ NE.last $ levenbergMarquardt1 scale 1 input2
-- Fit 1.9685 0.04735 4 0.27914
--
-- Not bad, but worse then linear fit which fits the intercept point too.
--
levenbergMarquardt1
    :: F.Foldable f
    => (Double -> a -> (Double, Double, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \partial_\beta f(\beta, x_i)\)
    -> Double                                     -- ^ initial parameter, \(\beta_0\)
    -> f a                                        -- ^ data, \(d\)
    -> NE.NonEmpty (Fit Double)                   -- ^ non-empty list of iteration results
levenbergMarquardt1 f b0 xs = loop lambda0 b0 acc0 where
    acc0 = calcAcc b0

    lambda0 = sqrt (c11 / fromIntegral n)
      where
        n   = lm1_n acc0
        c11 = getKBN $ lm1_c11 acc0

    calcAcc beta = F.foldl' (\acc p -> case f beta p of (y, g, d) -> addLM1Acc acc y g d) zeroLM1Acc xs

    loop lambda beta acc
        | lmStop lambda wssr wssr'
        = Fit beta errors ndf wssr NE.:| []

        | wssr' >= wssr
        = loop (lambda * 10) beta acc

        | otherwise
        = Fit beta errors ndf wssr `NE.cons` loop (lambda / 10) beta' acc'

      where
        lambda1 z = (1 + lambda) * z
        matrix    = inv (lambda1 c11)
        delta     = mult matrix z1

        beta' = add beta delta
        acc'  = calcAcc beta'

        a11    = inv c11
        errors = sa

        wssr  = max 0 $ getKBN $ lm1_wssr acc
        wssr' =         getKBN $ lm1_wssr acc'

        ndf   = lm1_n acc - 1
        ndf'  = fromIntegral ndf :: Double

        sa   = sqrt (a11 * wssr / ndf')

        c11 = getKBN $ lm1_c11 acc
        z1  = getKBN $ lm1_z1 acc

-- | 'levenbergMarquardt1' with weights.
levenbergMarquardt1WithWeights
    :: F.Foldable f
    => (Double -> a -> (Double, Double, Double, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \partial_\beta f(\beta, x_i), w_i\)
    -> Double                                             -- ^ initial parameter, \(\beta_0\)
    -> f a                                                -- ^ data, \(d\)
    -> NE.NonEmpty (Fit Double)                           -- ^ non-empty list of iteration results
levenbergMarquardt1WithWeights f b0 xs = loop lambda0 b0 acc0 where
    acc0 = calcAcc b0

    lambda0 = sqrt (c11 / fromIntegral n)
      where
        n   = lm1_n acc0
        c11 = getKBN $ lm1_c11 acc0

    calcAcc beta = F.foldl' (\acc p -> case f beta p of (y, g, d, w) -> addLM1AccW acc y g d w) zeroLM1Acc xs

    loop lambda beta acc
        | lmStop lambda wssr wssr'
        = Fit beta errors ndf wssr NE.:| []

        | wssr' >= wssr
        = loop (lambda * 10) beta acc

        | otherwise
        = Fit beta errors ndf wssr `NE.cons` loop (lambda / 10) beta' acc'

      where
        lambda1 z = (1 + lambda) * z
        matrix    = inv (lambda1 c11)
        delta     = mult matrix z1

        beta' = add beta delta
        acc'  = calcAcc beta'

        a11    = inv c11
        errors = sa

        wssr  = max 0 $ getKBN $ lm1_wssr acc
        wssr' =         getKBN $ lm1_wssr acc'

        ndf   = lm1_n acc - 1
        ndf'  = fromIntegral ndf :: Double

        sa   = sqrt (a11 * wssr / ndf')

        c11 = getKBN $ lm1_c11 acc
        z1  = getKBN $ lm1_z1 acc

-- | 'levenbergMarquardt1' with Y-errors.
levenbergMarquardt1WithYerrors
    :: F.Foldable f
    => (Double -> a -> (Double, Double, Double, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \partial_\beta f(\beta, x_i), \delta y_i\)
    -> Double                                             -- ^ initial parameter, \(\beta_0\)
    -> f a                                                -- ^ data, \(d\)
    -> NE.NonEmpty (Fit Double)                           -- ^ non-empty list of iteration results
levenbergMarquardt1WithYerrors f = levenbergMarquardt1WithWeights f' where
    f' beta x = case f beta x of (y, fbetax, grad, dy) -> (y, fbetax, grad, recip $ sq dy)

-- | 'levenbergMarquardt1' with XY-errors.
levenbergMarquardt1WithXYerrors
    :: F.Foldable f
    => (Double -> a -> (Double, Double, Double, Double, Double, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \partial_\beta f(\beta, x_i), \partial_x f(\beta, x_i), \delta x_i, \delta y_i\)
    -> Double                                                             -- ^ initial parameter, \(\beta_0\)
    -> f a                                                                -- ^ data, \(d\)
    -> NE.NonEmpty (Fit Double)                                           -- ^ non-empty list of iteration results
levenbergMarquardt1WithXYerrors g = levenbergMarquardt1WithWeights g' where
    g' beta x = case g beta x of (y, fbetax, grad, f', dx, dy) -> (y, fbetax, grad, recip $ sq (f' * dx) + sq dy)

data LM1Acc = LM1Acc
    { lm1_n    :: !Int
    , lm1_c11  :: !KBN
    , lm1_z1   :: !KBN
    , lm1_wssr :: !KBN
    }
  deriving Show

zeroLM1Acc :: LM1Acc
zeroLM1Acc = LM1Acc 0 zeroKBN zeroKBN zeroKBN

addLM1Acc :: LM1Acc -> Double -> Double -> Double -> LM1Acc
addLM1Acc LM1Acc {..} y f d1 = LM1Acc
    { lm1_n    = lm1_n + 1
    , lm1_c11  = addKBN lm1_c11  (d1 * d1)
    , lm1_z1   = addKBN lm1_z1   (d1 * res)
    , lm1_wssr = addKBN lm1_wssr (res * res)
    }
  where
    res = y - f

addLM1AccW :: LM1Acc -> Double -> Double -> Double -> Double -> LM1Acc
addLM1AccW LM1Acc {..} y f d1 w = LM1Acc
    { lm1_n    = lm1_n + 1
    , lm1_c11  = addKBN lm1_c11  (w * d1 * d1)
    , lm1_z1   = addKBN lm1_z1   (w * d1 * res)
    , lm1_wssr = addKBN lm1_wssr (w * res * res)
    }
  where
    res = y - f

-------------------------------------------------------------------------------
-- Levenberg–Marquardt 2
-------------------------------------------------------------------------------

-- | Levenberg–Marquardt for functions with two parameters.
--
-- You can use this sledgehammer to do a a linear fit:
--
-- >>> let lin (V2 a b) (x, y) = (y, a * x + b, V2 x 1)
--
-- We can then use 'levenbergMarquardt2' to find a fit:
--
-- >>> PP $ levenbergMarquardt2 lin (V2 1 1) input2
-- Fit (V2 1.00000 1.00000) (V2 1.0175 2.5385) 3 29.470
-- Fit (V2 1.2782 1.4831) (V2 0.57784 1.4416) 3 9.5041
-- Fit (V2 1.7254 1.4730) (V2 0.18820 0.46952) 3 1.0082
-- Fit (V2 1.9796 0.95226) (V2 0.09683 0.24157) 3 0.26687
-- Fit (V2 2.0060 0.88759) (V2 0.09550 0.23826) 3 0.25962
-- Fit (V2 2.0063 0.88685) (V2 0.09550 0.23826) 3 0.25962
--
-- This is the same result what 'linearFit' returns:
--
-- >>> PP $ linearFit id input2
-- Fit (V2 2.0063 0.88685) (V2 0.09550 0.23826) 3 0.25962
--
-- == Using AD
--
-- You can use @ad@ to calculate derivatives for you.
--
-- >>> import qualified Numeric.AD.Mode.Reverse.Double as AD
--
-- We need a ('Traversable') homogenic triple to represent the two parameters and @x@:
--
-- >>> data H3 a = H3 a a a deriving (Functor, Foldable, Traversable)
--
-- Then we define a function @ad@ can operate with:
--
-- >>> let linearF (H3 a b x) = a * x + b
--
-- which we can use to fit the curve in generic way:
--
-- >>> let lin' (V2 a b) (x, y) = case AD.grad' linearF (H3 a b x) of (f, H3 da db _f') -> (y, f, V2 da db)
-- >>> PP $ levenbergMarquardt2 lin' (V2 1 1) input2
-- Fit (V2 1.00000 1.00000) (V2 1.0175 2.5385) 3 29.470
-- Fit (V2 1.2782 1.4831) (V2 0.57784 1.4416) 3 9.5041
-- Fit (V2 1.7254 1.4730) (V2 0.18820 0.46952) 3 1.0082
-- Fit (V2 1.9796 0.95226) (V2 0.09683 0.24157) 3 0.26687
-- Fit (V2 2.0060 0.88759) (V2 0.09550 0.23826) 3 0.25962
-- Fit (V2 2.0063 0.88685) (V2 0.09550 0.23826) 3 0.25962
--
-- == Non-polynomial example
--
-- We can fit other curves too, for example an example from Wikipedia
-- https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm#Example
--
-- >>> let rateF (H3 vmax km s) = (vmax * s) / (km + s)
-- >>> let rateF' (V2 vmax km) (x, y) = case AD.grad' rateF (H3 vmax km x) of (f, H3 vmax' km' _) -> (y, f, V2 vmax' km')
-- >>> let input = zip [0.038,0.194,0.425,0.626,1.253,2.500,3.740] [0.050,0.127,0.094,0.2122,0.2729,0.2665,0.3317]
-- >>> PP $ levenbergMarquardt2 rateF' (V2 0.9 0.2) input
-- Fit (V2 0.90000 0.20000) (V2 0.43304 0.43936) 5 1.4455
-- Fit (V2 0.61786 0.36360) (V2 0.23270 0.50259) 5 0.26730
-- Fit (V2 0.39270 0.49787) (V2 0.05789 0.24170) 5 0.01237
-- Fit (V2 0.36121 0.54525) (V2 0.04835 0.23315) 5 0.00785
-- Fit (V2 0.36168 0.55530) (V2 0.04880 0.23790) 5 0.00784
-- Fit (V2 0.36182 0.55620) (V2 0.04885 0.23826) 5 0.00784
-- Fit (V2 0.36184 0.55626) (V2 0.04885 0.23829) 5 0.00784
--
-- We get the same result as in the article: 0.362 and 0.556
--
-- The algorithm terminates when a scaling parameter \(\lambda\) becomes larger than 1e20 or smaller than 1e-20, or relative WSSR change is smaller than 1e-10, or sum-of-squared-residuals candidate becomes @NaN@ (i.e. when it would start to produce garbage).
-- You may want to terminate sooner, Numerical Recipes suggest to stop when WSSR decreases by a neglible amount absolutely or fractionally.
--
levenbergMarquardt2
    :: F.Foldable f
    => (V2 -> a -> (Double, Double, V2))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \nabla_\beta f(\beta, x_i)\)
    -> V2                                 -- ^ initial parameters, \(\beta_0\)
    -> f a                                -- ^ data, \(d\)
    -> NE.NonEmpty (Fit V2)               -- ^ non-empty list of iteration results
levenbergMarquardt2 f b0 xs = loop lambda0 b0 acc0 where
    acc0 = calcAcc b0

    calcAcc beta = F.foldl' (\acc p -> case f beta p of (y, g, d) -> addLM2Acc acc y g d) zeroLM2Acc xs

    lambda0 = sqrt $ (c11 + c22) / fromIntegral n / 2
      where
        n   = lm2_n acc0
        c11 = getKBN $ lm2_c11 acc0
        c22 = getKBN $ lm2_c22 acc0

    loop lambda beta acc
        | lmStop lambda wssr wssr'
        = Fit beta errors ndf wssr NE.:| []

        | wssr' >= wssr
        = loop (lambda * 10) beta acc

        | otherwise
        = Fit beta errors ndf wssr `NE.cons` loop (lambda / 10) beta' acc'

      where
        lambda1 z = (1 + lambda) * z
        matrix    = inv (SM22 (lambda1 c11) c12 (lambda1 c22))
        delta     = mult matrix (V2 z1 z2)

        beta' = add beta delta
        acc'  = calcAcc beta'

        SM22 a11 _ a22 = inv (SM22 c11 c12 c22)
        errors = V2 sa sb

        wssr  = max 0 $ getKBN $ lm2_wssr acc
        wssr' =         getKBN $ lm2_wssr acc'

        ndf   = lm2_n acc - 2
        ndf'  = fromIntegral ndf :: Double

        sa   = sqrt (a11 * wssr / ndf')
        sb   = sqrt (a22 * wssr / ndf')

        c11 = getKBN $ lm2_c11 acc
        c12 = getKBN $ lm2_c12 acc
        c22 = getKBN $ lm2_c22 acc
        z1  = getKBN $ lm2_z1 acc
        z2  = getKBN $ lm2_z2 acc

-- | 'levenbergMarquardt2' with weights.
--
-- Because 'levenbergMarquardt2' is an iterative algorithm,
-- not only we can use it to fit curves with known y-errors ('levenbergMarquardt2WithYerrors'),
-- but also with both x and y-errors ('levenbergMarquardt2WithXYerrors').
--
levenbergMarquardt2WithWeights
    :: F.Foldable f
    => (V2 -> a -> (Double, Double, V2, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \nabla_\beta f(\beta, x_i), w_i\)
    -> V2                                         -- ^ initial parameters, \(\beta_0\)
    -> f a                                        -- ^ data, \(d\)
    -> NE.NonEmpty (Fit V2)                       -- ^ non-empty list of iteration results
levenbergMarquardt2WithWeights f b0 xs = loop lambda0 b0 acc0 where
    acc0 = calcAcc b0

    lambda0 = sqrt $ (c11 + c22) / fromIntegral n / 2
      where
        n   = lm2_n acc0
        c11 = getKBN $ lm2_c11 acc0
        c22 = getKBN $ lm2_c22 acc0

    calcAcc beta = F.foldl' (\acc p -> case f beta p of (y, g, d, w) -> addLM2AccW acc y g d w) zeroLM2Acc xs

    loop lambda beta acc
        | lmStop lambda wssr wssr'
        = Fit beta errors ndf wssr NE.:| []

        | wssr' >= wssr
        = loop (lambda * 10) beta acc

        | otherwise
        = Fit beta errors ndf wssr `NE.cons` loop (lambda / 10) beta' acc'

      where
        lambda1 z = (1 + lambda) * z
        matrix    = inv (SM22 (lambda1 c11) c12 (lambda1 c22))
        delta     = mult matrix (V2 z1 z2)

        beta' = add beta delta
        acc'  = calcAcc beta'

        SM22 a11 _ a22 = inv (SM22 c11 c12 c22)
        errors = V2 sa sb

        wssr  = max 0 $ getKBN $ lm2_wssr acc
        wssr' =         getKBN $ lm2_wssr acc'

        ndf   = lm2_n acc - 2
        ndf'  = fromIntegral ndf :: Double

        sa   = sqrt (a11 * wssr / ndf')
        sb   = sqrt (a22 * wssr / ndf')

        c11 = getKBN $ lm2_c11 acc
        c12 = getKBN $ lm2_c12 acc
        c22 = getKBN $ lm2_c22 acc
        z1  = getKBN $ lm2_z1 acc
        z2  = getKBN $ lm2_z2 acc

-- | 'levenbergMarquardt2' with Y-errors.
levenbergMarquardt2WithYerrors
    :: F.Foldable f
    => (V2 -> a -> (Double, Double, V2, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \nabla_\beta f(\beta, x_i), \delta y_i\)
    -> V2                                         -- ^ initial parameters, \(\beta_0\)
    -> f a                                        -- ^ data, \(d\)
    -> NE.NonEmpty (Fit V2)                       -- ^ non-empty list of iteration results
levenbergMarquardt2WithYerrors f = levenbergMarquardt2WithWeights f' where
    f' beta x = case f beta x of (y, fbetax, grad, dy) -> (y, fbetax, grad, recip $ sq dy)

-- | 'levenbergMarquardt2' with XY-errors.
levenbergMarquardt2WithXYerrors
    :: F.Foldable f
    => (V2 -> a -> (Double, Double, V2, Double, Double, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \nabla_\beta f(\beta, x_i), \partial_x f(\beta, x_i), \delta x_i, \delta y_i\)
    -> V2                                                         -- ^ initial parameters, \(\beta_0\)
    -> f a                                                        -- ^ data, \(d\)
    -> NE.NonEmpty (Fit V2)                                       -- ^ non-empty list of iteration results
levenbergMarquardt2WithXYerrors g = levenbergMarquardt2WithWeights g' where
    g' beta x = case g beta x of (y, fbetax, grad, f', dx, dy) -> (y, fbetax, grad, recip $ sq (f' * dx) + sq dy)

data LM2Acc = LM2Acc
    { lm2_n    :: !Int
    , lm2_c11  :: !KBN
    , lm2_c12  :: !KBN
    , lm2_c22  :: !KBN
    , lm2_z1   :: !KBN
    , lm2_z2   :: !KBN
    , lm2_wssr :: !KBN
    }
  deriving Show

zeroLM2Acc :: LM2Acc
zeroLM2Acc = LM2Acc 0 zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN

addLM2Acc :: LM2Acc -> Double -> Double -> V2 -> LM2Acc
addLM2Acc LM2Acc {..} y f (V2 d1 d2) = LM2Acc
    { lm2_n    = lm2_n + 1
    , lm2_c11  = addKBN lm2_c11  (d1 * d1)
    , lm2_c12  = addKBN lm2_c12  (d1 * d2)
    , lm2_c22  = addKBN lm2_c22  (d2 * d2)
    , lm2_z1   = addKBN lm2_z1   (d1 * res)
    , lm2_z2   = addKBN lm2_z2   (d2 * res)
    , lm2_wssr = addKBN lm2_wssr (res * res)
    }
  where
    res = y - f

addLM2AccW :: LM2Acc -> Double -> Double -> V2 -> Double -> LM2Acc
addLM2AccW LM2Acc {..} y f (V2 d1 d2) w = LM2Acc
    { lm2_n    = lm2_n + 1
    , lm2_c11  = addKBN lm2_c11  (w * d1 * d1)
    , lm2_c12  = addKBN lm2_c12  (w * d1 * d2)
    , lm2_c22  = addKBN lm2_c22  (w * d2 * d2)
    , lm2_z1   = addKBN lm2_z1   (w * d1 * res)
    , lm2_z2   = addKBN lm2_z2   (w * d2 * res)
    , lm2_wssr = addKBN lm2_wssr (w * res * res)
    }
  where
    res = y - f

-------------------------------------------------------------------------------
-- Levenberg–Marquardt 3
-------------------------------------------------------------------------------

-- | Levenberg–Marquardt for functions with three parameters.
--
-- See 'levenbergMarquardt2' for examples, this is very similar.
--
-- >>> let quad (V3 a b c) (x, y) = (y, a * x * x + b * x + c, V3 (x * x) x 1)
-- >>> PP $ NE.last $ levenbergMarquardt3 quad (V3 2 2 2) input3
-- Fit (V3 1.00000 (-0.00000) 2.0000) (V3 0.00000 0.00000 0.00000) 1 0.00000
--
-- Same as quadratic fit, just less direct:
--
-- >>> PP $ quadraticFit id input3
-- Fit (V3 1.00000 0.00000 2.0000) (V3 0.00000 0.00000 0.00000) 1 0.00000
--
levenbergMarquardt3
    :: F.Foldable f
    => (V3 -> a -> (Double, Double, V3))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \nabla_\beta f(\beta, x_i)\)
    -> V3                                 -- ^ initial parameters, \(\beta_0\)
    -> f a                                -- ^ data, \(d\)
    -> NE.NonEmpty (Fit V3)               -- ^ non-empty list of iteration results
levenbergMarquardt3 f b0 xs = loop lambda0 b0 acc0 where
    acc0 = calcAcc b0

    calcAcc beta = F.foldl' (\acc p -> case f beta p of (y, g, d) -> addLM3Acc acc y g d) zeroLM3Acc xs

    lambda0 = sqrt $ (c11 + c22 + c33) / fromIntegral n / 3
      where
        n   = lm3_n acc0
        c11 = getKBN $ lm3_c11 acc0
        c22 = getKBN $ lm3_c22 acc0
        c33 = getKBN $ lm3_c33 acc0

    loop lambda beta acc
        | lmStop lambda wssr wssr'
        = Fit beta errors ndf wssr NE.:| []

        | wssr' >= wssr
        = loop (lambda * 10) beta acc

        | otherwise
        = Fit beta errors ndf wssr `NE.cons` loop (lambda / 10) beta' acc'

      where
        lambda1 z = (1 + lambda) * z
        matrix    = inv (SM33 (lambda1 c11)
                               c12          (lambda1 c22)
                               c13          c23           (lambda1 c33))
        delta     = mult matrix (V3 z1 z2 z3)

        beta' = add beta delta
        acc'  = calcAcc beta'

        SM33 a11
             _  a22
             _  _   a33 = inv (SM33 c11
                                    c12 c22
                                    c13 c23 c33)
        errors = V3 sa sb sc

        wssr  = max 0 $ getKBN $ lm3_wssr acc
        wssr' =         getKBN $ lm3_wssr acc'

        ndf   = lm3_n acc - 3
        ndf'  = fromIntegral ndf :: Double

        sa   = sqrt (a11 * wssr / ndf')
        sb   = sqrt (a22 * wssr / ndf')
        sc   = sqrt (a33 * wssr / ndf')

        c11 = getKBN $ lm3_c11 acc
        c12 = getKBN $ lm3_c12 acc
        c13 = getKBN $ lm3_c13 acc
        c22 = getKBN $ lm3_c22 acc
        c23 = getKBN $ lm3_c23 acc
        c33 = getKBN $ lm3_c33 acc
        z1  = getKBN $ lm3_z1 acc
        z2  = getKBN $ lm3_z2 acc
        z3  = getKBN $ lm3_z3 acc

-- | 'levenbergMarquardt3' with weights.
levenbergMarquardt3WithWeights
    :: F.Foldable f
    => (V3 -> a -> (Double, Double, V3, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \nabla_\beta f(\beta, x_i), w_i\)
    -> V3                                         -- ^ initial parameters, \(\beta_0\)
    -> f a                                        -- ^ data, \(d\)
    -> NE.NonEmpty (Fit V3)                       -- ^ non-empty list of iteration results
levenbergMarquardt3WithWeights f b0 xs = loop lambda0 b0 acc0 where
    acc0 = calcAcc b0

    lambda0 = sqrt $ (c11 + c22 + c33) / fromIntegral n / 3
      where
        n   = lm3_n acc0
        c11 = getKBN $ lm3_c11 acc0
        c22 = getKBN $ lm3_c22 acc0
        c33 = getKBN $ lm3_c33 acc0

    calcAcc beta = F.foldl' (\acc p -> case f beta p of (y, g, d, w) -> addLM3AccW acc y g d w) zeroLM3Acc xs

    loop lambda beta acc
        | lmStop lambda wssr wssr'
        = Fit beta errors ndf wssr NE.:| []

        | wssr' >= wssr
        = loop (lambda * 10) beta acc

        | otherwise
        = Fit beta errors ndf wssr `NE.cons` loop (lambda / 10) beta' acc'

      where
        lambda1 z = (1 + lambda) * z
        matrix    = inv (SM33 (lambda1 c11)
                               c12          (lambda1 c22)
                               c13          c23           (lambda1 c33))
        delta     = mult matrix (V3 z1 z2 z3)

        beta' = add beta delta
        acc'  = calcAcc beta'

        SM33 a11 _ _ a22 _ a33 = inv (SM33 c11 c12 c13 c22 c23 c33)
        errors = V3 sa sb sc

        wssr  = max 0 $ getKBN $ lm3_wssr acc
        wssr' =         getKBN $ lm3_wssr acc'

        ndf   = lm3_n acc - 3
        ndf'  = fromIntegral ndf :: Double

        sa   = sqrt (a11 * wssr / ndf')
        sb   = sqrt (a22 * wssr / ndf')
        sc   = sqrt (a33 * wssr / ndf')

        c11 = getKBN $ lm3_c11 acc
        c12 = getKBN $ lm3_c12 acc
        c13 = getKBN $ lm3_c13 acc
        c22 = getKBN $ lm3_c22 acc
        c23 = getKBN $ lm3_c23 acc
        c33 = getKBN $ lm3_c33 acc
        z1  = getKBN $ lm3_z1 acc
        z2  = getKBN $ lm3_z2 acc
        z3  = getKBN $ lm3_z3 acc

-- | 'levenbergMarquardt3' with Y-errors.
levenbergMarquardt3WithYerrors
    :: F.Foldable f
    => (V3 -> a -> (Double, Double, V3, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \nabla_\beta f(\beta, x_i), \delta y_i\)
    -> V3                                         -- ^ initial parameters, \(\beta_0\)
    -> f a                                        -- ^ data, \(d\)
    -> NE.NonEmpty (Fit V3)                       -- ^ non-empty list of iteration results
levenbergMarquardt3WithYerrors f = levenbergMarquardt3WithWeights f' where
    f' beta x = case f beta x of (y, fbetax, grad, dy) -> (y, fbetax, grad, recip $ sq dy)

-- | 'levenbergMarquardt3' with XY-errors.
levenbergMarquardt3WithXYerrors
    :: F.Foldable f
    => (V3 -> a -> (Double, Double, V3, Double, Double, Double))  -- ^ \(\beta, d_i \mapsto y_i, f(\beta, x_i), \nabla_\beta f(\beta, x_i), \partial_x f(\beta, x_i), \delta x_i, \delta y_i\)
    -> V3                                                         -- ^ initial parameters, \(\beta_0\)
    -> f a                                                        -- ^ data, \(d\)
    -> NE.NonEmpty (Fit V3)                                       -- ^ non-empty list of iteration results
levenbergMarquardt3WithXYerrors g = levenbergMarquardt3WithWeights g' where
    g' beta x = case g beta x of (y, fbetax, grad, f', dx, dy) -> (y, fbetax, grad, recip $ sq (f' * dx) + sq dy)

data LM3Acc = LM3Acc
    { lm3_n    :: !Int
    , lm3_c11  :: !KBN
    , lm3_c12  :: !KBN
    , lm3_c13  :: !KBN
    , lm3_c22  :: !KBN
    , lm3_c23  :: !KBN
    , lm3_c33  :: !KBN
    , lm3_z1   :: !KBN
    , lm3_z2   :: !KBN
    , lm3_z3   :: !KBN
    , lm3_wssr :: !KBN
    }
  deriving Show

zeroLM3Acc :: LM3Acc
zeroLM3Acc = LM3Acc 0 zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN

addLM3Acc :: LM3Acc -> Double -> Double -> V3 -> LM3Acc
addLM3Acc LM3Acc {..} y f (V3 d1 d2 d3) = LM3Acc
    { lm3_n    = lm3_n + 1
    , lm3_c11  = addKBN lm3_c11  (d1 * d1)
    , lm3_c12  = addKBN lm3_c12  (d1 * d2)
    , lm3_c13  = addKBN lm3_c12  (d1 * d3)
    , lm3_c22  = addKBN lm3_c22  (d2 * d2)
    , lm3_c23  = addKBN lm3_c22  (d2 * d3)
    , lm3_c33  = addKBN lm3_c22  (d3 * d3)
    , lm3_z1   = addKBN lm3_z1   (d1 * res)
    , lm3_z2   = addKBN lm3_z2   (d2 * res)
    , lm3_z3   = addKBN lm3_z3   (d3 * res)
    , lm3_wssr = addKBN lm3_wssr (res * res)
    }
  where
    res = y - f

addLM3AccW :: LM3Acc -> Double -> Double -> V3 -> Double -> LM3Acc
addLM3AccW LM3Acc {..} y f (V3 d1 d2 d3) w = LM3Acc
    { lm3_n    = lm3_n + 1
    , lm3_c11  = addKBN lm3_c11  (w * d1 * d1)
    , lm3_c12  = addKBN lm3_c12  (w * d1 * d2)
    , lm3_c13  = addKBN lm3_c12  (w * d1 * d3)
    , lm3_c22  = addKBN lm3_c22  (w * d2 * d2)
    , lm3_c23  = addKBN lm3_c22  (w * d2 * d3)
    , lm3_c33  = addKBN lm3_c22  (w * d3 * d3)
    , lm3_z1   = addKBN lm3_z1   (w * d1 * res)
    , lm3_z2   = addKBN lm3_z2   (w * d2 * res)
    , lm3_z3   = addKBN lm3_z3   (w * d3 * res)
    , lm3_wssr = addKBN lm3_wssr (w * res * res)
    }
  where
    res = y - f

-------------------------------------------------------------------------------
-- Output
-------------------------------------------------------------------------------

-- | Result of a curve fit.
data Fit v = Fit
    { fitParams :: !v       -- ^ fit parameters
    , fitErrors :: !v       -- ^ asympotic standard errors, /assuming a good fit/
    , fitNDF    :: !Int     -- ^ number of degrees of freedom
    , fitWSSR   :: !Double  -- ^ sum of squares of residuals
    }
  deriving Show

instance (NFData v) => NFData (Fit v) where
    rnf (Fit p e _ _) = rnf p `seq` rnf e

-------------------------------------------------------------------------------
-- LinRegAcc
-------------------------------------------------------------------------------

-- | Linear regression accumulator.
data LinRegAcc = LinRegAcc
    { lra_n  :: {-# UNPACK #-} !Int  -- ^ \(n\)
    , lra_w  :: {-# UNPACK #-} !KBN  -- ^ \(\sum w_i\)
    , lra_x  :: {-# UNPACK #-} !KBN  -- ^ \(\sum x_i \)
    , lra_x2 :: {-# UNPACK #-} !KBN  -- ^ \(\sum x_i^2 \)
    , lra_y  :: {-# UNPACK #-} !KBN  -- ^ \(\sum y_i \)
    , lra_xy :: {-# UNPACK #-} !KBN  -- ^ \(\sum x_i y_i \)
    , lra_y2 :: {-# UNPACK #-} !KBN  -- ^ \(\sum y_i^2 \)
    }
  deriving Show

instance NFData LinRegAcc where
    rnf LinRegAcc {} = ()

-- | All-zeroes 'LinRegAcc'.
zeroLinRegAcc :: LinRegAcc
zeroLinRegAcc = LinRegAcc 0 zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN

-- | Add a point to linreg accumulator.
addLinReg
    :: LinRegAcc
    -> Double   -- ^ x
    -> Double   -- ^ y
    -> LinRegAcc
addLinReg LinRegAcc {..} x y = LinRegAcc
    { lra_n  = lra_n + 1
    , lra_w  = addKBN lra_w  1
    , lra_x  = addKBN lra_x  x
    , lra_x2 = addKBN lra_x2 (x * x)
    , lra_y  = addKBN lra_y  y
    , lra_xy = addKBN lra_xy (x * y)
    , lra_y2 = addKBN lra_y2 (y * y)
    }

-- | Add a weighted point to linreg accumulator.
addLinRegW
    :: LinRegAcc
    -> Double   -- ^ x
    -> Double   -- ^ y
    -> Double   -- ^ w
    -> LinRegAcc
addLinRegW LinRegAcc {..} x y w = LinRegAcc
    { lra_n  = lra_n + 1
    , lra_w  = addKBN lra_w  w
    , lra_x  = addKBN lra_x  (w * x)
    , lra_x2 = addKBN lra_x2 (w * x * x)
    , lra_y  = addKBN lra_y  (w * y)
    , lra_xy = addKBN lra_xy (w * x * y)
    , lra_y2 = addKBN lra_y2 (w * y * y)
    }

-------------------------------------------------------------------------------
-- QuadRegAcc
-------------------------------------------------------------------------------

-- | Quadratic regression accumulator.
data QuadRegAcc = QuadRegAcc
    { qra_n   :: {-# UNPACK #-} !Int  -- ^ \(n\)
    , qra_w   :: {-# UNPACK #-} !KBN  -- ^ \(\sum w_i\)
    , qra_x   :: {-# UNPACK #-} !KBN  -- ^ \(\sum x_i \)
    , qra_x2  :: {-# UNPACK #-} !KBN  -- ^ \(\sum x_i^2 \)
    , qra_x3  :: {-# UNPACK #-} !KBN  -- ^ \(\sum x_i^3 \)
    , qra_x4  :: {-# UNPACK #-} !KBN  -- ^ \(\sum x_i^4 \)
    , qra_y   :: {-# UNPACK #-} !KBN  -- ^ \(\sum y_i \)
    , qra_xy  :: {-# UNPACK #-} !KBN  -- ^ \(\sum x_i y_i \)
    , qra_x2y :: {-# UNPACK #-} !KBN  -- ^ \(\sum x_i^2 y_i \)
    , qra_y2  :: {-# UNPACK #-} !KBN  -- ^ \(\sum y_i^2 \)
    }
  deriving Show

instance NFData QuadRegAcc where
    rnf QuadRegAcc {} = ()

-- | All-zeroes 'QuadRegAcc'.
zeroQuadRegAcc :: QuadRegAcc
zeroQuadRegAcc = QuadRegAcc 0 zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN zeroKBN

-- | Add a point to quadreg accumulator.
addQuadReg
    :: QuadRegAcc
    -> Double  -- ^ x
    -> Double  -- ^ y
    -> QuadRegAcc
addQuadReg QuadRegAcc {..} x y = QuadRegAcc
    { qra_n    = qra_n + 1
    , qra_w    = addKBN qra_w   1
    , qra_x    = addKBN qra_x   x
    , qra_x2   = addKBN qra_x2  x2
    , qra_x3   = addKBN qra_x3  (x * x2)
    , qra_x4   = addKBN qra_x4  (x2 * x2)
    , qra_y    = addKBN qra_y   y
    , qra_xy   = addKBN qra_xy  (x * y)
    , qra_x2y  = addKBN qra_x2y (x2 * y)
    , qra_y2   = addKBN qra_y2  (y * y)
    }
  where
    x2 = x * x

-- | Add a weighted point to quadreg accumulator.
addQuadRegW
    :: QuadRegAcc
    -> Double  -- ^ x
    -> Double  -- ^ y
    -> Double  -- ^ w
    -> QuadRegAcc
addQuadRegW QuadRegAcc {..} x y w = QuadRegAcc
    { qra_n    = qra_n + 1
    , qra_w    = addKBN qra_w   w
    , qra_x    = addKBN qra_x   (w * x)
    , qra_x2   = addKBN qra_x2  (w * x2)
    , qra_x3   = addKBN qra_x3  (w * x * x2)
    , qra_x4   = addKBN qra_x4  (w * x2 * x2)
    , qra_y    = addKBN qra_y   (w * y)
    , qra_xy   = addKBN qra_xy  (w * x * y)
    , qra_x2y  = addKBN qra_x2y (w * x2 * y)
    , qra_y2   = addKBN qra_y2  (w * y * y)
    }
  where
    x2 = x * x

-- | Convert 'QuadRegAcc' to 'LinRegAcc'.
--
-- Using this we can try quadratic and linear fits with a single data scan.
--
quadRegAccToLin :: QuadRegAcc -> LinRegAcc
quadRegAccToLin QuadRegAcc {..} = LinRegAcc
    { lra_n  = qra_n
    , lra_w  = qra_w
    , lra_x  = qra_x
    , lra_x2 = qra_x2
    , lra_y  = qra_y
    , lra_xy = qra_xy
    , lra_y2 = qra_y2
    }

-------------------------------------------------------------------------------
-- utils
-------------------------------------------------------------------------------

sq :: Num a => a -> a
sq x = x * x
{-# INLINE sq #-}

iterate1 :: (b -> b) -> b -> NE.NonEmpty b
iterate1 g x = NE.cons x (iterate1 g (g x))

-- | Levenberg-Marquard stop condition
lmStop :: Double -> Double -> Double -> Bool
lmStop lambda wssr wssr' =
    lambda < 1e-20 || lambda > 1e20 || isNaN wssr' || relDiff < 1e-10
  where
    relDiff = abs (wssr' - wssr) / wssr
