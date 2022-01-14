{-# LANGUAGE CPP                    #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE GADTs                  #-}
{-# LANGUAGE RecordWildCards        #-}
module Math.Regression.Simple (
    -- * Linear regression
    linear,
    linearWithErrors,
    linearWithWeights,
    -- ** Step-by-step interface
    linearFit,
    LinRegAcc (..),
    zeroLinRegAcc,
    addLinReg,
    addLinRegW,
    -- * Quadratic regression
    quadratic,
    quadraticWithErrors,
    quadraticWithWeights,
    -- ** Step-by-step interface
    quadraticFit,
    QuadRegAcc (..),
    zeroQuadRegAcc,
    addQuadReg,
    addQuadRegW,
    quadRegAccToLin,
    -- * Auxiliary types
    V2 (..),
    V3 (..),
    Fit (..),
) where

import Control.DeepSeq (NFData (..))

import qualified Data.Foldable as F

import Math.Regression.Simple.LinAlg
import Numeric.Kahan


-- $setup
-- >>> :set -XTypeApplications
--
-- >>> import Numeric (showFFloat)
--
-- Don't show too much decimal digits
--
-- >>> showDouble x = showFFloat @Double (Just (min 5 (5 - ceiling (logBase 10 x)))) x
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
-- >>> instance (Show' a, Show' b) => Show' (a, b) where showsPrec' d (x, y) = showParen True $ showsPrec' 0 x . showString ", " . showsPrec' 0 y
-- >>> instance Show' v => Show' (Fit v) where showsPrec' d (Fit p e e') = showParen (d > 10) $ showString "Fit " . showsPrec' 11 p . showChar ' ' . showsPrec' 11 e . showChar ' ' . showsPrec' 11 e'
-- >>> instance Show' V2 where showsPrec' d (V2 x y)   = showParen (d > 10) $ showString "V2 " . showsPrec' 11 x . showChar ' ' . showsPrec' 11 y
-- >>> instance Show' V3 where showsPrec' d (V3 x y z) = showParen (d > 10) $ showString "V3 " . showsPrec' 11 x . showChar ' ' . showsPrec' 11 y . showChar ' ' . showsPrec' 11 z
--
-- Inputs:
-- >>> let input1 = [(0, 1), (1, 3), (2, 5)]
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> let input3 = [(0, 2), (1, 3), (2, 6), (3, 11)]
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
linear f = fitParams . linearWithErrors f

-- | Like 'linear' but also return parameters' standard errors.
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
-- >>> PP $ linearWithErrors id input1
-- Fit (V2 2.0000 1.00000) (V2 0.00000 0.00000) 0.00000
--
-- The second input is quite good:
--
-- >>> PP $ linearWithErrors id input2
-- Fit (V2 2.0063 0.88685) (V2 0.09550 0.23826) 0.25962
--
-- But the third input isn't so much,
-- standard error of a slope argument is 20%.
--
-- >>> let input3 = [(0, 2), (1, 3), (2, 6), (3, 11)]
-- >>> PP $ linearWithErrors id input3
-- Fit (V2 3.0000 1.00000) (V2 0.63246 1.1832) 4.0000
--
linearWithErrors :: F.Foldable f => (a -> (Double, Double)) -> f a -> Fit V2
linearWithErrors f = linearFit . linRegAcc f

-- | Weighted linear regression.
--
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> PP $ linearWithErrors id input2
-- Fit (V2 2.0063 0.88685) (V2 0.09550 0.23826) 0.25962
--
-- >>> let input2w = [(0.1, 1.2, 1), (1.3, 3.1, 1), (1.9, 4.9, 1), (3.0, 7.1, 1/4), (4.1, 9.0, 1/4)]
-- >>> PP $ linearWithWeights id input2w
-- Fit (V2 2.0060 0.86993) (V2 0.18280 0.33512) 0.22074
--
linearWithWeights :: F.Foldable f => (a -> (Double, Double, Double)) -> f a -> Fit V2
linearWithWeights f = linearFit . linRegAccW f

-- | Calculate linear fit from 'LinRegAcc'.
linearFit :: LinRegAcc -> Fit V2
linearFit (LinRegAcc (Kahan n _) (Kahan x _) (Kahan x2 _) (Kahan y _) (Kahan xy _) (Kahan y2 _)) =
    Fit params errors err
  where
    matrix@(SM22 a11 _ a22) = inv (SM22 x2 x n)
    params@(V2 a b)         = mult matrix (V2 xy y)

    errors = V2 sa sb

    -- ensure that error is always non-negative.
    -- Due rounding errors, in perfect fit situations it can be slightly negative.
    err = max 0 (y2 - a * xy - b * y)
    sa  = sqrt (a11 * err / (n - 2))
    sb  = sqrt (a22 * err / (n - 2))

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
quadratic f = fitParams . quadraticWithErrors f

-- | Like 'quadratic' but also return parameters' standard errors.
--
-- >>> PP $ quadraticWithErrors id input2
-- Fit (V3 (-0.00589) 2.0313 0.87155) (V3 0.09281 0.41070 0.37841) 0.25910
--
-- >>> PP $ quadraticWithErrors id input3
-- Fit (V3 1.00000 0.00000 2.0000) (V3 0.00000 0.00000 0.00000) 0.00000
--
quadraticWithErrors :: F.Foldable f => (a -> (Double, Double)) -> f a -> Fit V3
quadraticWithErrors f = quadraticFit . quadRegAcc f

-- | Weighted quadratic regression.
--
-- >>> let input2w = [(0.1, 1.2, 1), (1.3, 3.1, 1), (1.9, 4.9, 1), (3.0, 7.1, 1/4), (4.1, 9.0, 1/4)]
-- >>> PP $ quadraticWithWeights id input2w
-- Fit (V3 0.02524 1.9144 0.91792) (V3 0.21549 0.84212 0.70414) 0.21484
--
quadraticWithWeights :: F.Foldable f => (a -> (Double, Double, Double)) -> f a -> Fit V3
quadraticWithWeights f = quadraticFit . quadRegAccW f

-- | Calculate quadratic fit from 'QuadRegAcc'.
quadraticFit :: QuadRegAcc -> Fit V3
quadraticFit (QuadRegAcc (Kahan n _) (Kahan x _) (Kahan x2 _) (Kahan x3 _) (Kahan x4 _) (Kahan y _) (Kahan xy _) (Kahan x2y _) (Kahan y2 _)) =
    Fit params errors err
  where
    matrix@(SM33 a11 
                 _   a22
                 _   _   a33) = inv (SM33 x4 
                                          x3 x2
                                          x2  x n)

    params@(V3 a b c) = mult matrix (V3 x2y xy y)

    errors = V3 sa sb sc

    err = max 0 (y2 - a * x2y - b * xy - c * y)
    sa  = sqrt (a11 * err / (n - 3))
    sb  = sqrt (a22 * err / (n - 3))
    sc  = sqrt (a33 * err / (n - 3))

quadRegAcc :: F.Foldable f => (a -> (Double, Double)) -> f a -> QuadRegAcc
quadRegAcc f = F.foldl' (\acc a -> case f a of (x,y) -> addQuadReg acc x y) zeroQuadRegAcc

quadRegAccW :: F.Foldable f => (a -> (Double, Double, Double)) -> f a -> QuadRegAcc
quadRegAccW f = F.foldl' (\acc a -> case f a of (x,y,w) -> addQuadRegW acc x y w) zeroQuadRegAcc

-------------------------------------------------------------------------------
-- Output
-------------------------------------------------------------------------------

-- | Result of a curve fit.
data Fit v = Fit
    { fitParams :: !v
    , fitErrors :: !v
    , fitR2     :: !Double
    }
  deriving Show

instance NFData v => NFData (Fit v) where
    rnf (Fit p e _) = rnf p `seq` rnf e

-------------------------------------------------------------------------------
-- LinRegAcc
-------------------------------------------------------------------------------

-- | Linear regression accumulator.
data LinRegAcc = LinRegAcc
    { lra_w  :: {-# UNPACK #-} !Kahan  -- ^ \(\sum w_i\)
    , lra_x  :: {-# UNPACK #-} !Kahan  -- ^ \(\sum x_i \)
    , lra_x2 :: {-# UNPACK #-} !Kahan  -- ^ \(\sum x_i^2 \)
    , lra_y  :: {-# UNPACK #-} !Kahan  -- ^ \(\sum y_i \)
    , lra_xy :: {-# UNPACK #-} !Kahan  -- ^ \(\sum x_i y_i \)
    , lra_y2 :: {-# UNPACK #-} !Kahan  -- ^ \(\sum y_i^2 \)
    }

instance NFData LinRegAcc where
    rnf LinRegAcc {} = ()

-- | All-zeroes 'LinRegAcc'.
zeroLinRegAcc :: LinRegAcc
zeroLinRegAcc = LinRegAcc zeroKahan zeroKahan zeroKahan zeroKahan zeroKahan zeroKahan

-- | Add a point to linreg accumulator.
addLinReg
    :: LinRegAcc
    -> Double   -- ^ x
    -> Double   -- ^ y
    -> LinRegAcc
addLinReg LinRegAcc {..} x y = LinRegAcc
    { lra_w  = addKahan lra_w  1
    , lra_x  = addKahan lra_x  x
    , lra_x2 = addKahan lra_x2 (x * x)
    , lra_y  = addKahan lra_y  y
    , lra_xy = addKahan lra_xy (x * y)
    , lra_y2 = addKahan lra_y2 (y * y)
    }

-- | Add a weighted point to linreg accumulator.
addLinRegW
    :: LinRegAcc
    -> Double   -- ^ x
    -> Double   -- ^ y
    -> Double   -- ^ w
    -> LinRegAcc
addLinRegW LinRegAcc {..} x y w = LinRegAcc
    { lra_w  = addKahan lra_w  w
    , lra_x  = addKahan lra_x  (w * x)
    , lra_x2 = addKahan lra_x2 (w * x * x)
    , lra_y  = addKahan lra_y  (w * y)
    , lra_xy = addKahan lra_xy (w * x * y)
    , lra_y2 = addKahan lra_y2 (w * y * y)
    }

-------------------------------------------------------------------------------
-- QuadRegAcc
-------------------------------------------------------------------------------

-- | Quadratic regression accumulator.
data QuadRegAcc = QuadRegAcc
    { qra_w   :: {-# UNPACK #-} !Kahan  -- ^ \(\sum w_i\)
    , qra_x   :: {-# UNPACK #-} !Kahan  -- ^ \(\sum x_i \)
    , qra_x2  :: {-# UNPACK #-} !Kahan  -- ^ \(\sum x_i^2 \)
    , qra_x3  :: {-# UNPACK #-} !Kahan  -- ^ \(\sum x_i^3 \)
    , qra_x4  :: {-# UNPACK #-} !Kahan  -- ^ \(\sum x_i^4 \)
    , qra_y   :: {-# UNPACK #-} !Kahan  -- ^ \(\sum y_i \)
    , qra_xy  :: {-# UNPACK #-} !Kahan  -- ^ \(\sum x_i y_i \)
    , qra_x2y :: {-# UNPACK #-} !Kahan  -- ^ \(\sum x_i^2 y_i \)
    , qra_y2  :: {-# UNPACK #-} !Kahan  -- ^ \(\sum y_i^2 \)
    }

instance NFData QuadRegAcc where
    rnf QuadRegAcc {} = ()

-- | All-zeroes 'QuadRegAcc'.
zeroQuadRegAcc :: QuadRegAcc
zeroQuadRegAcc = QuadRegAcc zeroKahan zeroKahan zeroKahan zeroKahan zeroKahan zeroKahan zeroKahan zeroKahan zeroKahan

-- | Add a point to quadreg accumulator.
addQuadReg
    :: QuadRegAcc
    -> Double  -- ^ x
    -> Double  -- ^ y
    -> QuadRegAcc
addQuadReg QuadRegAcc {..} x y = QuadRegAcc
    { qra_w    = addKahan qra_w   1
    , qra_x    = addKahan qra_x   x
    , qra_x2   = addKahan qra_x2  x2
    , qra_x3   = addKahan qra_x3  (x * x2)
    , qra_x4   = addKahan qra_x4  (x2 * x2)
    , qra_y    = addKahan qra_y   y
    , qra_xy   = addKahan qra_xy  (x * y)
    , qra_x2y  = addKahan qra_x2y (x2 * y)
    , qra_y2   = addKahan qra_y2  (y * y)
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
    { qra_w    = addKahan qra_w   w
    , qra_x    = addKahan qra_x   (w * x)
    , qra_x2   = addKahan qra_x2  (w * x2)
    , qra_x3   = addKahan qra_x3  (w * x * x2)
    , qra_x4   = addKahan qra_x4  (w * x2 * x2)
    , qra_y    = addKahan qra_y   (w * y)
    , qra_xy   = addKahan qra_xy  (w * x * y)
    , qra_x2y  = addKahan qra_x2y (w * x2 * y)
    , qra_y2   = addKahan qra_y2  (w * y * y)
    }
  where
    x2 = x * x

-- | Convert 'QuadRegAcc' to 'LinRegAcc'.
--
-- Using this we can try quadratic and linear fits with a single data scan.
--
quadRegAccToLin :: QuadRegAcc -> LinRegAcc
quadRegAccToLin QuadRegAcc {..} = LinRegAcc
    { lra_w  = qra_w
    , lra_x  = qra_x
    , lra_x2 = qra_x2
    , lra_y  = qra_y
    , lra_xy = qra_xy
    , lra_y2 = qra_y2
    }
