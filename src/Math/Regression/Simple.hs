{-# LANGUAGE CPP                    #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE GADTs                  #-}
module Math.Regression.Simple (
    -- * Regressions
    linear,
    linearWithErrors,
    quadratic,
    quadraticWithErrors,
    quadraticAndLinear,
    quadraticAndLinearWithErrors,
    -- * Operations
    Add (..),
    Eye (..),
    Mult (..),
    Det (..),
    Inv (..),
    -- * Zeros
    zerosLin,
    zerosQuad,
    optimaQuad,
    -- * Two dimensions
    V2 (..),
    M22 (..),
    -- * Three dimensions
    V3 (..),
    M33 (..),
    -- * Auxiliary classes
    Foldable' (..),
    IsDoublePair (..),
    ) where

import Data.Complex (Complex (..))

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

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
-- >>> instance Show' V2 where showsPrec' d (V2 x y)   = showParen (d > 10) $ showString "V2 " . showsPrec' 11 x . showChar ' ' . showsPrec' 11 y
-- >>> instance Show' V3 where showsPrec' d (V3 x y z) = showParen (d > 10) $ showString "V3 " . showsPrec' 11 x . showChar ' ' . showsPrec' 11 y . showChar ' ' . showsPrec' 11 z
--
-- Inputs:
-- >>> let input1 = [(0, 1), (1, 3), (2, 5)]
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> let input3 = [(0, 2), (1, 3), (2, 6), (3, 11)]
--

-------------------------------------------------------------------------------
-- Classes
-------------------------------------------------------------------------------

-- | Addition
class Add a where
    zero :: a
    add  :: a -> a -> a

-- | Identity
class Eye a where
    eye :: a

-- | Multiplication of different things.
class Eye a => Mult a b c | a b -> c where
    mult :: a -> b -> c

-- | Determinant
class Eye a => Det a where
    det :: a -> Double

-- | Inverse
class Det a => Inv a where
    inv :: a -> a

infixl 6 `add`
infixl 7 `mult`

instance Eye Double where
    eye = 1

instance Add Double where
    zero = 0
    add = (+)

instance Det Double where
    det = id

instance Inv Double where
    inv = recip

-------------------------------------------------------------------------------
-- Zeros
-------------------------------------------------------------------------------

-- | Solve linear equation.
--
-- >>> zerosLin (V2 1 2)
-- -2.0
--
zerosLin :: V2 -> Double
zerosLin (V2 a b) = negate (b / a)

-- | Solve quadratic equation.
--
-- >>> zerosQuad (V3 2 0 (-1))
-- Right (-0.7071067811865476,0.7071067811865476)
--
-- >>> zerosQuad (V3 2 0 1)
-- Left ((-0.0) :+ (-0.7071067811865476),(-0.0) :+ 0.7071067811865476)
--
-- Double root is not treated separately:
--
-- >>> zerosQuad (V3 1 0 0)
-- Right (-0.0,0.0)
--
-- >>> zerosQuad (V3 1 (-2) 1)
-- Right (1.0,1.0)
--
zerosQuad :: V3 -> Either (Complex Double, Complex Double) (Double, Double)
zerosQuad (V3 a b c)
    | delta < 0 = Left ((-b/da) :+ (-sqrtNDelta/da), (-b/da) :+ (sqrtNDelta/da))
    | otherwise = Right ((- b - sqrtDelta) / da, (-b + sqrtDelta) / da)
  where
    delta = b*b - 4 * a * c
    sqrtDelta = sqrt delta
    sqrtNDelta = sqrt (- delta)
    da = 2 * a

-- | Find an optima point.
--
-- >>> optimaQuad (V3 1 (-2) 0)
-- 1.0
--
-- compare to
--
-- >>> zerosQuad (V3 1 (-2) 0)
-- Right (0.0,2.0)
--
optimaQuad :: V3 -> Double
optimaQuad (V3 a b _) = zerosLin (V2 (2 * a) b)

-------------------------------------------------------------------------------
-- 2 dimensions
-------------------------------------------------------------------------------

-- | 2d vector. Strict pair of 'Double's.
--
-- Also used to represent linear polynomial: @V2 a b@  \(= a x + b\).
--
data V2 = V2 !Double !Double
  deriving (Eq, Show)

instance Add V2 where
    zero = V2 0 0
    add (V2 x y) (V2 x' y') = V2 (x + x') (y + y')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Mult Double V2 V2 where
    mult k (V2 x y) = V2 (k * x) (k * y)
    {-# INLINE mult #-}

-- | 2×2 matrix.
data M22 = M22 !Double !Double !Double !Double
  deriving (Eq, Show)

instance Add M22 where
    zero = M22 0 0 0 0
    add (M22 a b c d) (M22 a' b' c' d') = M22 (a + a') (b + b') (c + c') (d + d')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Eye M22 where
    eye = M22 1 0 0 1
    {-# INLINE eye #-}

instance Det M22 where det = det2
instance Inv M22 where inv = inv2

instance Mult Double M22 M22 where
    mult k (M22 a b c d) = M22 (k * a) (k * b) (k * c) (k * d)
    {-# INLINE mult #-}

instance Mult M22 V2 V2 where
    mult (M22 a b c d) (V2 u v) = V2 (a * u + b * v) (c * u + d * v)
    {-# INLINE mult #-}

-- | >>> M22 1 2 3 4 `mult` eye @M22
-- M22 1.0 2.0 3.0 4.0
instance Mult M22 M22 M22 where
    mult (M22 a b c d) (M22 x y z w) = M22
        (a * x + b * z) (a * y + b * w)
        (c * x + d * z) (c * y + d * w)
    {-# INLINE mult #-}

det2 :: M22 -> Double
det2 (M22 a b c d) = a * d - b * c
{-# INLINE det2 #-}

inv2 :: M22 -> M22
inv2 m@(M22 a b c d) = M22
    (  d / det) (- b / det)
    (- c / det) (  a / det)
  where
    det = det2 m
{-# INLINE inv2 #-}

-------------------------------------------------------------------------------
-- 3 dimensions
-------------------------------------------------------------------------------

-- | 3d vector. Strict triple of 'Double's.
--
-- Also used to represent quadratic polynomial: @V3 a b c@  \(= a x^2 + b x + c\).
data V3 = V3 !Double !Double !Double
  deriving (Eq, Show)

instance Add V3 where
    zero = V3 0 0 0
    add (V3 x y z) (V3 x' y' z') = V3 (x + x') (y + y') (z + z')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Mult Double V3 V3 where
    mult k (V3 x y z) = V3 (k * x) (k * y) (k * z)
    {-# INLINE mult #-}

-- | 3×3 matrix.
data M33 = M33
    !Double !Double !Double
    !Double !Double !Double
    !Double !Double !Double
  deriving (Eq, Show)

instance Add M33 where
    zero = M33 0 0 0 0 0 0 0 0 0

    add (M33 a b c d e f g h i) (M33 a' b' c' d' e' f' g' h' i') = M33
        (a + a') (b + b') (c + c')
        (d + d') (e + e') (f + f')
        (g + g') (h + h') (i + i')
    {-# INLINE zero #-}
    {-# INLINE add #-}

instance Eye M33 where
    eye = M33 1 0 0 0 1 0 0 0 1
    {-# INLINE eye #-}

instance Det M33 where det = det3
instance Inv M33 where inv = inv3

instance Mult Double M33 M33 where
    mult k (M33 a b c d e f g h i) = M33
        (k * a) (k * b) (k * c)
        (k * d) (k * e) (k * f)
        (k * g) (k * h) (k * i)
    {-# INLINE mult #-}

instance Mult M33 V3 V3 where
    mult (M33 a b c
           d e f
           g h i) (V3 u v w) = V3
        (a * u + b * v + c * w)
        (d * u + e * v + f * w)
        (g * u + h * v + i * w)
    {-# INLINE mult #-}

-- TODO: instance Mult M33 M33 M33 where

det3 :: M33 -> Double
det3 (M33 a b c
          d e f
          g h i)
    = a * (e*i-f*h) - d * (b*i-c*h) + g * (b*f-c*e)
{-# INLINE det3 #-}

inv3 :: M33 -> M33
inv3 m@(M33 a b c
            d e f
            g h i)
    = M33 a' b' c'
          d' e' f'
          g' h' i'
  where
    a' = cofactor e f h i / det
    b' = cofactor c b i h / det
    c' = cofactor b c e f / det
    d' = cofactor f d i g / det
    e' = cofactor a c g i / det
    f' = cofactor c a f d / det
    g' = cofactor d e g h / det
    h' = cofactor b a h g / det
    i' = cofactor a b d e / det
    cofactor q r s t = det2 (M22 q r s t)
    det = det3 m
{-# INLINE inv3 #-}

-------------------------------------------------------------------------------
-- Regressions
-------------------------------------------------------------------------------

-- | Linear regression.
--
-- The type is
--
-- @
-- 'linear' :: [('Double', 'Double')] -> 'V2'
-- @
--
-- but overloaded to work with boxed and unboxed 'Vector's.
--
-- >>> let input1 = [(0, 1), (1, 3), (2, 5)]
-- >>> PP $ linear input1
-- V2 2.0000 1.00000
--
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> PP $ linear input2
-- V2 2.0063 0.88685
--
linear :: (Foldable' xs x, IsDoublePair x) => xs -> V2
linear = fst . linearWithErrors

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
-- >>> PP $ linearWithErrors input1
-- (V2 2.0000 1.00000, V2 0.00000 0.00000)
--
-- The second input is quite good:
--
-- >>> PP $ linearWithErrors input2
-- (V2 2.0063 0.88685, V2 0.09550 0.23826)
--
-- But the third input isn't so much,
-- standard error of a slope argument is 20%.
--
-- >>> let input3 = [(0, 2), (1, 3), (2, 6), (3, 11)]
-- >>> PP $ linearWithErrors input3
-- (V2 3.0000 1.00000, V2 0.63246 1.1832)
--
-- @since 0.1.1
--
linearWithErrors :: (Foldable' xs x, IsDoublePair x) => xs -> (V2, V2)
linearWithErrors = linearImpl . kahan2

linearImpl :: Kahan2 -> (V2, V2)
linearImpl (K2 n' (V2 x _) (V2 x2 _) (V2 y _) (V2 y2 _) (V2 xy _)) =
    (params, errors)
  where
    n :: Double
    n = fromIntegral n'

    matrix@(M22 a11 _ _ a22) = inv2 (M22 x2 x x n)
    params@(V2 a b)          = mult matrix (V2 xy y)

    errors = V2 sa sb

    -- ensure that error is always non-negative.
    -- Due rounding errors, in perfect fit situations it can be slightly negative.
    err = max 0 (y2 - a * xy - b * y)
    sa  = sqrt (a11 * err / (n - 2))
    sb  = sqrt (a22 * err / (n - 2))

-- | Quadratic regression.
--
-- The type is
--
-- @
-- 'quadratic' :: [('Double', 'Double')] -> 'V3'
-- @
--
-- but overloaded to work with boxed and unboxed 'Vector's.
--
-- >>> let input1 = [(0, 1), (1, 3), (2, 5)]
-- >>> quadratic input1
-- V3 0.0 2.0 1.0
--
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> PP $ quadratic input2
-- V3 (-0.00589) 2.0313 0.87155
--
-- >>> let input3 = [(0, 2), (1, 3), (2, 6), (3, 11)]
-- >>> PP $ quadratic input3
-- V3 1.00000 0.00000 2.0000
--
quadratic :: (Foldable' xs x, IsDoublePair x) => xs -> V3
quadratic data_ = mult (inv3 (M33 x2 x n x3 x2 x x4 x3 x2)) (V3 y xy x2y)
  where
    K3 n' (V2 x _) (V2 x2 _) (V2 x3 _) (V2 x4 _) (V2 y _) (V2 xy _) (V2 x2y _) (V2 _y2 _) = kahan3 data_

    n :: Double
    n = fromIntegral n'

-- | Like 'quadratic' but also return parameters' standard errors.
--
-- >>> PP $ quadraticWithErrors input2
-- (V3 (-0.00589) 2.0313 0.87155, V3 0.09281 0.41070 0.37841)
--
-- >>> PP $ quadraticWithErrors input3
-- (V3 1.00000 0.00000 2.0000, V3 0.00000 0.00000 0.00000)
--
-- @since 0.1.1
--
quadraticWithErrors :: (Foldable' xs x, IsDoublePair x) => xs -> (V3, V3)
quadraticWithErrors = quadraticImpl . kahan3

quadraticImpl :: Kahan3 -> (V3, V3)
quadraticImpl (K3 n' (V2 x _) (V2 x2 _) (V2 x3 _) (V2 x4 _) (V2 y _) (V2 xy _) (V2 x2y _) (V2 y2 _)) =
    (params, errors)
  where
    n :: Double
    n = fromIntegral n'

    matrix@(M33 a11 _   _
                _   a22 _
                _   _   a33) = inv3 (M33 x4 x3 x2 x3 x2 x x2 x n)

    params@(V3 a b c) = mult matrix (V3 x2y xy y)

    errors = V3 sa sb sc

    err = max 0 (y2 - a * x2y - b * xy - c * y)
    sa  = sqrt (a11 * err / (n - 3))
    sb  = sqrt (a22 * err / (n - 3))
    sc  = sqrt (a33 * err / (n - 3))

-- | Do both linear and quadratic regression in one data scan.
--
-- >>> PP $ quadraticAndLinear input2
-- (V3 (-0.00589) 2.0313 0.87155, V2 2.0063 0.88685)
--
quadraticAndLinear :: (Foldable' xs x, IsDoublePair x) => xs -> (V3, V2)
quadraticAndLinear = fst . quadraticAndLinearWithErrors

-- | Like 'quadraticAndLinear' but also return parameters' standard errors
--
-- >>> PP $ quadraticAndLinearWithErrors input2
-- ((V3 (-0.00589) 2.0313 0.87155, V2 2.0063 0.88685), (V3 0.09281 0.41070 0.37841, V2 0.09550 0.23826))
--
-- @since 0.1.1
--
quadraticAndLinearWithErrors :: (Foldable' xs x, IsDoublePair x) => xs -> ((V3, V2), (V3, V2))
quadraticAndLinearWithErrors data_ =
    ((paramsQ, paramsL), (errorsQ, errorsL))
  where
    k3@(K3 n x x2 x3 x4 y xy x2y y2) = kahan3 data_
    k2 = K2 n x x2 y y2 xy

    (paramsL, errorsL) = linearImpl k2
    (paramsQ, errorsQ) = quadraticImpl k3

-------------------------------------------------------------------------------
-- Input
-------------------------------------------------------------------------------

-- | Like 'Foldable' but with element in the class definition.
class Foldable' xs x | xs -> x where
    foldl' :: (b -> x -> b) -> b -> xs -> b

instance              Foldable' [a]          a where foldl' = L.foldl'
instance              Foldable' (V.Vector a) a where foldl' = V.foldl'
instance U.Unbox a => Foldable' (U.Vector a) a where foldl' = U.foldl'

-- | Class witnessing that @dp@ has a pair of 'Double's.
class IsDoublePair dp where
    withDP :: dp -> (Double -> Double -> r) -> r
    makeDP :: Double -> Double -> dp

instance IsDoublePair V2 where
    withDP (V2 x y) k = k x y
    makeDP = V2

instance (a ~ Double, b ~ Double) => IsDoublePair (a, b) where
    withDP ~(x, y) k = k x y
    makeDP = (,)

-------------------------------------------------------------------------------
-- Kahan2
-------------------------------------------------------------------------------

data Kahan2 = K2
    { k2n  :: {-# UNPACK #-} !Int
    , k2x  :: {-# UNPACK #-} !V2
    , k2x2 :: {-# UNPACK #-} !V2
    , k2y  :: {-# UNPACK #-} !V2
    , k2y2 :: {-# UNPACK #-} !V2
    , k2xy :: {-# UNPACK #-} !V2
    }

zeroKahan2 :: Kahan2
zeroKahan2 = K2 0 zero zero zero zero zero

-- | https://en.wikipedia.org/wiki/Kahan_summation_algorithm
addKahan :: V2 -> Double -> V2
addKahan (V2 acc c) i =
    let y = i - c
        t = acc + y
    in V2 t ((t - acc) - y)

kahan2 :: (Foldable' xs x, IsDoublePair x) => xs -> Kahan2
kahan2 = foldl' f zeroKahan2 where
    f (K2 n x x2 y y2 xy) uv = withDP uv $ \u v -> K2
        (succ n)
        (addKahan x u)
        (addKahan x2 (u * u))
        (addKahan y v)
        (addKahan y2 (v * v))
        (addKahan xy (u * v))

-------------------------------------------------------------------------------
-- Kahan3
-------------------------------------------------------------------------------

data Kahan3 = K3
    { k3n   :: {-# UNPACK #-} !Int
    , k3x   :: {-# UNPACK #-} !V2
    , k3x2  :: {-# UNPACK #-} !V2
    , k3x3  :: {-# UNPACK #-} !V2
    , k3x4  :: {-# UNPACK #-} !V2
    , k3y   :: {-# UNPACK #-} !V2
    , k3xy  :: {-# UNPACK #-} !V2
    , k3x2y :: {-# UNPACK #-} !V2
    , ky2   :: {-# UNPACK #-} !V2
    }

zeroKahan3 :: Kahan3
zeroKahan3 = K3 0 zero zero zero zero zero zero zero zero

kahan3 :: (Foldable' xs x, IsDoublePair x) => xs -> Kahan3
kahan3 = foldl' f zeroKahan3 where
    f (K3 n x x2 x3 x4 y xy x2y y2) uv = withDP uv $ \u v ->
        let u2 = u * u
        in K3
            (succ n)
            (addKahan x u)
            (addKahan x2 u2)
            (addKahan x3 (u * u2))
            (addKahan x4 (u2 * u2))
            (addKahan y v)
            (addKahan xy (u * v))
            (addKahan x2y (u2 * v))
            (addKahan y2 (v * v))
