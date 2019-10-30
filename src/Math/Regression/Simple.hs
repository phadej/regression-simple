{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE GADTs                  #-}
module Math.Regression.Simple (
    -- * Regressions
    linear,
    quadratic,
    quadraticAndLinear,
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
    HasDoublePair (..),
    ) where

import Data.Complex (Complex (..))

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

-- $setup
-- >>> :set -XTypeApplications

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
zerosLin (V2 a b) = - b / a

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
-- >>> linear input1
-- V2 2.0 1.0
--
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> linear input2
-- V2 2.0063237774030345 0.8868465430016883
--
linear :: (Foldable' xs x, HasDoublePair x) => xs -> V2
linear data_ = mult (inv2 (M22 x n x2 x)) (V2 y xy)
  where
    K2 n' (V2 x _) (V2 x2 _) (V2 y _) (V2 xy _) = kahan2 data_
    n :: Double
    n = fromIntegral n'

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
-- >>> quadratic input2
-- V3 (-5.886346291028133e-3) 2.0312938469708826 0.8715454176158062
--
-- >>> let input3 = [(0, 2), (1, 3), (2, 6), (3, 11)]
-- >>> quadratic input3
-- V3 1.0 0.0 1.999999999999993
--
quadratic :: (Foldable' xs x, HasDoublePair x) => xs -> V3
quadratic data_ = mult (inv3 (M33 x2 x n x3 x2 x x4 x3 x2)) (V3 y xy x2y)
  where
    K3 n' (V2 x _) (V2 x2 _) (V2 x3 _) (V2 x4 _) (V2 y _) (V2 xy _) (V2 x2y _) = kahan3 data_
    n :: Double
    n = fromIntegral n'

-- | Do both linear and quadratic regression in one data scan.
--
-- >>> let input2 = [(0.1, 1.2), (1.3, 3.1), (1.9, 4.9), (3.0, 7.1), (4.1, 9.0)]
-- >>> quadraticAndLinear input2
-- (V3 (-5.886346291028133e-3) 2.0312938469708826 0.8715454176158062,V2 2.0063237774030345 0.8868465430016883)
--
quadraticAndLinear :: (Foldable' xs x, HasDoublePair x) => xs -> (V3, V2)
quadraticAndLinear data_ =
    ( mult (inv3 (M33 x2 x n x3 x2 x x4 x3 x2)) (V3 y xy x2y)
    , mult (inv2 (M22 x n x2 x)) (V2 y xy)
    )
  where
    K3 n' (V2 x _) (V2 x2 _) (V2 x3 _) (V2 x4 _) (V2 y _) (V2 xy _) (V2 x2y _) = kahan3 data_
    n :: Double
    n = fromIntegral n'

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
class HasDoublePair dp where
    withDP :: dp -> (Double -> Double -> r) -> r

instance HasDoublePair V2 where
    withDP (V2 x y) k = k x y

instance (a ~ Double, b ~ Double) => HasDoublePair (a, b) where
    withDP ~(x, y) k = k x y

-------------------------------------------------------------------------------
-- Kahan2
-------------------------------------------------------------------------------

data Kahan2 = K2
    { k2n  :: {-# UNPACK #-} !Int
    , k2x  :: {-# UNPACK #-} !V2
    , k2x2 :: {-# UNPACK #-} !V2
    , k2y  :: {-# UNPACK #-} !V2
    , k2xy :: {-# UNPACK #-} !V2
    }

zeroKahan2 :: Kahan2
zeroKahan2 = K2 0 zero zero zero zero

-- | https://en.wikipedia.org/wiki/Kahan_summation_algorithm
addKahan :: V2 -> Double -> V2
addKahan (V2 acc c) i =
    let y = i - c
        t = acc + y
    in V2 t ((t - acc) - y)

kahan2 :: (Foldable' xs x, HasDoublePair x) => xs -> Kahan2
kahan2 = foldl' f zeroKahan2 where
    f (K2 n x x2 y xy) uv = withDP uv $ \u v -> K2
        (succ n)
        (addKahan x u)
        (addKahan x2 (u * u))
        (addKahan y v)
        (addKahan xy (u * v))

-------------------------------------------------------------------------------
-- Kahan3
-------------------------------------------------------------------------------

data Kahan3 = K3
    { k3n  :: {-# UNPACK #-} !Int
    , k3x  :: {-# UNPACK #-} !V2
    , k3x2 :: {-# UNPACK #-} !V2
    , k3x3 :: {-# UNPACK #-} !V2
    , k3x4 :: {-# UNPACK #-} !V2
    , k3y  :: {-# UNPACK #-} !V2
    , k3xy :: {-# UNPACK #-} !V2
    , k3x2y :: {-# UNPACK #-} !V2
    }

zeroKahan3 :: Kahan3
zeroKahan3 = K3 0 zero zero zero zero zero zero zero

kahan3 :: (Foldable' xs x, HasDoublePair x) => xs -> Kahan3
kahan3 = foldl' f zeroKahan3 where
    f (K3 n x x2 x3 x4 y xy x2y) uv = withDP uv $ \u v ->
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
