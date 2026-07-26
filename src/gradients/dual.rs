//! Tiny fixed-width forward-mode scalar used by analytical kernels.

use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::core::Scalar;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Dual<const N: usize> {
    pub value: Scalar,
    pub derivative: [Scalar; N],
}

impl<const N: usize> Dual<N> {
    pub const fn constant(value: Scalar) -> Self {
        Self {
            value,
            derivative: [0.0; N],
        }
    }

    pub fn variable(value: Scalar, index: usize) -> Self {
        let mut derivative = [0.0; N];
        derivative[index] = 1.0;
        Self { value, derivative }
    }

    pub fn sqrt(self) -> Self {
        let value = self.value.sqrt();
        let scale = 0.5 / value;
        Self {
            value,
            derivative: self.derivative.map(|item| item * scale),
        }
    }
}

impl<const N: usize> From<Scalar> for Dual<N> {
    fn from(value: Scalar) -> Self {
        Self::constant(value)
    }
}

impl<const N: usize> Add for Dual<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value + rhs.value,
            derivative: std::array::from_fn(|index| self.derivative[index] + rhs.derivative[index]),
        }
    }
}

impl<const N: usize> Sub for Dual<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value - rhs.value,
            derivative: std::array::from_fn(|index| self.derivative[index] - rhs.derivative[index]),
        }
    }
}

impl<const N: usize> Mul for Dual<N> {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value * rhs.value,
            derivative: std::array::from_fn(|index| {
                self.derivative[index] * rhs.value + self.value * rhs.derivative[index]
            }),
        }
    }
}

impl<const N: usize> Div for Dual<N> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let denominator = rhs.value * rhs.value;
        Self {
            value: self.value / rhs.value,
            derivative: std::array::from_fn(|index| {
                (self.derivative[index] * rhs.value - self.value * rhs.derivative[index])
                    / denominator
            }),
        }
    }
}

impl<const N: usize> Neg for Dual<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            value: -self.value,
            derivative: self.derivative.map(|item| -item),
        }
    }
}
