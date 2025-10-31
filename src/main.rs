use clap::Parser;
use num_bigint::{BigUint, RandBigInt};
use num_integer::Integer;
use num_traits::{One, Zero};
use rand::rngs::OsRng;
use std::convert::TryFrom;

/// Default Miller–Rabin rounds. Increase for extra certainty.
const DEFAULT_MR_ROUNDS: usize = 64;

/// CLI arguments parsed via clap.
#[derive(Parser, Debug)]
#[command(name = "find-big-prime", about = "Generate large probable primes and safe primes")]
struct Args {
    /// Number of bits for the generated prime (e.g. 2048, 3072, 4096).
    #[arg(short = 'b', long = "bits", default_value_t = 2048)]
    bits: usize,

    /// Generate a safe prime p where p = 2q + 1 and q is also prime.
    #[arg(long = "safe")]
    safe: bool,

    /// Miller–Rabin rounds to run when testing primality.
    #[arg(long = "rounds", default_value_t = DEFAULT_MR_ROUNDS)]
    rounds: usize,
}

fn main() {
    let args = Args::parse();
    assert!(
        args.bits >= 512,
        "At least 512 bits are recommended; use >= 2048 bits for production."
    );

    if args.safe {
        let p = generate_safe_prime(args.bits, args.rounds);
        println!("safe_prime_bits={}", p.bits());
        println!("{p}");
    } else {
        let p = generate_probable_prime(args.bits, args.rounds);
        println!("prime_bits={}", p.bits());
        println!("{p}");
    }
}

/// Generate a random probable prime with the requested bit length.
fn generate_probable_prime(bits: usize, rounds: usize) -> BigUint {
    let mut rng = OsRng;
    let bits_u64 = u64::try_from(bits).expect("bit size must fit in u64");
    loop {
        let mut n = rng.gen_biguint(bits_u64);
        let one = BigUint::one();

        // Force highest bit to ensure bit length and make the candidate odd.
        n.set_bit(bits_u64 - 1, true);
        if n.is_even() {
            n |= &one;
        }

        if !small_prime_precheck(&n) {
            continue;
        }

        if is_probable_prime(&n, rounds) {
            return n;
        }
    }
}

/// Generate a safe prime p = 2q + 1 where both p and q are probable primes.
fn generate_safe_prime(bits: usize, rounds: usize) -> BigUint {
    assert!(bits >= 3, "Safe primes require at least 3 bits.");
    let q_bits = bits - 1;
    loop {
        let q = generate_probable_prime(q_bits, rounds);
        let p = (&q << 1usize) + BigUint::one();
        if is_probable_prime(&p, rounds) {
            return p;
        }
    }
}

/// Miller–Rabin probabilistic primality test.
fn is_probable_prime(n: &BigUint, rounds: usize) -> bool {
    let two = BigUint::from(2u32);

    if *n < two {
        return false;
    }
    if *n == two {
        return true;
    }
    if n.is_even() {
        return false;
    }

    let one = BigUint::one();
    let n_minus_one = n - &one;
    let (s, d) = factor_out_twos(&n_minus_one);

    let mut rng = OsRng;
    'witness: for _ in 0..rounds {
        let a = random_range(&two, &n_minus_one, &mut rng);
        let mut x = a.modpow(&d, n);

        if x == one || x == n_minus_one {
            continue 'witness;
        }

        for _ in 1..s {
            x = x.modpow(&two, n);
            if x == n_minus_one {
                continue 'witness;
            }
            if x == one {
                return false;
            }
        }

        return false;
    }

    true
}

/// Express n as d * 2^s with d odd, returning (s, d).
fn factor_out_twos(n: &BigUint) -> (u32, BigUint) {
    let mut s = 0u32;
    let mut d = n.clone();
    while d.is_even() {
        d >>= 1;
        s += 1;
    }
    (s, d)
}

/// Sample a random value in the inclusive range [low, high].
fn random_range(low: &BigUint, high: &BigUint, rng: &mut OsRng) -> BigUint {
    if low == high {
        return low.clone();
    }
    let high_exclusive = high + BigUint::one();
    rng.gen_biguint_range(low, &high_exclusive)
}

/// Filter out obvious composites using a small set of primes.
fn small_prime_precheck(n: &BigUint) -> bool {
    const SMALLS: [u32; 16] = [
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    ];

    if n == &BigUint::one() {
        return false;
    }

    for &p in SMALLS.iter() {
        let p_big = BigUint::from(p);
        if n == &p_big {
            return true;
        }
        if (n % &p_big).is_zero() {
            return false;
        }
    }

    true
}
