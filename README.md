# dmmbuild
![](https://img.shields.io/badge/license-MIT-brightgreen.svg)

## What is DUMMER?

DUMMER (Dumb Uncomplicated Match ModelER) aims to find distant
relationships between genetic sequences (nucleotide or protein).  It's
similar to [HMMER](http://hmmer.org), but much simpler, and aspires
to be better.

In more detail, it finds similar regions between sequences and
"profiles".  A profile is a set of position-specific letter, deletion,
and insertion probabilities: typically made from a family of related
sequences.

This is a proof-of-principle for the paper [A simple way to find
related sequences with position-specific probabilities](https://doi.org/10.1101/2025.03.14.643233).

## What is dmmbuild?

Originally, DUMMER used HMMER-generated "profiles" as an input.
However, HMMER defines insertion and deletion probabilities slightly
differently than DUMMER, so they are suboptimal for its algorithms.

**dmmbuild** is supposed to change that. It's based on HMMER's
**hmmbuild**, and has been adapted to build profiles for DUMMER.
Essentially, as opposed to the original *plan-7* profile Hidden
Markov Model, we aim to build the following *fig-4* HMM:

<img width="872" height="826" alt="image" src="https://github.com/user-attachments/assets/489fcb3e-a2c0-4db7-9074-6f820757a32d" />

## How do I use hmmbuild?

As the original HMMER project, we depend on Easel. You'll need to clone both the
DUMMER and Easel repositories, configure Easel and build both, as follows:

```bash
   % git clone https://github.com/Padraig20/hmmbuild
   % cd hmmbuild
   % git clone https://github.com/EddyRivasLab/easel
   % cd easel
   % autoconf
   % ./configure
   % cd ..
   % make
```

Afterwards, you should be good to go! (If not, please feel free to reach out!)
