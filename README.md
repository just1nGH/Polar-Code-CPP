# Polar-Code-CPP
Polar Code C++ Implementation

try out *`example_run.m`  which shows how to use the implemented polar codes.*

## Features

C++ 17 needed to run the examples and use the code.  No additional package is used, only standard library.

The functions are encapsulated in a class for easy use.

### code construction ：

- channel polarization using Huawei approximation [3]

### Encoding

- 3GPP encoding (without bit-reversal permutation)  $\mathbf x = \mathbf u \mathbf F_n$

### Decoding

- SC decoder
- CRC-aided SC decoder

### Rate matching and rate recovery

- to support different codeword length and code rate
