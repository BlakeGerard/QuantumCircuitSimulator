/*
        Title: gate_enum.h
        Author: Blake Gerard
        Date: 05/16/2020
        Description: Enumeration of basic quantum gates. Functions
        requiring gates of higher dimension (CNOT, CCNOT, etc.) 
        simply refer to their corresponding single-qubit versions.
        These values act as indices into the static array of 
        matrix-representations of gates that all circuits utilize.
*/

#ifndef GATE_ENUM_H
#define GATE_ENUM_H
enum Gate {
        H,      // Hadamard
        X,      // Used for X, CNOT, and CCNOT
        Y,      // Used for Y, CY
        Z,      // Used for Z, CZ
        S,      // S gate
        T,      // T gate
        SWAP,   // Swap gate
};
#endif