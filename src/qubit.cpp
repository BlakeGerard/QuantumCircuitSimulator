/*
        Title: qubit.cpp
        Author: Blake Gerard
        Date: 05/16/2020
        Description: Implementation of the Qubit class.
        See "qubit.h" for more information.
*/

#include <complex>
#include "qubit.h"

/*
    Default consuctor. Qubit initialized to [1 0]^T state.
*/
Qubit::Qubit() {
    state << 1.0, 0.0;
};

/*
    Specified probabilities consuctor. Qubit initialized to [a b]^T state.
*/
Qubit::Qubit(std::complex<double> a, std::complex<double> b) {
    state << a, b;
    
    if (state.norm() != 1.0) {
        state.normalize();
    };
};

/*
    Return the qubit's current state.

    Params:
        None
    Return:
        state: Current state of the qubit.
*/
Eigen::Vector2cd Qubit::get_state() {
    return state;
};

/*
    Set the qubit's current state with new a and b values.

    Params:
        a: probability amplitude
        b: probability amplitude
    Return:
        None
*/
void Qubit::set_state(std::complex<double> a, std::complex<double> b) {
    state[0] = a;
    state[1] = b;
};