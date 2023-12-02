# SQTV2
Qiskit-like quantum algorithm tool
This is an academic project i undertook to better understand the mathematical framework of quantum algorithms.

Implements some of the basic functionality of qiskit circuit objects with similar syntax.
Has support for any 2D unitary operation and the application of such an operation based on one control bit.
A few commonly used gates are implemented as standard and can be quickly accesed. These are hadamard (h), phase (p) and pauli gates (xyz). Others can be created and then saved for quick access.

Create a cicruit using the quantum_circuit (qc) object and specify the number of bits.
Use qc.set_register to decide the initial register state. It must be a computational basis vector. e.g. "000111" for a 6 qubit circuit.

Gates are added left to write creating a new layer if the user places a gate on a qubit that is already populated on the current layer. Control gates always create a new layer.

Add basic gates to the quantum circuit using the functions qc.X. E.g. qc.h(0) to add a hadamard gate on the first quebit.

Save custom gates using the qc.create_unitary function, specify the ascii representation using the character option. Apply them to the circuit using the qc.su function. Alternatively apply a new unitary directly without saving it by using the qc.u function.

Apply basic control gates using the qc.cX method. e.g. qc.ch(0,1) for a controlled application of hadamard on the second qubit.

Similarly use qc.cu or qc.scu for controlled application of not saved / saved quantum gates.

use qc.print_circuit to output an ascii representation of the circuit.

use qc.show_results() to output the probability of measuring the different computational basis states after running the circuit. Use qc.get_results to get these probabilities in a vector format.

LIMITATIONS
The code does not directly support gates with multiple control bits e.g. toffoli. This makes it tedious to implement many commonly used quantum algorithms, since these gates have to be built manually out of simpler ones. 

The code works by representing the circuit using a complex matrix of size NxN where N is number of possible states of the register (N=2^n where n is the number of qubits). Thus the computational complexity scales exponentially with the number of qubits. The code could likely be made much more efficient, but the scaling issue is not fixable. This means that even if the code was made twice as efficient, this would only allow for the simulation of one extra qubit. If quantum emulation was a tractable problem, then there would be no benefit to building quantum computers in the first place! From my testing upwards of 9 qubits is doable on a personal computer.

Also there might be memory leaks, and i am not at all proud of my tensor product code.
