#ifndef QUANTUM_GATES_H
#define QUANTUM_GATES_H
#include "complex_matrix.h"
#include <vector>
#include <complex>
#include <memory>
#include <stdio.h>
#include <wchar.h>
#include <windows.h>
#include <string>

class quantum_gate
{
public:
	virtual matrix get_matrix() = 0; //any quantum gate has a matrix representation
};

class single_bit_gate : public quantum_gate //this class is for single bit gate e.g. hadamard and so on
{
protected:
	matrix matrix_rep{ 2,2 };
	std::string symbol;
public:
	virtual ~single_bit_gate() {}
	matrix get_matrix() { return matrix_rep; }
	std::string get_symbol() { return symbol; }
};

class full_lane_gate : public quantum_gate //this class is meant for composite gates and control gates that return matrices for operating on the full tensor product space
{
public:
	virtual void add_gate(int qbit_number, const std::shared_ptr<single_bit_gate> & gate_pointer) {} //these two are used for interacting with the parallel gate object through the quantum_circuit object
	virtual std::shared_ptr<single_bit_gate> read_gate(int qbit_number) { return nullptr; }
	virtual void print_gate() = 0;
};


class control_gate : public full_lane_gate
{
private:
	int number_of_qbits; //This object has a lot to keep track of
	int control_bit;
	int controlled_bit;
	std::shared_ptr<matrix> on_projection_pointer{ nullptr };
	std::shared_ptr<matrix> off_projection_pointer{ nullptr };
	std::shared_ptr<single_bit_gate> identity_pointer{ nullptr };
	std::shared_ptr<single_bit_gate> control_pointer{ nullptr };
public:
	control_gate(int number_of_qbits, int control_bit, int controlled_bit, //A bit clunky but readable
		std::shared_ptr<matrix> on_projection_pointer, std::shared_ptr<matrix> off_projection_pointer,
		std::shared_ptr<single_bit_gate> identity_pointer, std::shared_ptr<single_bit_gate> control_pointer) :
		number_of_qbits{ number_of_qbits }, control_bit{ control_bit }, controlled_bit{ controlled_bit },
		on_projection_pointer{ on_projection_pointer }, off_projection_pointer{ off_projection_pointer },
		identity_pointer{ identity_pointer }, control_pointer{ control_pointer }
	{
		if (control_bit >= number_of_qbits || control_bit < 0) { std::cout << "Control bit out of range" << std::endl; exit(1); }
		if (controlled_bit >= number_of_qbits || controlled_bit < 0) { std::cout << "Controlled bit out of range" << std::endl; exit(1); }
		if (control_bit == controlled_bit) { std::cout << "bit cannot control itself!" << std::endl; exit(1); }
	}
	matrix get_matrix()
	{
		matrix off_projection;
		matrix on_projection;
		for (int i{ number_of_qbits - 1 }; i >= 0; i--) {
			if (i == control_bit) {
				off_projection = tensor_product(off_projection, *off_projection_pointer);
				on_projection = tensor_product(on_projection, *on_projection_pointer);
			}
			else if (i == controlled_bit) {
				off_projection = tensor_product(off_projection, identity_pointer->get_matrix());
				on_projection = tensor_product(on_projection, control_pointer->get_matrix());
			}
			else {
				off_projection = tensor_product(off_projection, identity_pointer->get_matrix());
				on_projection = tensor_product(on_projection, identity_pointer->get_matrix());
			}
		}
		return off_projection + on_projection;
	}

	void print_gate() {
		for (int i{ 0 }; i < number_of_qbits; i++) {
			if ((i < controlled_bit && i < control_bit) || (i > controlled_bit && i > control_bit)) { std::cout << "--"; }
			else if (i == controlled_bit) { std::cout << control_pointer->get_symbol() << '-'; }
			else if (i == control_bit) { std::cout << "C-"; }
			else { std::cout << "|-"; }
			wprintf(L"\x1b[1B");
			wprintf(L"\x1b[2D");
		}
	}
};

class wall : public full_lane_gate // its called wall instead of barrier because barrier is a standard template.
{
	int number_of_qbits;
public:
	wall(int number_of_qbits) : number_of_qbits{ number_of_qbits } {}

	matrix get_matrix() //the barrier does nothing, so we can just represent it with an identity matrix
	{
		matrix matrix_rep{ (int)std::pow(2,number_of_qbits), (int)std::pow(2,number_of_qbits) };
		for (int i = 1; i <= matrix_rep.get_rows(); i++) {
			for (int j = 1; j <= matrix_rep.get_cols(); j++) {
				if (i == j) matrix_rep(i, j) = 1;
				else matrix_rep(i, j) = 0;
			}
		}
		return matrix_rep;
	}

	void print_gate() {
		for (int i{ 0 }; i < number_of_qbits; i++) {
			std::cout << '#' << '-';
			wprintf(L"\x1b[1B");
			wprintf(L"\x1b[2D");

		}
	}

};

//This object is used to keep track of single qbit gates applied in parallel.
class parallel_gates : public full_lane_gate
{
public:
	int number_of_qbits;
	std::shared_ptr<single_bit_gate>* gate_array{ nullptr };
public:
	//Has no default constructor since a quantum circuit always has a size 
	//constructor always takes a pointer to an identity matrix as these are used to populate any empty spaces in the gate array.
	//You can also give reference to other gate types, e.g. to quickly make a wall of hadamards.
	parallel_gates(int number_of_qbits, std::shared_ptr<single_bit_gate> norm_pointer) : number_of_qbits{ number_of_qbits }
	{
		gate_array = new std::shared_ptr<single_bit_gate> [number_of_qbits];
		for (int i{ 0 }; i < number_of_qbits; i++) {
			gate_array[i] = norm_pointer;
		}
	}

	~parallel_gates() { delete[] gate_array; }

	//function for populating specific gate types to certain elements
	void add_gate(int qbit_number, const std::shared_ptr<single_bit_gate> & gate_pointer) {
		if (qbit_number < 0 || qbit_number >= number_of_qbits) { std::cout << "qbit doesn't exist" << std::endl; exit(1); }
		else gate_array[qbit_number] = gate_pointer;
	}

	//This matrix manually grinds out tesnor product iteratively
	matrix get_matrix() {

		matrix matrix_rep{ 0,0 };
		for (int i{ number_of_qbits - 1 }; i >= 0; i--) {
			matrix_rep = tensor_product(matrix_rep, gate_array[i]->get_matrix());
		}
		return matrix_rep;

	}

	std::shared_ptr<single_bit_gate> read_gate(int qbit_number) { return gate_array[qbit_number]; }

	void print_gate() {
		for (int i{ 0 }; i < number_of_qbits; i++) {
			std::cout << gate_array[i]->get_symbol() << '-';
			wprintf(L"\x1b[1B");
			wprintf(L"\x1b[2D");

		}
	}
};

//basic gates
class hadamard : public single_bit_gate
{
public:
	hadamard()
	{
		matrix_rep(1, 1) = 1; matrix_rep(1, 2) = 1;
		matrix_rep(2, 1) = 1; matrix_rep(2, 2) = -1;
		matrix_rep /= sqrt(2);

		symbol = "H";
	}
};
class identity : public single_bit_gate
{

public:
	identity()
	{
		matrix_rep(1, 1) = 1; matrix_rep(1, 2) = 0;
		matrix_rep(2, 1) = 0; matrix_rep(2, 2) = 1;

		symbol = "-";
	}
};
class pauli_x : public single_bit_gate
{

public:
	pauli_x()
	{
		matrix_rep(1, 1) = 0; matrix_rep(1, 2) = 1;
		matrix_rep(2, 1) = 1; matrix_rep(2, 2) = 0;

		symbol = "X";
	}
};
class pauli_y : public single_bit_gate
{

public:
	pauli_y()
	{
		matrix_rep(1, 1) = 0; matrix_rep(1, 2) = std::complex<double>{ 0,-1 };
		matrix_rep(2, 1) = std::complex<double>{ 0,1 }; matrix_rep(2, 2) = 0;

		symbol = "Y";
	}
};
class pauli_z : public single_bit_gate
{

public:
	pauli_z()
	{
		matrix_rep(1, 1) = 1; matrix_rep(1, 2) = 0;
		matrix_rep(2, 1) = 0; matrix_rep(2, 2) = -1;

		symbol = "Z";
	}
};
class phase : public single_bit_gate
{
private:
	std::string symbol{ "P" };
public:
	phase()
	{
		matrix_rep(1, 1) = 1; matrix_rep(1, 2) = 0;
		matrix_rep(2, 1) = 0; matrix_rep(2, 2) = std::complex<double>{ 0,1 };

		symbol = "P";
	}
};

class u_gate : public single_bit_gate { //This is a completely generalized unitary operator
private:
	double theta{ 0 }; // If no angles are specified
	double phi{ 0 };   // Then this just becomes the identity
	double lambda{ 0 };// Angles need to be in radians
public:
	u_gate(double theta, double phi, double lambda, char character = 'U') :
		theta{ theta }, phi{ phi }, lambda{ lambda }
	{
		matrix_rep(1, 1) = cos(theta / 2);	matrix_rep(1, 2) = std::polar(1.0, lambda) * sin(theta / 2) * (-1.0);
		matrix_rep(2, 1) = std::polar(1.0, phi) * sin(theta / 2);	matrix_rep(2, 2) = std::polar(1.0, phi + lambda) * cos(theta / 2);

		symbol = character;
	}
};



#endif
