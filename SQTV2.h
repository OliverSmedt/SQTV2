#ifndef SQTV2_H
#define SQTV2_H

#include <complex>
#include <memory>
#include "complex_matrix.h"
#include "quantum_gates.h"
#include <map>
#include <string>
#include <bitset>
#include <vector>
#include <stdio.h>
#include <wchar.h>
#include <windows.h>

class quantum_circuit : full_lane_gate{
public:
	const int number_of_qbits;
	matrix state_vector;

	//basic single qbit gates
	std::shared_ptr<single_bit_gate> ptrH = std::shared_ptr<single_bit_gate>(new hadamard);
	std::shared_ptr<single_bit_gate> ptrX = std::shared_ptr<single_bit_gate>(new pauli_x);
	std::shared_ptr<single_bit_gate> ptrY = std::shared_ptr<single_bit_gate>(new pauli_y);
	std::shared_ptr<single_bit_gate> ptrZ = std::shared_ptr<single_bit_gate>(new pauli_z);
	std::shared_ptr<single_bit_gate> ptrI = std::shared_ptr<single_bit_gate>(new identity);
	std::shared_ptr<single_bit_gate> ptrP = std::shared_ptr<single_bit_gate>(new phase);


	//This vector contains every step in the circuit emulation it consists of parallel and control gates
	std::vector<std::shared_ptr<full_lane_gate>> gate_vector;

	//This vector is for any custom unitary gates the user may define
	std::vector<std::shared_ptr<single_bit_gate>> custom_gates;
	std::map<std::string, int> gate_map{}; //this is used to name user defined unitary gates

	//This vector is for circuit components
	std::vector<std::shared_ptr<full_lane_gate>> custom_circuits;
	std::map<std::string, int> circuit_map{};

	//defining projection matrices

	std::shared_ptr<matrix> ptrON = std::shared_ptr<matrix>(new matrix{2,2});
	std::shared_ptr<matrix> ptrOFF = std::shared_ptr<matrix>(new matrix{2,2});

	char character;

public:
	quantum_circuit(int number_of_qbits, char character = 'C') : number_of_qbits{ number_of_qbits }, character{character}
	{
		//making pointers to all the standard gate types

		//building projection matrices, these are needed for generating the matrix representations of arbitrary control gates
		ptrON->operator()(1, 1) = 0; ptrON->operator()(1, 2) = 0;
		ptrON->operator()(2, 1) = 0; ptrON->operator()(2, 2) = 1;

		ptrOFF->operator()(1, 1) = 1; ptrOFF->operator()(1, 2) = 0;
		ptrOFF->operator()(2, 1) = 0; ptrOFF->operator()(2, 2) = 0;

		//creating state vector
		state_vector = matrix{ (int)pow(2,number_of_qbits), 1 }; //pow return double type, so i'm casting them to int
		state_vector(1, 1) = 1;
		for (int i{ 2 }; i <= state_vector.get_size(); i++) {
			state_vector(i, 1) = 0;
		}

	}

	matrix get_matrix() { //computes a full matrix representation of the entire circuit

		//start with identity
		matrix matrix_rep{ (int)pow(2, number_of_qbits),(int)pow(2, number_of_qbits) };
		for (int i = 1; i <= matrix_rep.get_rows(); i++) {
			for (int j = 1; j <= matrix_rep.get_cols(); j++) {
				if (i == j) matrix_rep(i, j) = 1;
				else matrix_rep(i, j) = 0;

			}
		}
		// then iteratively compute matrix multiplication
		for (int i{ static_cast<int>(gate_vector.size()) - 1 }; i >= 0; i--) { // needs a static cast as gate_vector.size returns unsigned
			matrix_rep *= gate_vector[i]->get_matrix();
		}
		return matrix_rep;
	}

	void print_gate() {
		for (int i{ 0 }; i < number_of_qbits; i++) {
			std::cout << character << '-';
			wprintf(L"\x1b[1B");
			wprintf(L"\x1b[2D");

		}
	}

	//void save_subcircuit( std::string name, quantum_circuit circuit){
	//	if (circuit.number_of_qbits != number_of_qbits) { std::cout << "size of circuits don't match"; exit(1); }
	//	custom_circuits.push_back(new quantum_circuit{ circuit });
	//}

	//void add_subcircuit(quantum_circuit circuit) {
	//	if (circuit.number_of_qbits != number_of_qbits) { std::cout << "size of circuits don't match"; exit(1); }

	//}

	void set_register(std::string setting) { // this function allows the user to initialize the register to any basis vector
		if (setting.size() != number_of_qbits) { std::cout << "register setting has wrong number of qbits"; exit(1); }
		int state_vector_pos = 1;
		std::reverse(setting.begin(), setting.end()); // flip order so we read from the least significant bit
		for (int i{0}; i < number_of_qbits; i++) { // using a plain for loop instead of iterator as I need to know the index
			if (setting[i] == '1') {
				state_vector_pos += (int)pow(2, i); } //calculating index of given basis vector
			else if (setting[i] != '0') { std::cout << "register setting must be in binary"; exit(1); }
		}

		//creating state vector
		for (int i{ 1 }; i <= state_vector.get_size(); i++) {
			if (i == state_vector_pos) { state_vector(i, 1) = 1; }
			else { state_vector(i, 1) = 0; }
		}
	}

	void show_result() {
		//apply circuit matrix then print results
		matrix output = (get_matrix() * state_vector).get_modulus();
		std::string bit_representation;
		for (int i{ 1 }; i <= output.get_size(); i++) {
			if (output(i, 1).real() > 1.5494e-32) { //When using unitary gates, floating point precission can make zeros become very tiny numbers. We don't want to print those
				bit_representation = std::bitset<64>(i - 1).to_string();
				std::cout << bit_representation.substr(64 - number_of_qbits, -1) << ' ' << output(i, 1) << std::endl;
			}

		}
	}
	void add_parallel() {
		gate_vector.push_back(std::shared_ptr<full_lane_gate>(new parallel_gates{ number_of_qbits, ptrI }));
	}
	void barrier() {

		gate_vector.push_back(std::shared_ptr<full_lane_gate>(new wall{ number_of_qbits })); //makes a gap in the circuit to seperate out steps of a computation
	}
	void add_control(int control_bit, int controlled_bit, const std::shared_ptr<single_bit_gate> & control_pointer) {
		gate_vector.push_back(std::shared_ptr<full_lane_gate>(new control_gate{ number_of_qbits, control_bit, controlled_bit, ptrON, ptrOFF, ptrI, control_pointer }));
	}
	void add_gate(int qbit_number, const std::shared_ptr<single_bit_gate> & gate_pointer) {
		if (gate_vector.size() == 0 || gate_vector.back()->read_gate(qbit_number) != ptrI) { //check if new parallel gate needs to be added
			add_parallel();
		}
		gate_vector.back()->add_gate(qbit_number, gate_pointer);
	}

	//standard gates
	void h(int qbit_number)
	{
		add_gate(qbit_number, ptrH);
	}
	void x(int qbit_number)
	{
		add_gate(qbit_number, ptrX);
	}
	void y(int qbit_number)
	{
		add_gate(qbit_number, ptrY);
	}
	void z(int qbit_number)
	{
		add_gate(qbit_number, ptrZ);
	}
	void p(int qbit_number)
	{
		add_gate(qbit_number, ptrP);
	}

	//standard control gates
	void ch(int control_bit, int controlled_bit) {
		add_control(control_bit, controlled_bit, ptrH);
	}
	void cx(int control_bit, int controlled_bit) {
		add_control(control_bit, controlled_bit, ptrX);
	}
	void cy(int control_bit, int controlled_bit) {
		add_control(control_bit, controlled_bit, ptrY);
	}
	void cz(int control_bit, int controlled_bit) {
		add_control(control_bit, controlled_bit, ptrZ);
	}
	void cp(int control_bit, int controlled_bit) {
		add_control(control_bit, controlled_bit, ptrP);
	}

	//functions for creating custom gates to add multiple times
	void create_unitary(std::string name, double theta, double phi, double lambda, char character = 'U') {
		gate_map.insert({name, custom_gates.size()});
		custom_gates.push_back(std::shared_ptr<u_gate>(new  u_gate{ theta, phi, lambda }));
	}
	void su(std::string name, int qbit_number) {
		add_gate(qbit_number, custom_gates[gate_map[name]]);
	} // add saved unitary
	void scu(std::string name,int control_bit, int controlled_bit) {
		add_control(control_bit, controlled_bit, custom_gates[gate_map[name]]);
	} // add saved unitary as controlled gate

	//functions for adding a single custom gate
	void u(double theta, double phi, double lambda, int qbit_number, char character = 'U') {
		custom_gates.push_back(std::shared_ptr<u_gate>(new  u_gate{ theta, phi, lambda }));
		add_gate(qbit_number, custom_gates.back());
	}
	void cu(double theta, double phi, double lambda, int control_bit, int controlled_bit) {
		custom_gates.push_back(std::shared_ptr<u_gate>(new  u_gate{ theta, phi, lambda }));
		add_control(control_bit, controlled_bit, custom_gates.back());
	}

	void print_circuit() {
		for (int i{ 0 }; i < gate_vector.size(); i++) {
			gate_vector[i]->print_gate();
			printf("\x1b[%dA", number_of_qbits);
			printf("\x1b[2C");

			if ((i + 1) % 30 == 0)
			{
				wprintf(L"\x1b[60D");
				wprintf(L"\x1b[%dB", number_of_qbits + 1);
			}
		}
		wprintf(L"\x1b[%dB", number_of_qbits);
	}
};


#endif