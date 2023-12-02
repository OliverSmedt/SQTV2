
#include "SQTV2.h"

int main()
{
	quantum_circuit qc{ 5 };
	qc.h(1);
	qc.h(2);
	qc.h(3);
	qc.h(4);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);
	qc.h(0);


	qc.barrier();

	//std::cout << qc.get_matrix();

	qc.show_result();
	qc.print_circuit();
}